import pysam
import pandas
import cPickle as pickle
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
import argparse
import os


class VCFLine:
    def __init__(self,VCF_name,region):
        '''
            VCF_name is the name of the gzipped vcf file
            region is the position of the SNP, in format chr4:12345
        '''
        VCF_object = pysam.Tabixfile(VCF_name)

        VCF_header = list(VCF_object.header)
        last_header_line = VCF_header[len(VCF_header)-1]
        self.samples = last_header_line.split()[9:]

        region_list = region.split(':')

        self.contig = region_list[0]
        self.position = int(region_list[1])
        self.id = None
        self.ref = None
        self.alt = None
        self.genotype_calls = {}
        
        try:
            feature_iterator = VCF_object.fetch(region)
            for feature in feature_iterator:
                vcf_line_array = feature.strip('\n').split()
                
                contig_name = vcf_line_array[0]
                position = int(vcf_line_array[1])
                if contig_name == self.contig and position == self.position:
                    self.id = vcf_line_array[2]
                    self.ref = vcf_line_array[3]
                    self.alt = filter(lambda x: x != '.', vcf_line_array[4].split(','))
                    
                    genotype_calls_list = vcf_line_array[9:]

                    for i, indiv_id in enumerate(self.samples):
                        if '.' not in genotype_calls_list[i].split(':')[0]:
                            self.genotype_calls[indiv_id] = genotype_calls_list[i]
                    break

            if self.id == None:
                print "There is no variant at {0}".format(region)
                raise Exception

        except ValueError:
            print "{0} is not a valid SNP position for this VCF file".format(region)
            raise Exception
                

    def __getitem__(self,key):
        try:
            return self._determine_genotype_bases(self.genotype_calls[key])
        except KeyError:
            return None

    def __contains__(self,key):
        return key in self.genotype_calls

    def alleles_list(self):
        alleles = [self.ref]
        alleles.extend(self.alt)
        return alleles

    def genotypes_list(self):
        genotypes_list = []
        alleles_list = self.alleles_list()
        for i in range(len(alleles_list)):
            for j in range(i,len(alleles_list)):
                genotypes_list.append('{0}{1}'.format(alleles_list[i],alleles_list[j]))
        return genotypes_list

    def _determine_genotype_bases(self,genotype_string):
        calls = genotype_string.split(':')[0]
        calls_as_numeric = None
        if '|' in calls:
            calls_as_numeric = sorted(map(int,calls.split('|')))
        else:
            calls_as_numeric = sorted(map(int,calls.split('/')))

        genotype_string = ''.join(map(lambda x: self.alleles_list()[x],calls_as_numeric))
        return genotype_string
        
        

def average_read_depth_by_genotype(read_depth_dict,vcf_file_name,var_pos):
    '''
        average_read_depth_by_genotype averages the ReadDepth objects in read_depth_dict based on the genotype
            in a .vcf file

        read_depth_dict is a dict of ReadDepth objects, where the keys are the individual IDs and the values are
            the corresponding ReadDepth objects

        vcf_file_name is the location of the .vcf file

        var_pos is the position of the SNP, with format chr1:12345

        return values:
            A dict containing the average read depths. The keys are the genotypes, and the values are the
                corresponding ReadDepth objects
            
            A dict which maps indiv IDs to genotype. The keys are the indiv IDs, and the values are the
                corresponding genotypes

            A list of genotypes which occur in the dataset, where homozygous reference is first

    '''

    try:
        vcf_line = VCFLine(vcf_file_name,var_pos)

        possible_genotypes_bucket_counts = {}
        average_read_depths_dict = {}
        genotypes_in_data = set()
        genotype_by_id = {}
        for indiv_id, read_depth_object in read_depth_dict.items():
            if indiv_id in vcf_line:
                indiv_genotype = vcf_line[indiv_id]
                genotype_by_id[indiv_id] = indiv_genotype
                genotypes_in_data.add(indiv_genotype)

                if indiv_genotype not in possible_genotypes_bucket_counts:
                    possible_genotypes_bucket_counts[indiv_genotype] = 1
                    average_read_depths_dict[indiv_genotype] = read_depth_object
                else:
                    possible_genotypes_bucket_counts[indiv_genotype] += 1
                    average_read_depths_dict[indiv_genotype] = average_read_depths_dict[indiv_genotype] + read_depth_object

        for genotype, counts in possible_genotypes_bucket_counts.items():
            average_read_depths_dict[genotype].divide_by_constant(counts)

        filtered_genotypes_list = filter(lambda x: x in genotypes_in_data, vcf_line.genotypes_list())
        return average_read_depths_dict, genotype_by_id, filtered_genotypes_list

    except IOError:
        print 'There is no .vcf file at {0}'.format(vcf_file_name)
        raise Exception


def map_indiv_id_to_bam_name(id_map_file):

    '''
        map_indiv_id_to_bam_name maps the individual ids in the vcf file to the corresponding .bam files

        id_map_file is the name of the file containing the id to bam file correspondences

        return value:
            a dict that maps individual id to .bam file location

            a list of all the .bam file locations to process
    '''

    bam_to_id = {}
    bam_list = []

    try:
        lines = open(id_map_file,'r').readlines()
        for line in lines:
            line = line.strip('\n').split()

            file_path = line[1]
            if os.path.exists(file_path):
                bam_to_id[line[1]] = line[0]
                bam_list.append(line[1])
            else:
                print 'Skipping {0}. Invalid file path'.format(file_path)

        return bam_to_id, bam_list

    except IOError:
        print 'There is no mapping file at {0}'.format(id_map_file)
        raise Exception


class Exon:
    def __init__(self,chrm,low,high,strand):
        self.chrm = chrm
        self.low = low
        self.high = high
        self.strand = strand

    @classmethod
    def create_from_gtf(cls,line):
        '''
            create_from_gtf creates an Exon object from a single line of a gtf file

            line is a single line of text (a string) from a gtf file
        '''
        info = line.strip('\n').split()
        return cls(info[0],int(info[3]),int(info[4]),info[6])

    def determine_proportion_covered(self,read_depth):
        if read_depth.chrm != self.chrm or self.low < read_depth.low or self.high > read_depth.high:
            return 0
        
        # determine the top and bottom indices
        bottom_index = self.low - read_depth.low
        top_index = bottom_index + (self.high - self.low)

        bases_covered = 0
        for item in range(bottom_index, top_index + 1):
            if read_depth.wiggle[item] > 0:
                bases_covered += 1.0

        return bases_covered / (top_index + 1 - bottom_index)

    def determine_average_coverage(self,read_depth):
        return self.determine_proportion_covered(read_depth) / (self.high - self.low + 1.0)

    def __str__(self):
        return '{0}:{1}-{2}{3}'.format(self.chrm,self.low,self.high,self.strand)

class EvaluatedExon:
    def __init__(self,exon,value):
        self.exon = exon
        self.value = value
        self.length = exon.high - exon.low + 1

    def __cmp__(self,other):
        if self.value - other.value != 0:
            return self.value - other.value

        return self.length - other.length

    def __str__(self):
        return '{0},{1}'.format(self.exon.__str__(), self.value)

def get_splice_junc_coordinates(junction_name):

    '''
        get_splice_junc_coordinates takes in a junction name and returns the chromosome name,
            the shared splicing site, and the other splicing sites

        junction_name is a string containing the junction name, in the format "chr1:100-200,chr1:100-300"

        return values:
            chrom: a string containing the chromosome name
            shared_site: an int containing the base number of the shared splice site
            upper_ss: a list of ints containing the base numbers of the other splice sites
    '''

    junctions_list = junction_name.split(',')
    
    chromosome_name_list = map(lambda x: x.split(':')[0], junctions_list)
    if not all(x==chromosome_name_list[0] for x in chromosome_name_list):
        print '{0} is not a valid junction_name'.format(junction_name)
        raise Exception

    chrom = chromosome_name_list[0]

    # figure out the shared splice site
    splice_junc_coordinate_list = []
    shared_site = None

    lower_ss, upper_ss = [], []
    for intronic_region in junctions_list:
        low, high = map(int, intronic_region.split(':')[1].split('-'))
        if low not in lower_ss:
            lower_ss.append(low)
        if high not in upper_ss:
            upper_ss.append(high)

    if len(upper_ss) == 1:
        shared_site = upper_ss[0]
        return chrom, shared_site, lower_ss

    elif len(lower_ss) == 1:
        shared_site = lower_ss[0]
        return chrom, shared_site, upper_ss
    else:
        print '{0} is not a valid junction name'.format(junction_name)
        raise Exception


def determine_exons_and_coordinates(gtf_file_name,chrom_name,shared_site,other_sites):

    '''
        determine_exons_and_coordinates creates a list of possible exons that could have been involved
            the creating the junction, and gives the lower and upper coordinates of the region to search
            for in the bam file

        gtf_file_name is the file path for the gtf file
        chrom_name is a string containing the chromosome name
        shared_site is an integer representing the base number of the shared splice site
        other_sites is a list of integers containing the base numbers of the other splice sites


        return values:
            min_coordinate is an int representing the lower border of the region to search in the bam files
            max_coordinate is an int representing the upper border of the region to search in the bam files
            filtered_exons_list is a list of Exon objects containing exons which could have been involved in
                the splice junction

    '''

    exons_found_buckets = {shared_site:False}
    for coordinate in other_sites:
        exons_found_buckets[coordinate] = False

    splice_junc_coordinate_list = [shared_site]
    splice_junc_coordinate_list.extend(other_sites)

    try:
        gtf_tabix = pysam.Tabixfile(gtf_file_name,'r')
        relevant_exons_iterator = gtf_tabix.fetch(chrom_name,min(splice_junc_coordinate_list)-1,max(splice_junc_coordinate_list)+1)

        filtered_exons_list = []
        min_coordinate = float('inf')
        max_coordinate = float('-inf')
        for exon_line in relevant_exons_iterator:
            # shared_site coordinate < other_sites coordinates case
            exon = Exon.create_from_gtf(exon_line)
            if shared_site > other_sites[0]:
                if exon.low == shared_site:
                    filtered_exons_list.append(exon)
                    exons_found_buckets[shared_site] = True

                    if exon.high > max_coordinate:
                        max_coordinate = exon.high

                elif exon.high in other_sites:
                    filtered_exons_list.append(exon)
                    exons_found_buckets[exon.high] = True
                
                    if exon.low < min_coordinate:
                        min_coordinate = exon.low
            else:
                if exon.high == shared_site:
                    filtered_exons_list.append(exon)
                    exons_found_buckets[shared_site] = True

                    if exon.low < min_coordinate:
                        min_coordinate = exon.low
                elif exon.low in other_sites:
                    filtered_exons_list.append(exon)
                    exons_found_buckets[exon.low] = True
                    
                    if exon.high > max_coordinate:
                        max_coordinate = exon.high

        if min_coordinate == float('inf') or max_coordinate == float('-inf') or not all(exons_found_buckets[x] for x in exons_found_buckets):
            print 'The given junction coordinates do not correspond to exons in the annotation'
            raise Exception

        return min_coordinate, max_coordinate, filtered_exons_list
    except IOError:
        print 'There is no gtf file at {0}'.format(gtf_file_name)
        raise Exception
                


def initialize_read_depths_and_determine_exons(junction_name,gtf_file_name,bam_list,bam_to_id_dict):

    '''
        initialize_read_depths_and_determine_exons creates a dictionary of read depths and assembles a possible set of exons

        junction_name is the name of the junction being examined. It is a string with the 
            format chr1:17055-17915,chr1:17055-17606,chr1:17055-17233,
            where the numbers represent the genomic coordinates of the splice sites

        gtf_file_name is a string representing the path to the gtf file containing known exons

        bam_list is a list of strings containing the file paths to the bam files

        return values:
            A dictionary of ReadDepth objects, where the key is the file path (as a string) to the bam file,
                and where the value is a ReadDepth object corresponding to the bam file

            A mRNAsObject, which represents the set of possible mRNA segments determined from the junctions
    '''

    chrom, shared_site, other_sites = get_splice_junc_coordinates(junction_name)
    minimum_coordinate, maximum_coordinate, filtered_exons_list = determine_exons_and_coordinates(gtf_file_name,chrom,shared_site,other_sites)

    splice_junc_coordinate_list = [shared_site]
    splice_junc_coordinate_list.extend(other_sites)

    read_depth_dict = {}
    total_depth = ReadDepth.create_blank()
    for bam_file in bam_list:
        current_read_depth = ReadDepth.determine_depth(bam_file,chrom,minimum_coordinate,maximum_coordinate)

        read_depth_dict[bam_to_id_dict[bam_file]] = current_read_depth
        total_depth = total_depth + current_read_depth

    current_best_exon_by_coordinate = {}

    for exon in filtered_exons_list:
        prop_covered = exon.determine_average_coverage(total_depth)
        current_exon_evaluated = EvaluatedExon(exon,prop_covered)

        if exon.low in splice_junc_coordinate_list:
            if exon.low not in current_best_exon_by_coordinate or current_exon_evaluated > current_best_exon_by_coordinate[exon.low]:
                current_best_exon_by_coordinate[exon.low] = current_exon_evaluated
        else:
            if exon.high not in current_best_exon_by_coordinate or current_exon_evaluated > current_best_exon_by_coordinate[exon.high]:
                current_best_exon_by_coordinate[exon.high] = current_exon_evaluated

    # convert current_best_exon_by_coordinate into list format needed by plotter
    mRNAs = []

    # determine which exon is in all transcripts. Also determine boundaries for resizing the read depth objects to fit the exons
    shared_exon = current_best_exon_by_coordinate[shared_site].exon

    resize_lower_bound = float('inf')
    resize_upper_bound = float('-inf')

    for key in sorted(current_best_exon_by_coordinate.keys()):
        exon = current_best_exon_by_coordinate[key].exon

        if exon != shared_exon:
            mRNA_segment = sorted([sorted([shared_exon.low,shared_exon.high]),sorted([exon.low,exon.high])],
                                    cmp=lambda x, y: x[0] - y[0])
            mRNAs.append(mRNA_segment)

        if exon.low < resize_lower_bound:
            resize_lower_bound = exon.low
        if exon.high > resize_upper_bound:
            resize_upper_bound = exon.high


    # resize each of the read depth objects so that their lengths correspond to the possible mRNAs
    for key in read_depth_dict:
        read_depth_dict[key].shrink(resize_lower_bound,resize_upper_bound)

    mRNAs_info = mRNAsObject(exon.chrm,exon.strand,resize_lower_bound,resize_upper_bound,mRNAs)

    return read_depth_dict, mRNAs_info



def create_data_frame(read_depth_dict,junction_name,var_pos,genotype_lookup_dict,filtered_genotypes_list):

    '''
        create_data_frame creates a pandas.DataFrame object in the format required by the hive
            and structure plotting functions

        read_depth_dict is a dict, where the keys are the bam file paths and the values are ReadDepth objects

        junction_name is the name of the junction, in the usual format

        genotype_lookup_dict is a dict, where the key is the indiv ID (as specified in .vcf file)
            and where the value is the individual's genotype

        var_pos is the position of the SNP, in the format chr1:12345

        return values:
            A pandas.DataFrame. Each row represents an individual, and the row name is the indiv id.
                The first column contains the genotypes of each individual. All remaining columns contain
                the splicing ratios corresponding to a splice junction

    '''

    precursor_dict = {}
    data_frame_index = []
    junctions_list = junction_name.split(',')
    for key, value in read_depth_dict.items():

        indiv_genotype = genotype_lookup_dict[key]

        data_frame_index.append(key)

        # add the genotype to precursor_dict
        if var_pos not in precursor_dict:
            precursor_dict[var_pos] = []
        precursor_dict[var_pos].append(indiv_genotype)

        # calculate the splicing ratios
        total_junc_read_count = 0
        for junction in junctions_list:
            if junction in value.junctions_dict:
                total_junc_read_count += value.junctions_dict[junction]

        # add splicing ratios to precursor_dict
        for junction in junctions_list:
            write_value = 0.0

            if junction in value.junctions_dict:
                write_value = value.junctions_dict[junction] * 1.0 / total_junc_read_count

            if junction not in precursor_dict:
                precursor_dict[junction] = []

            precursor_dict[junction].append(write_value)

    # rearrange the data frame so that genotype is the first column
    df = pandas.DataFrame(precursor_dict,index=data_frame_index)
    new_col_order = [var_pos]
    new_col_order.extend(junctions_list)
    return df.reindex(columns=new_col_order)

def calculate_average_expression_and_data_frame(var_pos,junction_name,vcf,annotation,map_file):
    bam_to_id_dict, bam_list = map_indiv_id_to_bam_name(map_file)
    
    read_depths_dict, mRNA_info_object = initialize_read_depths_and_determine_exons(junction_name,
        annotation,bam_list,bam_to_id_dict)

    new_read_depths_dict = {}
    for indiv_id, read_depth_object in read_depths_dict.items():
        new_read_depths_dict[indiv_id] = read_depth_object.filter_junctions_dict_for_event(junction_name)

    genotype_averages_dict, genotype_by_id, filtered_genotypes_list = average_read_depth_by_genotype(new_read_depths_dict,vcf,var_pos)


    data_frame = create_data_frame(new_read_depths_dict,junction_name,var_pos,genotype_by_id,filtered_genotypes_list)

    return genotype_averages_dict, data_frame, mRNA_info_object, filtered_genotypes_list


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Initialize and pickle alternative splice junction data')

    parser.add_argument('varpos',type=str,help='string describing position of SNP. should have format chr_name:base_number')
    parser.add_argument('junc',type=str,help="string representing junction. should have format chr_name:lower_base-upper_base,chr_name:lower_base-upper_base, where lower_base and upper_base represent possible intronic regions")
    parser.add_argument('--vcf',type=str,required=True,help='location of the vcf file')
    parser.add_argument('--gtf',type=str,required=True,help='location of the gtf file')
    parser.add_argument('--mf',type=str,required=True,help='location of the map file')
    parser.add_argument('--output',type=str,required=False,default=None,help='location of output pickle file. Optional parameter')

    args = parser.parse_args()
    try:
        genotype_averages_dict, data_frame, mRNA_info_object, filtered_genotypes_list = calculate_average_expression_and_data_frame(args.varpos,args.junc,args.vcf,args.gtf,args.mf)

        output_file_path = '{0}/pickle_files/'.format(os.path.dirname(os.path.abspath(__file__)))
        if args.output is not None:
            output_file_path = args.output

        stem = output_file_path
        tail = '{0}@{1}.p'.format(args.varpos,args.junc)
        
        if output_file_path[len(output_file_path)-2:] == '.p':
            stem, tail = os.path.split(output_file_path)

        try:
            os.makedirs(stem)
        except OSError:
            if os.path.isdir(stem):
                pass
            else:
                print 'Cannot create directory {0}'.format(stem)
                raise Exception

        if stem != '' and stem[len(stem)-1] != '/':
            stem = stem + '/'

        # pickle the data
        pickle_file = open('{0}{1}'.format(stem,tail),'wb')
        pickle.dump(args.varpos,pickle_file)
        pickle.dump(args.junc,pickle_file)
        pickle.dump(genotype_averages_dict,pickle_file)
        pickle.dump(mRNA_info_object,pickle_file)
        pickle.dump(data_frame,pickle_file)
        pickle.dump(filtered_genotypes_list,pickle_file)

        pickle_file.close()
    except:
        print 'Failed'
