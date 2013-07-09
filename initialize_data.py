import pysam
import pandas
import cPickle as pickle
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
import argparse
import os


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

    '''

    chrm = var_pos.split(':')[0]
    position = int(var_pos.split(':')[1])

    vcf_object = pysam.VCF()

    try:
        vcf_object.connect(vcf_file_name)

        
        variant_line = list(vcf_object.fetch('{0}:{1}-{1}'.format(chrm,position)))
        individual_id_list = vcf_object.getsamples()

        # the outer for loop should run at most once
        variant = None
        for variant in variant_line:

            alleles_list = list(variant.ref)
            alt_alleles = variant.alt

            alleles_list.extend(alt_alleles)

            read_depths_by_genotype = {}
            genotype_by_id = {}

            for indiv_id, read_depth_object in read_depth_dict.items():
                indiv_genotype_string = variant[indiv_id]['GT'][0]


                indiv_genotype_list_numeric = map(int,[indiv_genotype_string[0], indiv_genotype_string[2]])
                
                genotype_list_bases = []
                for item in indiv_genotype_list_numeric:
                    genotype_list_bases.append(alleles_list[item])
                genotype_list_bases.sort()
            
                genotype_key = ''.join(genotype_list_bases)

                if genotype_key not in read_depths_by_genotype:
                    read_depths_by_genotype[genotype_key] = []

                read_depths_by_genotype[genotype_key].append(read_depth_object)

                genotype_by_id[indiv_id] = genotype_key

            # compute the averages for each genotype
            average_read_depths_dict = {}
            for genotype, depth_list in read_depths_by_genotype.items():
                total_depth = ReadDepth.create_blank()
                for read_depth in depth_list:
                    total_depth += read_depth
                average_depth = total_depth.divide_by_constant(len(depth_list))
                average_read_depths_dict[genotype] = average_depth

            return average_read_depths_dict, genotype_by_id

        if variant == None:
            print 'There is no variant at {0}'.format(var_pos)
            raise Exception

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
    id_to_bam = {}
    bam_list = []

    try:
        lines = open(id_map_file,'r').readlines()
        for line in lines:
            line = line.strip('\n').split()
            bam_to_id[line[1]] = line[0]
            bam_list.append(line[1])

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

    junctions_list = junction_name.split(',')
    
    chrom = junctions_list[0].split(':')[0]

    # figure out the shared site

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
        splice_junc_coordinate_list.extend(lower_ss)
        splice_junc_coordinate_list.extend(upper_ss)
        shared_site = upper_ss[0]
    elif len(lower_ss) == 1:
        splice_junc_coordinate_list.extend(upper_ss)
        splice_junc_coordinate_list.extend(lower_ss)
        shared_site = lower_ss[0]

    else:
        print '{0} is not a valid junction name'.format(junction_name)
        raise Exception
            
    splice_junc_coordinate_list = set(splice_junc_coordinate_list)

    gtf_tabix = pysam.Tabixfile(gtf_file_name,'r')
    relevant_exons_iterator = gtf_tabix.fetch(chrom,min(splice_junc_coordinate_list)-1,max(splice_junc_coordinate_list)+1)

    filtered_exons_list = []
    minimum_coordinate = float('inf')
    maximum_coordinate = float('-inf')
    for exon_line in relevant_exons_iterator:
        exon = Exon.create_from_gtf(exon_line)
        if exon.low in splice_junc_coordinate_list or exon.high in splice_junc_coordinate_list:
            filtered_exons_list.append(exon)
            
            if exon.low < minimum_coordinate:
                minimum_coordinate = exon.low
            if exon.high > maximum_coordinate:
                maximum_coordinate = exon.high

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



def create_data_frame(read_depth_dict,junction_name,var_pos,genotype_lookup_dict):

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

    genotype_averages_dict, genotype_by_id = average_read_depth_by_genotype(new_read_depths_dict,vcf,var_pos)


    data_frame = create_data_frame(new_read_depths_dict,junction_name,var_pos,genotype_by_id)

    return genotype_averages_dict, data_frame, mRNA_info_object


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Initialize and pickle alternative splice junction data')

    parser.add_argument('--varpos',type=str,required=True,help='string describing position of SNP. should have format chr_name:base_number')
    parser.add_argument('--junc',type=str,required=True,help="string representing junction. should have format chr_name:lower_base-upper_base,chr_name:lower_base-upper_base, where lower_base and upper_base represent possible intronic regions")
    parser.add_argument('--vcf',type=str,required=True,help='location of the vcf file')
    parser.add_argument('--gtf',type=str,required=True,help='location of the gtf file')
    parser.add_argument('--mf',type=str,required=True,help='location of the map file')
    parser.add_argument('--output',type=str,required=False,default=None,help='location of output pickle file. Optional parameter')

    args = parser.parse_args()

    genotype_averages_dict, data_frame, mRNA_info_object = calculate_average_expression_and_data_frame(args.varpos,args.junc,args.vcf,args.gtf,args.mf)

    output_file_path = args.output
    if args.output == None:
        output_file_path = '{0}/pickle_files/{1}@{2}.p'.format(os.path.dirname(os.path.abspath(__file__)),args.varpos,args.junc)

    # pickle the data
    pickle_file = open(output_file_path,'wb')
    pickle.dump(args.varpos,pickle_file)
    pickle.dump(args.junc,pickle_file)
    pickle.dump(genotype_averages_dict,pickle_file)
    pickle.dump(mRNA_info_object,pickle_file)
    pickle.dump(data_frame,pickle_file)

    pickle_file.close()
