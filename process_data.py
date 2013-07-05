import pysam
import numpy
import pandas
import cPickle as pickle


class ReadDepth:

    def __init__(self,chrm,low,high,wiggle,junctions_dict):
        assert chrm == None or high - low + 1 == len(wiggle), 'Wiggle, lower bound, and upper bound do not correspond'
        self.low = low
        self.high = high
        self.chrm = chrm
        self.wiggle = wiggle
        self.junctions_dict = junctions_dict

    @classmethod
    def determine_depth(cls,bam_file_path,chrm,start_coord,end_coord):

        '''
            determine_depth determines the coverage at each base between start_coord and end_coord, inclusive.

            bam_file_path is the path to the bam file used to determine the depth and junctions on chrm between start_coord and end_coord

            return values:
                depth_vector, which is a Numpy array which contains the coverage at each base position between start_coord and end_coord
                spanned_junctions, which is a dictionary containing the junctions supported by reads. The keys in spanned_junctions are the
                    names of the junctions, with the format chromosome:lowerBasePosition-higherBasePosition
        '''
        try:
            bam_file = pysam.Samfile(bam_file_path,'rb')
            relevant_reads = bam_file.fetch(reference=chrm,start=start_coord,end=end_coord)

            depth_vector = numpy.zeros(end_coord-start_coord+1,dtype='f')
            spanned_junctions = {}

            for read in relevant_reads:
                
                # make sure that the read can be used
                cigar_string = read.cigar

                # each read must have a cigar string
                if cigar_string == None:
                    continue

                # read cannot have insertions or deletions
                contains_indel = False
                spans_more_than_one_junction = False
                for cigar_event in cigar_string:
                    if cigar_event[0] == 1 or cigar_event[0] == 2:
                        contains_indel = True
                        break
                    
                if contains_indel:
                    continue
                    
                for index, base_position in enumerate(read.positions):
                    if base_position >= start_coord and base_position <= end_coord:
                        depth_vector[base_position-start_coord] += 1

                    # junction spanning case
                    if (index+1) < len(read.positions) and base_position + 1 != read.positions[index+1]:
                        junction_name = '{0}:{1}-{2}'.format(chrm,base_position+1,read.positions[index+1]+1)
                        if junction_name not in spanned_junctions:
                            spanned_junctions[junction_name] = 0

                        spanned_junctions[junction_name] = spanned_junctions[junction_name] + 1

            return cls(chrm,start_coord,end_coord,depth_vector,spanned_junctions)
        except IOError:
            print 'There is no .bam file at {0}'.format(bam_file_path)
            raise Exception

    @classmethod
    def create_blank(cls):

        '''
            create_blank creates an instance of ReadDepth where all of the attributes are None
        '''
        return cls(None,None,None,None,None)

    def is_invalid(self):
        '''
            is_invalid determines whether any of the attributes are None
        '''
        return self.chrm == None or self.low == None or self.high == None or self.wiggle == None or self.junctions_dict == None
        
    def shrink(self,new_low,new_high):

        '''
            shrink changes the boundaries of the ReadDepth object

            new_low is the new lower genomic coordinate boundary for the ReadDepth object
            new_high is the new upper genomic coordinate boundary for the ReadDepth object

            This method also changes self.wiggle and self.junctions_dict so that they only contain data between new_low and new_high

            return value:
                Nothing. Method changes the ReadDepth object
        '''
        if new_low < self.low or new_high > self.high:
            raise Exception, 'New boundaries are not valid'

        # filter through junctions_dict to remove junctions which are no longer in the region
        new_junctions_dict = {}
        for key, value in self.junctions_dict.items():
            ss_low, ss_high = map(int,key.split(':')[1].split('-'))
            
            if ss_low >= new_low and ss_high <= new_high:
                new_junctions_dict[key] = value

        self.junctions_dict = new_junctions_dict

        # replace the wiggle
        bottom_index = new_low - self.low
        top_index = bottom_index + (new_high - new_low)

        self.wiggle = self.wiggle[bottom_index:top_index+1]

        # change the lower and upper bound (last step)
        self.low = new_low
        self.high = new_high

    def __add__(self,other):

        '''
            __add__ allows two ReadDepth objects to be added together using the + symbol

            Both self and other must have the same low and high attributes

            return value:
                A new ReadDepth object containing the sum of the two original ReadDepth objects
        '''

        if self.is_invalid():
            return other
        if other.is_invalid():
            return self

        assert self.chrm == other.chrm, 'Cannot add depths from different chromosomes'
        assert self.low == other.low and self.high == other.high, 'Cannot add depths with different start and end points'
        new_wiggle = self.wiggle + other.wiggle

        new_junctions_dict = {}
        
        for key, value in self.junctions_dict.items():
            if key in other.junctions_dict:
                new_junctions_dict[key] = value + other.junctions_dict[key]
            else:
                new_junctions_dict[key] = value

        for key, value in other.junctions_dict.items():
            if key not in self.junctions_dict:
                new_junctions_dict[key] = value

        return ReadDepth(self.chrm,self.low,self.high,new_wiggle,new_junctions_dict)

    def __str__(self):
        return '{0}:{1}-{2},{3},{4}'.format(self.chrm,self.low,self.high,self.wiggle,self.junctions_dict)

    def divide_by_constant(self,constant):

        '''
            divide_by_constant divides self.wiggle and self.junctions_dict by a constant value

            constant is a number

            return value:
                A new ReadDepth object containing the divided values. Method leaves the original ReadDepth object unchanged
        '''
        new_wiggle = self.wiggle / constant

        new_junctions_dict = {}
        for key, value in self.junctions_dict.items():
            new_junctions_dict[key] = value * 1.0 / constant

        return ReadDepth(self.chrm,self.low,self.high,new_wiggle,new_junctions_dict)

    def filter_junctions_dict_for_event(self,splice_event_name):
        '''
            filter_junctions_dict_for_event removes all entries frm junctions_dict that cannot possibly be
                involved in the alternative splicing event splice_event_name

            splice_event_name is the name of the alternative splicing event, in the format
                format chr1:17055-17915,chr1:17055-17606,chr1:17055-17233,
                where the numbers represent the genomic coordinates of the splice sites

            return values:
                A new ReadDepth object containing only the relevant junctions
        '''
        junction_names_list = splice_event_name.split(',')
        new_junctions_dict = {}
        for junction_name in junction_names_list:
            if junction_name in self.junctions_dict:
                new_junctions_dict[junction_name] = self.junctions_dict[junction_name]

        return ReadDepth(self.chrm,self.low,self.high,self.wiggle,new_junctions_dict)



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
    shared_site = int(junctions_list[0].split(':')[1].split('-')[0])
    splice_junc_coordinate_list = map(lambda x: int(x.split(':')[1].split('-')[1]),junctions_list)
    splice_junc_coordinate_list.append(shared_site)
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
        prop_covered = exon.determine_proportion_covered(total_depth)
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


class mRNAsObject:

    '''
        mRNAsObject represents a set of possible mRNA segments

        chrm is the chromosome number
        strand is the strand
        low is the lowest possible coordinate in the set of possible mRNAs
        high is the highest possible coordinate in the set of possible mRNAs
        mRNAs is a list containing possible mRNAs. Each possible mRNA is represented by a list of exons, and
            each exon is represented by a list containing its lowest coordinate and its highest coordinate
    '''

    def __init__(self,chrm,strand,low,high,mRNAs):
        self.chrm = chrm
        self.strand = strand
        self.low = low
        self.high = high
        self.mRNAs = mRNAs


        exon_starts = []
        exon_ends = []

        for possible_mRNA in mRNAs:
            for exon_coordinate_list in possible_mRNA:
                exon_starts.append(exon_coordinate_list[0])
                exon_ends.append(exon_coordinate_list[1])

        self.exon_starts = sorted(list(set(exon_starts)))
        self.exon_ends = sorted(list(set(exon_ends)))

    def __str__(self):
        return '{0}:{1}-{2}{3},{4}'.format(self.chrm,self.low,self.high,self.strand,self.mRNAs)


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
    
    
# debugging code
if __name__ == '__main__':
    
    bam_to_id_dict, bam_list = map_indiv_id_to_bam_name('/home/wueric/bam/map_file.txt')
    read_depths_dict, mRNA_info_object = initialize_read_depths_and_determine_exons('chr1:17055-17915,chr1:17055-17606,chr1:17055-17233','/home/wueric/no_duplicates.unzip.exons.gencode.v17.annotation.gtf.gz',bam_list,bam_to_id_dict)
    print read_depths_dict
    print mRNA_info_object

    
    genotype_averages_dict, genotype_by_id = average_read_depth_by_genotype(read_depths_dict,'/home/wueric/bam/sample.vcf.gz','chr1:10583')

    print genotype_averages_dict

    print genotype_by_id

    df = create_data_frame(read_depths_dict,'chr1:17055-17915,chr1:17055-17606,chr1:17055-17233','chr1:10583',genotype_by_id)
    print df.to_string()

    pickle.dump(read_depths_dict,open('read_depths_dict.p','wb'))
    pickle.dump(mRNA_info_object,open('possible_mRNAs.p','wb'))
    pickle.dump(df,open('df_pickled.p','wb'))
    pickle.dump(genotype_averages_dict,open('genotype_averages_pickled.p','wb'))
