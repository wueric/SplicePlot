import numpy
import pysam

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
        return self.chrm is None or self.low is None or self.high is None or self.wiggle is None or self.junctions_dict is None
        
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
