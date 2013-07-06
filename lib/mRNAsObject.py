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
