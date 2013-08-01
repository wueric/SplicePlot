import pysam
import argparse

def get_bam_locations(map_file):
    f1 = open(map_file,'r').readlines()
    bam_list = []
    for line in f1:
        bam_path = line.strip('\n').split()[1]
        bam_list.append(bam_path)
    return bam_list

def get_alternative_junctions(window,bam_list):
    
    alternative_junctions_list = []

    broken_window = window.split(':')
    contig_name = broken_window[0]
    low_coordinate = int(broken_window[1].split('-')[0])-1
    high_coordinate = int(broken_window[1].split('-')[1])+1

    junctions_by_low_ss = {}
    junctions_by_high_ss = {}
    
    for bam_file_name in bam_list:
        try:
            bam_file = pysam.Samfile(bam_file_name,'rb')
            relevant_reads = bam_file.fetch(contig_name,low_coordinate,high_coordinate)
            
            for read in relevant_reads:
                for index, base_position in enumerate(read.positions):
                    if (index+1) < len(read.positions) and base_position + 1 != read.positions[index+1]:
                        junction_low_coordinate = base_position + 1
                        junction_high_coordinate = read.positions[index+1]+1
                        
                        if junction_low_coordinate not in junctions_by_low_ss:
                            junctions_by_low_ss[junction_low_coordinate] = set()
                        if junction_high_coordinate not in junctions_by_high_ss:
                            junctions_by_high_ss[junction_high_coordinate] = set()

                        junctions_by_low_ss[junction_low_coordinate].add(junction_high_coordinate)
                        junctions_by_high_ss[junction_high_coordinate].add(junction_low_coordinate)
        except IOError:
            pass

    for key, value in junctions_by_low_ss.items():
        if len(value) > 1:
            upper_coordinates = sorted(list(value))
            junction_string = ''
            for upper_coordinate in upper_coordinates:
                junction_string += '{0}:{1}-{2},'.format(contig_name,key,upper_coordinate)

            alternative_junctions_list.append(junction_string[0:len(junction_string)-1])

    for key, value in junctions_by_high_ss.items():
        if len(value) > 1:
            lower_coordinates = sorted(list(value))
            junction_string = ''
            for lower_coordinate in lower_coordinates:
                junction_string += '{0}:{1}-{2},'.format(contig_name,lower_coordinate,key)

            alternative_junctions_list.append(junction_string[0:len(junction_string)-1])

    return alternative_junctions_list

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get potential alternative junctions in a genomic coordinate window')

    parser.add_argument('window',type=str,help='Genomic window, in the format chr1:10000-20000')
    parser.add_argument('mf',type=str,help='map file containing locations of .bam files')

    args = parser.parse_args()
    
    bam_list = get_bam_locations(args.mf)

    alternative_junctions_list = get_alternative_junctions(args.window,bam_list)
    for junction in alternative_junctions_list:
        print junction
