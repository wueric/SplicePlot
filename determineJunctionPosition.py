import sys
import optparse
import subprocess
import os
def processBedFiles(sourcePath,destPath):

    # open the files in sourcePath
    listOfBedFilesToOpen = os.listdir(sourcePath)

    for name in listOfBedFilesToOpen:

        print 'processing ', name
        # get the actual path of the file
        fileName = sourcePath + '/' + name
        newFileName = destPath +'/' + name

        # open the bed file and process
        # skip the first line because the first line has no data
        # read from f1
        # write to f2
        f1 = open(fileName,'r')
        f2 = open(newFileName,'w')
        
        next(f1)
        for line in f1:

            line = line.strip('\n').split()

            chrm = line[0]
            readStart = int(line[1])
            readEnd = int(line[2])
            junctionName = line[3]
            readDepth = line[4]
	    strand = line[5]
            segmentLengths = line[10]
            firstSegmentLength = int(segmentLengths.split(',')[0])
            secondSegmentLength = int(segmentLengths.split(',')[1])
            

            intronStart = readStart + firstSegmentLength
            intronEnd = readEnd - secondSegmentLength + 1

            f2.write('%s\t%s\t%s\t%s\t%s\t%s' % (str(chrm),str(intronStart),str(intronEnd),junctionName,str(readDepth),strand))
            f2.write('\n')

        f2.close()
        f1.close()

        


def main_modify_bed(argv=None):

        if not argv: argv = sys.argv

        usage = "usage: junctionstk --junctions_enrichment [options]"

        parser = optparse.OptionParser(usage=usage)
        parser.add_option("--sloc", type="string", dest = "sloc", help="folder where source bed files are")
        parser.add_option("--wloc", type="string", dest = "wloc", help="folder where processed bed files are written")
        parser.set_defaults(sloc=None,wloc=None)
        (options, args)=parser.parse_args(argv[1:])

        if options.sloc==None:
                print "no source of bed files given"
                sys.exit(1)

        if options.wloc==None:
                print "no destination of bed files given"
                sys.exit(1)

        processBedFiles(sourcePath=options.sloc,destPath=options.wloc)


if __name__ == "__main__":
    main_modify_bed(sys.argv)
