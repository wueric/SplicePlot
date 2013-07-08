import sys
import argparse
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
import cPickle as pickle
import initialize_data
from lib.parse_settings import parse_settings

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw hive plots, structure plots, and sashimi plots from genomic data')
    parser.add_argument('-p','--pickle',type=str,help='Location of pickle file for plotting')
    parser.add_argument('--gtf',type=str,help='Location of the gtf annotation file')
    parser.add_argument('--vcf',type=str,help='Location of the vcf file')
    parser.add_argument('--mf',type=str,help='Location of the map file')
    parser.add_argument('--varpos',type=str,help='Chromosome and base number of SNP')
    parser.add_argument('--junc',type=str,help='Name of the junction')
    parser.add_argument('--settings',type=str,required=True,help='Location of the settings file')

    args = parser.parse_args()

    if not args.pickle and not args.gtf and not args.vcf \
        and not args.mf and not args.varpos and not args.junc:
        print "Inadequate number of inputs. Must include either a pickle file location or gtf, vcf, mapping file, variant name, and junction name"
        sys.exit(1)

    # parse the settings file
    hive_plot_settings, struct_plot_settings, sashimi_plot_settings = parse_settings(args.settings)

    if hive_plot_settings['draw_hive_plot']:
        

    if struct_plot_settings['draw_struct_plot']:

    if sashimi_plot_settings['draw_sashimi_plot']:
