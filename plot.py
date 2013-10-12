import sys
import argparse
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
from lib.plot_settings import parse_settings
from lib.hive_struct_utils import draw_hive_plot, draw_population_structure_graph
from lib.sashimi_plot_utils import draw_sashimi_plot
import os
import cPickle as pickle
import pandas

def create_data_frame_from_file(file_name):

    """ Creates a pandas.DataFrame object containing splicing expression information by reading a text file.

    file_name is the name of the file containing the splicing expression
    """

    try:
        f1 = open(file_name,'r').readlines()
        header = f1[0].strip('\n').split('\t')

        col_names = header[1:]

        construction_data_dict = {}
        
        for i in range(1,len(f1)):
            line = f1[i].strip('\n').split('\t')

            expression_values = map(float,line[2:])
            line[2:] = expression_values

            indiv_name = line[0]
            construction_data_dict[indiv_name] = pandas.Series(line[1:],col_names)

        data_frame = pandas.DataFrame(construction_data_dict).T
        return data_frame
            

    except IOError:
        print 'Invalid data file'
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw hive plots, structure plots, and sashimi plots from genomic data')

    parser.add_argument('pickle',type=str,help='Location of pickle file for plotting')
    parser.add_argument('settings',type=str,help='Location of the settings file')
    parser.add_argument('--output',type=str,default=None,help='Location of output folder. Optional argument')

    args = parser.parse_args()

    # parse the settings file
    hive_plot_settings, struct_plot_settings, sashimi_plot_settings = parse_settings(args.settings)


    try:

        # if an output directory is specified, check if it exists and create it if it doesn't
        if args.output is not None:
            try:
                os.makedirs(args.output)
            except OSError:
                if os.path.isdir(args.output):
                    pass
                else:
                    print 'Could not create directory {0}'.format(args.output)
                    raise Exception


        # extract data from pickle file
        pickle_file = open(args.pickle,'rb')
        varpos = pickle.load(pickle_file)
        junc_name = pickle.load(pickle_file)
        genotype_averages_dict = pickle.load(pickle_file)
        mRNA_info_object = pickle.load(pickle_file)
        data_frame = pickle.load(pickle_file)
        ordered_genotypes_list = pickle.load(pickle_file)

        pickle_file.close()

        plot_name_stem = args.pickle.split('/')
        plot_name_stem = plot_name_stem[len(plot_name_stem)-1]


        file_write_stem = '{0}/plots'.format(os.path.dirname(os.path.abspath(__file__)))
        if args.output is not None:
            file_write_stem = args.output.strip('/')

        if hive_plot_settings['draw_hive_plot']:
            print 'Drawing hive plot...'
            draw_hive_plot(file_name='{0}/{1}_hive.svg'.format(file_write_stem,plot_name_stem),
                    data=data_frame,
                    hive_plot_settings=hive_plot_settings,
                    genotype_ordering=ordered_genotypes_list)

        if struct_plot_settings['draw_struct_plot']:
            print 'Drawing structure plot...'
            draw_population_structure_graph(output_file_name='{0}/{1}_structure.svg'.format(file_write_stem,plot_name_stem),
                                        data=data_frame,
                                        genotypes_ordering=ordered_genotypes_list,
                                        struct_plot_settings=struct_plot_settings)

        if sashimi_plot_settings['draw_sashimi_plot']:
            print 'Drawing sashimi plot...'
            draw_sashimi_plot(output_file_path='{0}/{1}_sashimi.svg'.format(file_write_stem,plot_name_stem),
                        settings=sashimi_plot_settings,
                        var_pos=varpos,
                        average_depths_dict=genotype_averages_dict,
                        mRNAs_object=mRNA_info_object,
                        ordered_genotypes_list=ordered_genotypes_list)
        

        print 'Done!'
    except IOError:
        print 'There is no pickle file at {0}'.format(args.pickle)
        print 'Failed'
    except EOFError:
        print 'Pickle file is invalid'
        print 'Failed'
    except Exception:
        print 'Failed'
