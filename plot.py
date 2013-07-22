import sys
import argparse
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
from lib.plot_settings import parse_settings
from lib.hive_struct_utils import draw_hive_plot, draw_population_structure_graph
from lib.sashimi_plot_utils import draw_sashimi_plot
import os
import cPickle as pickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw hive plots, structure plots, and sashimi plots from genomic data')

    parser.add_argument('pickle',type=str,help='Location of pickle file for plotting')
    parser.add_argument('--settings',type=str,required=True,help='Location of the settings file')

    args = parser.parse_args()

    # parse the settings file
    hive_plot_settings, struct_plot_settings, sashimi_plot_settings = parse_settings(args.settings)


    try:
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

        if hive_plot_settings['draw_hive_plot']:
            print 'Drawing hive plot...'
            draw_hive_plot(file_name='{0}/plots/{1}_hive.svg'.format(os.path.dirname(os.path.abspath(__file__)),plot_name_stem),
                    data=data_frame,
                    hive_plot_settings=hive_plot_settings)

        if struct_plot_settings['draw_struct_plot']:
            print 'Drawing structure plot...'
            draw_population_structure_graph(output_file_name='{0}/plots/{1}_structure.svg'.format(os.path.dirname(os.path.abspath(__file__)),plot_name_stem),
                                        data=data_frame,
                                        genotypes_ordering=ordered_genotypes_list,
                                        struct_plot_settings=struct_plot_settings)

        if sashimi_plot_settings['draw_sashimi_plot']:
            print 'Drawing sashimi plot...'
            draw_sashimi_plot(output_file_path='{0}/plots/{1}_sashimi.svg'.format(os.path.dirname(os.path.abspath(__file__)),plot_name_stem),
                        settings=sashimi_plot_settings,
                        var_pos=varpos,
                        average_depths_dict=genotype_averages_dict,
                        mRNAs_object=mRNA_info_object,
                        ordered_genotypes_list=ordered_genotypes_list)
        

        print 'Done!'
    except EOFError:
        print 'Pickle file is invalid'
    except Exception:
        print 'Failed'
