import sys
import argparse
from lib.ReadDepth import ReadDepth
from lib.mRNAsObject import mRNAsObject
import anydbm as dbm
from lib.parse_settings import parse_settings
from lib.hive_struct_utils import draw_hive_plot, draw_population_structure_graph
from lib.sashimi_plot_utils import 
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw hive plots, structure plots, and sashimi plots from genomic data')

    parser.add_argument('db',type=str,help='Location of dbm file for plotting')
    parser.add_argument('--settings',type=str,required=True,help='Location of the settings file')

    args = parser.parse_args()

    # parse the settings file
    hive_plot_settings, struct_plot_settings, sashimi_plot_settings = parse_settings(args.settings)

    # extract data from database file
    db = dbm.open(args.db,'r')
    varpos = db['varpos']
    junc_name = db['junc_name']
    genotype_averages_dict = db['averages']
    data_frame = db['df']
    mRNA_info_object = db['mRNA']
    db.close()

    if hive_plot_settings['draw_hive_plot']:
        draw_hive_plot(file_name='{0}/{1}_hive.svg'.format(os.path.dirname(__file__),db),
                data=data_frame,
                hive_plot_settings=hive_plot_settings)

    if struct_plot_settings['draw_struct_plot']:
        draw_population_structure_graph(output_file_name='{0}/{1}_structure.svg'.format(os.path.dirname(__file__),db),
                                    data=data_frame,
                                    struct_plot_settings=struct_plot_settings)

    if sashimi_plot_settings['draw_sashimi_plot']:
        draw_sashimi_plot(output_file_path='{0}/{1}_structure.svg'.format(os.path.dirname(__file__),db),
                    settings=sashimi_plot_settings,
                    var_pos=varpos,
                    average_depths_dict=genotype_averages_dict,
                    mRNAs_object=mRNA_info_object)
    


