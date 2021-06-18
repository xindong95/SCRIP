#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   start.py
@Time    :   2021/04/16 12:35:34
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import argparse
import os
import sys
import ruamel.yaml
from pkg_resources import resource_filename
from SCRIPT.Constants import *
from SCRIPT.utilities.utils import read_config
yaml = ruamel.yaml.YAML()

CONFIG, CONFIG_PATH = read_config()


def main():
    argparser = prepare_argparser()
    args = argparser.parse_args()

    subcommand  = args.subcommand

    if subcommand == "enrich":
        from SCRIPT.enrichment.enrich import run
        try:
            run( args )
        except MemoryError:
            sys.exit( "MemoryError occurred.")
    elif subcommand == "impute":
        pass
    elif subcommand == "target":
        pass
    elif subcommand == "setting":
        from SCRIPT.conf.settings import update_setting
        try:
            update_setting( args )
        except:
            pass
    elif subcommand == "format":
        pass

    return

def prepare_argparser():
    description = "%(prog)s"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    # top-level parser
    argparser = argparse.ArgumentParser( description = description, epilog = epilog )
    argparser.add_argument( "--version", action="version", version="%(prog)s "+SCRIPT_VERSION )
    subparsers = argparser.add_subparsers( dest = 'subcommand' )
    subparsers.required = True

    add_enrich_parser( subparsers )

    add_setting_parser( subparsers )

    return argparser

def add_setting_parser( subparsers ):
    argparser_setting = subparsers.add_parser("setting", help="Configuration.")
    argparser_setting.add_argument( "--show", dest = "show", action = 'store_true', default=False, help = "" )
    argparser_setting.add_argument( "--human_chip_index", dest = "human_chip_index", type = str, help = "" )
    argparser_setting.add_argument( "--human_motif_index", dest = "human_motif_index", type = str, help = "" )
    argparser_setting.add_argument( "--mouse_chip_index", dest = "mouse_chip_index", type = str, help = "" )
    argparser_setting.add_argument( "--mouse_motif_index", dest = "mouse_motif_index", type = str, help = "" )

    return

def add_enrich_parser( subparsers ):
    """Add main function 'enrich' argument parsers.
    """
    argparser_enrich = subparsers.add_parser("enrich", help="Main function.")
    
    # group for input files
    group_input = argparser_enrich.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-e", "--processed_experiment", dest = "processed_experiment", type = str, required = True, 
                              help = 'A processed single cell experiment anndata. Seurat or scanpy or MAESTRO. REQUIRED.' )
    group_input.add_argument( "-i", "--input_feature_matrix", dest = "feature_matrix", type = str, required = True, 
                              help = 'A cell by peak matrix . REQUIRED.' )
    group_input.add_argument( "-s", "--species", dest = "species", choices= ['hs', 'mm'], required = True, 
                              help = 'Species. "hs"(human) or "mm"(mouse). REQUIRED.' )
    # group for output files
    group_output = argparser_enrich.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str, default = "" ,
                               help = 'Project name, which will be used to generate output files folder. DEFAULT: Random generate.')
    
    # # group for preprocessing
    # group_preprocessing = argparser_enrich.add_argument_group( "Preprocessing paramater arguments" )
    # group_preprocessing.add_argument( "--cell_number", dest = "cell_number_per_group", type = str, default = 'auto',
    #                                     help='Number of cell of each group when imputation. DEFAULT: "auto".')
    # group_preprocessing.add_argument("--peak_confidence", dest="peak_confidence", type=str, default='auto',
    #                                     help='Remove peak that confidence less than this number. (Not including equal to.). Recommand 0.1*cell_number_per_group. DEFAULT: "auto".')
    
    # group for impute
    group_impute = argparser_enrich.add_argument_group( "Peak imputation paramater arguments" )
    group_impute.add_argument( "--cell_number", dest = "cell_number_per_group", type = str, default = 'auto',
                               help='Number of cell of each group when imputation. DEFAULT: "auto".')
    group_impute.add_argument("--peak_confidence", dest="peak_confidence", type=str, default='auto',
                              help='Remove peak that confidence less than this number. (Not including equal to.). Recommand 0.1*cell_number_per_group. DEFAULT: "auto".')
    # processing options
    group_processing = argparser_enrich.add_argument_group( "Processing options" )
    group_processing.add_argument( "--bg_iter", dest = "bg_iter", type = str, default = 'auto',
                                   help='Number of background iteraton. Recommand greater than 1000. DEFAULT: "auto".')
    group_processing.add_argument( "--reference", dest = "reference", choices=['integration','both','chip','motif'], default = "integration",
                                   help='Choise of integration / both / chip / motif, "integration" means integrate motif and ChIP result by standard deviation, "both" means keep the both result without merge. DEFAULT: integration. ')

    group_other = argparser_enrich.add_argument_group( "Other options" )
    group_other.add_argument( "-t", '--thread', dest='n_cores', type = int, default = 16,
                              help="Number of cores use to run SCRIPT. DEFAULT: 16.")
    group_other.add_argument( "-y", "--yes", dest='yes', action = 'store_true', default = False,
                              help="Whether ask for confirmation. DEFAULT: False.")
    group_other.add_argument( "--clean", dest = "clean", action = 'store_true', default = False,
                              help="Whether delete tmp files(including bed and search results) generated by SCRIPT. DEFAULT: False.")
    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)
