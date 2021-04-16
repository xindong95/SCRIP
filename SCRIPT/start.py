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
                              help = "A processed single cell experiment anndata. Seurat or scanpy or MAESTRO. REQUIRED." )
    group_input.add_argument( "-i", "--input_feature_matrix", dest = "feature_matrix", type = str, required = True, 
                              help = "A cell by peak matrix . REQUIRED." )
    group_input.add_argument( "-s", "--species", dest = "species", choices= ['hs', 'mm'], required = True, 
                              help = "species. REQUIRED." )
    # group for output files
    group_output = argparser_enrich.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str,
                               help = "Project name, which will be used to generate output files folder. DEFAULT: Random generate", 
                               default = "" )
    # group for impute
    group_impute = argparser_enrich.add_argument_group( "Peak aggregation paramater arguments" )
    group_impute.add_argument( "--peak_method", dest = "aggregate_peak_method", type = str, default = 'group',
                                help = "group or nearest" )
    group_impute.add_argument( "--cell_number", dest = "cell_number_per_group", type = int, default = 50,
                                help = "default 50 " )
    group_impute.add_argument( "--cell_cutoff", dest = "cell_cutoff", type = int, default = 20,
                                help = "default 20 " )
    group_impute.add_argument( "--peak_confidence", dest = "peak_confidence", type = int, default = 5,
                                help = "default 5 " )
    # processing options
    group_processing = argparser_enrich.add_argument_group( "Processing options" )
    group_processing.add_argument( "--bg_iter", dest = "bg_iter", type = str, default = 'auto',
                                help = "default auto" )
    group_processing.add_argument( "--reference", dest = "reference", choices=['integration','both','chip','motif'], default = "integration",
                                help = "integration / both / chip / motif" )

    group_other = argparser_enrich.add_argument_group( "Other options" )
    group_other.add_argument( "-t", '--thread', dest='n_cores', type = int, default = 16,
                              help = "default 16 " )
    group_other.add_argument( "-y", "--yes", dest='yes', action = 'store_true', default = False,
                              help = "default False" )
    group_other.add_argument( "--clean", dest = "clean", action = 'store_true', default = False,
                                help = "default False" )
    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)