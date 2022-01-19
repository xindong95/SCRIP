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
import sys
import ruamel.yaml
from SCRIP.Constants import SCRIP_VERSION
from SCRIP.enrichment.enrich import run_enrich
from SCRIP.imputation.impute import run_impute
from SCRIP.targets.target import run_target
from SCRIP.conf.config import update_setting
from SCRIP.utilities.utils import read_config
yaml = ruamel.yaml.YAML()

CONFIG, CONFIG_PATH = read_config()

def main():
    argparser = prepare_argparser()
    args = argparser.parse_args()

    subcommand  = args.subcommand

    if subcommand == "enrich":
        try:
            run_enrich(args)
        except MemoryError:
            sys.exit( "MemoryError occurred.")
    elif subcommand == "impute":
        try:
            run_impute(args)
        except MemoryError:
            sys.exit("MemoryError occurred.")
    elif subcommand == "target":
        try:
            run_target(args)
        except MemoryError:
            sys.exit("MemoryError occurred.")
    elif subcommand == "config":
        try:
            update_setting( args )
        except:
            sys.exit("Setting set wrong.")
    elif subcommand == "format":
        pass

    return

def prepare_argparser():
    description = "%(prog)s"
    epilog = "For command line options of each command, type: %(prog)s COMMAND -h"

    # top-level parser
    argparser = argparse.ArgumentParser( description = description, epilog = epilog )
    argparser.add_argument( "--version", action="version", version="%(prog)s "+SCRIP_VERSION )
    subparsers = argparser.add_subparsers( dest = 'subcommand' )
    subparsers.required = True

    add_enrich_parser(subparsers)
    add_impute_parser(subparsers)
    add_target_parser(subparsers)
    add_config_parser(subparsers)

    return argparser


def add_config_parser(subparsers):
    argparser_setting = subparsers.add_parser("config", help="Configuration.")
    argparser_setting.add_argument( "--show", dest = "show", action = 'store_true', default=False, help = "" )
    argparser_setting.add_argument("--human_tf_index", dest="human_tf_index", type=str, help="")
    argparser_setting.add_argument("--human_hm_index", dest="human_hm_index", type=str, help="")
    argparser_setting.add_argument("--mouse_tf_index", dest="mouse_tf_index", type=str, help="")
    argparser_setting.add_argument("--mouse_hm_index", dest="mouse_hm_index", type=str, help="")


def add_enrich_parser( subparsers ):
    """Add main function 'enrich' argument parsers.
    """
    argparser_enrich = subparsers.add_parser("enrich", help="Main function.")

    # group for input files
    group_input = argparser_enrich.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-i", "--input_feature_matrix", dest = "feature_matrix", type = str, required = True,
                              help = 'A cell by peak matrix . REQUIRED.' )
    group_input.add_argument( "-s", "--species", dest = "species", choices= ['hs', 'mm'], required = True,
                              help = 'Species. "hs"(human) or "mm"(mouse). REQUIRED.' )
    # group for output files
    group_output = argparser_enrich.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str, default = "" ,
                               help = 'Project name, which will be used to generate output files folder. DEFAULT: Random generate.')

    # group for preprocessing
    group_preprocessing = argparser_enrich.add_argument_group( "Preprocessing paramater arguments" )
    group_preprocessing.add_argument( "--min_cells", dest = "min_cells", type = str, default = 'auto',
                                      help='Minimal cell cutoff for features. Auto will take 0.05%% of total cell number.DEFAULT: "auto".')
    group_preprocessing.add_argument("--min_peaks", dest="min_peaks", type=str, default='auto',
                                     help='Minimal peak cutoff for cells. Auto will take the mean-3*std of all feature number (if less than 500 is 500). DEFAULT: "auto".')
    group_preprocessing.add_argument("--max_peaks", dest="max_peaks", type=str, default='auto',
                                     help='Max peak cutoff for cells. This will help you to remove the doublet cells. Auto will take the mean+5*std of all feature number. DEFAULT: "auto".')

    group_other = argparser_enrich.add_argument_group( "Other options" )
    group_other.add_argument( "-t", '--thread', dest='n_cores', type = int, default = 16,
                              help="Number of cores use to run SCRIP. DEFAULT: 16.")
    group_other.add_argument( "-y", "--yes", dest='yes', action = 'store_true', default = False,
                              help="Whether ask for confirmation. DEFAULT: False.")
    group_other.add_argument( "--clean", dest = "clean", action = 'store_true', default = False,
                              help="Whether delete tmp files(including bed and search results) generated by SCRIP. DEFAULT: False.")


def add_impute_parser(subparsers):
    """Add main function 'impute' argument parsers.
    """
    argparser_impute = subparsers.add_parser("impute", help="Imputation Factor function.")

    # group for input files
    group_input = argparser_impute.add_argument_group("Input files arguments")
    group_input.add_argument("-i", "--input_feature_matrix", dest="feature_matrix", type=str, required=True,
                             help='A cell by peak matrix. h5 or h5ad supported. REQUIRED.')
    group_input.add_argument("-s", "--species", dest="species", choices=['hs', 'mm'], required=True,
                             help='Species. "hs"(human) or "mm"(mouse). REQUIRED.')

    # group for output files
    group_output = argparser_impute.add_argument_group("Output arguments")
    group_output.add_argument("-p", "--project", dest="project", type=str, default="",
                              help='Project name, which will be used to generate output files folder. DEFAULT: Random generate.')
    group_output.add_argument("-f", "--format", dest="format", choices=['h5ad', 'mtx'], default="h5ad",
                              help='Format generate for output RP count. DEFAULT: h5ad.')

    # group for impute
    group_impute = argparser_impute.add_argument_group("Peak imputation paramater arguments")
    group_impute.add_argument("--factor", dest="factor", type=str, default='', required=True,
                              help='The factor you want to impute. REQUIRED.')
    group_impute.add_argument("--ref_baseline", dest="ref_baseline", type=int, default=500,
                              help='Remove dataset which peaks number less than this value. DEFAULT: 500.')
    group_impute.add_argument("--remove_others", dest="remove_others", action='store_true',  default=False,
                              help='Remove signal not from best match. DEFAULT: False.')

    group_other = argparser_impute.add_argument_group("Other options")
    group_other.add_argument("--min_cells", dest="min_cells", type=str, default='auto',
                                     help='Minimal cell cutoff for features. Auto will take 0.05%% of total cell number.DEFAULT: "auto".')
    group_other.add_argument("--min_peaks", dest="min_peaks", type=str, default='auto',
                                     help='Minimal peak cutoff for cells. Auto will take the mean-3*std of all feature number (if less than 500 is 500). DEFAULT: "auto".')
    group_other.add_argument("--max_peaks", dest="max_peaks", type=str, default='auto',
                                     help='Max peak cutoff for cells. This will help you to remove the doublet cells. Auto will take the mean+5*std of all feature number. DEFAULT: "auto".')
    group_other.add_argument("-t", '--thread', dest='n_cores', type=int, default=16,
                             help="Number of cores use to run SCRIP. DEFAULT: 16.")


def add_target_parser(subparsers):
    """Add main function 'target' argument parsers.
    """
    argparser_impute = subparsers.add_parser("target", help="Calculate targets based on factor peak count.")

    # group for input files
    group_input = argparser_impute.add_argument_group("Input files arguments")
    group_input.add_argument("-i", "--input_feature_matrix", dest="feature_matrix", type=str, required=True,
                             help='A cell by peak matrix. h5 or h5ad supported. REQUIRED.')
    group_input.add_argument("-s", "--species", dest="species", choices=['hs', 'mm'], required=True,
                             help='Species. "hs"(human) or "mm"(mouse). REQUIRED.')

    # group for output files
    group_output = argparser_impute.add_argument_group("Output arguments")
    group_output.add_argument("-o", "--output", dest="output", type=str, default="RP.h5ad",
                              help='output h5ad file. DEFAULT: RP.h5ad')

    group_other = argparser_impute.add_argument_group("Other options")
    group_other.add_argument("-d", '--decay', dest='decay', type=int, default=10000,
                             help="Range to the effect of peaks. DEFAULT: 10000.")

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)
