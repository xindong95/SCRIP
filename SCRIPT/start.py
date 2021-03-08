import argparse
import os
import sys

from SCRIPT.Constants import *

def main():
    argparser = prepare_argparser()
    args = argparser.parse_args()

    subcommand  = args.subcommand

    if subcommand == "enrich":
        from SCRIPT.enrich import run
        try:
            run( args )
        except MemoryError:
            sys.exit( "MemoryError occurred.")
    elif subcommand == "setting":
        pass
    elif subcommand == "impute":
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

    add_enrich_parser(subparsers)

    return argparser

def add_enrich_parser( subparsers ):
    """Add main function 'enrich' argument parsers.
    """
    argparser_enrich = subparsers.add_parser("enrich", help="Main function.")
    
    # group for input files
    group_input = argparser_enrich.add_argument_group( "Input files arguments" )
    group_input.add_argument( "-e", "--processed_experiment", dest = "processed_experiment", type = str, required = True, 
                              help = "A processed single cell experiment anndata. Seurat or scanpy or MAESTRO. REQUIRED." )
    group_input.add_argument( "-m", "--feature_matrix", dest = "feature_matrix", type = str, required = True, 
                              help = "A cell by peak matrix. REQUIRED." )        
    # group for output files
    group_output = argparser_enrich.add_argument_group( "Output arguments" )
    group_output.add_argument( "-p", "--project", dest = "project", type = str,
                               help = "Project name, which will be used to generate output file names. DEFAULT: Random generate", 
                               default = "" )
    group_output.add_argument( "--result_store_path", dest = "result_store_path", type = str, default = 'SCRIPT_computed.h5ad',
                                help = "default SCRIPT_computed.h5ad" )
    # group for impute
    group_impute = argparser_enrich.add_argument_group( "model paramater arguments" )
    group_impute.add_argument( "--peak_method", dest = "aggregate_peak_method", type = str, default = 'group',
                                help = "group or nearest" )
    group_impute.add_argument( "--cell_number", dest = "cell_number_per_group", type = int, default = 50,
                                help = "default 50 " )
    group_impute.add_argument( "--cell_cutoff", dest = "cell_cutoff", type = int, default = 20,
                                help = "default 20 " )
    group_impute.add_argument( "--peak_confidence", dest = "cell_cutoff", type = int, default = 5,
                                help = "default 5 " )
    group_impute.add_argument( "--bg_iteration", dest = "bg_iteration", type = str, default = 'auto',
                                help = "default auto " )
    # post-processing options
    group_postprocessing = argparser_enrich.add_argument_group( "Post-processing options" )
    group_postprocessing.add_argument( "--integration", dest = "bg_iteration", action = 'store_true', default = False,
                                help = "default False" )

    group_other = argparser_enrich.add_argument_group( "Other options" )
    group_other.add_argument( "-t", '--thread', dest='n_cores', type = int, default = 16,
                              help = "default 5 " )
    group_other.add_argument( "-y", '--confirm', dest='confirm', action = 'store_false', default = True,
                              help = "default True " )
    group_other.add_argument( "--clean", dest = "clean", action = 'store_true', default = False,
                                help = "default False" )
    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted!\n")
        sys.exit(0)