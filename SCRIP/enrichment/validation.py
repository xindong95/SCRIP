#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   validation.py
@Time    :   2021/04/16 12:34:48
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

from SCRIP.utilities.utils import print_log
import re
import os
import sys


def check_para(cell_feature_adata, project,
               chip_index,  
               beds_path, chip_result_path, result_store_path,
               yes, clean, n_cores):
    # if processed_adata.shape[0] > cell_feature_adata.shape[0]:
    #     print_log('WARNING: There are more cells in anndata than input feature matrix.')
    if re.match(r'chr.*_\d*_\d*', cell_feature_adata.var_names[0]) == None:
        print_log("cell_feature_adata's index/rownames should be like this: chr1_222222_333333. Each row is a feature and each column is a cell.")
        sys.exit()
    # if reference_method in ['chip' , 'integration' , 'both']:
    #     if not os.path.exists(chip_index) or chip_index == '':
    #         print_log("chip_index do not exist!")
    #         sys.exit()
    # if reference_method in ['motif', 'integration', 'both']:
    #     if not os.path.exists(motif_index) or motif_index == '':
    #         print_log("motif_index do not exist!")
    #         sys.exit()
    # prepare information
    information = '\n~~~~~\nWelcome to use SCRIP. Here are your settings:\n\n'
    information += 'There are UMAP information and cluster information in your processed_adata. '
    information += 'The whole output file will be stored in {project}. '.format(project=project)
    # information += 'This will borrow peaks from {cell_number_impute} nearest cells for every cell. '.format(cell_number_impute=cell_number_impute)
    # information += 'For each cell, peaks that confidence is less than {peak_confidence_impute} will be deleted.\n'.format(peak_confidence_impute=peak_confidence_impute)
    information += 'The peaks will be stored in {beds_path}.\n'.format(beds_path=beds_path)
    # if reference_method == 'integration' or reference_method == 'both':
    #     information += 'SCRIP will search ChIP-seq and motif index both.\n'
    #     if reference_method == 'integration':
    #         information += 'After that, SCRIP will integrate into one result based on variability from each result.\n'
    # elif reference_method == 'chip':
    #     information += 'SCRIP will search ChIP-seq index.\n'
    # elif reference_method == 'motif':
    #     information += 'SCRIP will search motif index.\n'
    # if reference_method != 'motif':
    #     information += 'ChIP-seq index locates at {chip_index}.\n'.format(chip_index=chip_index)
    #     information += 'ChIP-seq intervals intersection result will be stored in {chip_result_path}'.format(chip_result_path=chip_result_path)
    # if reference_method != 'chip':
    #     information += 'Motif index locates at {motif_index}.\n'.format(motif_index=motif_index)
    #     information += 'Motif intervals intersection result will be stored in {motif_result_path}, '.format(motif_result_path=motif_result_path)
    if result_store_path != '':
        information += 'Computed result will be stored in {result_store_path}.\n'.format(result_store_path=result_store_path)
    else:
        information += 'Computed result will not be stored.\n'
    information += 'All folders will{not_string} be removed after processing. '.format(not_string='' if clean == True else ' not')
    information += 'All processes will use {n_cores} cores.\n~~~~~\n'.format(n_cores=n_cores)
    print(information)
