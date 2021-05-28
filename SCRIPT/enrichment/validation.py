#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   validation.py
@Time    :   2021/04/16 12:34:48
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

from SCRIPT.utilities.utils import print_log
import re
import os
import sys


def check_para(processed_adata, cell_feature_adata, project,
                cell_number_per_group, peak_confidence, bg_iteration, 
                chip_index, motif_index, reference_method, result_adata_path,
                yes, clean, n_cores):
    if processed_adata.shape[0] > cell_feature_adata.shape[0]:
        print_log('WARNING: There are more cells in anndata than input feature matrix.')
    if bg_iteration * cell_number_per_group < cell_feature_adata.shape[0] * 3:
        print_log('WARNING: Background iteration number less than 3 times number of cell. May be not enough to estimate background.')
    if re.match(r'chr.*_\d*_\d*', cell_feature_adata.var_names[0]) == None:
        print_log("cell_feature_adata's index/rownames should be like this: chr1_222222_333333. Each row is a feature and each column is a cell.")
        sys.exit()
    if reference_method in ['chip' , 'integration' , 'both']:
        if not os.path.exists(chip_index) or chip_index == '':
            print_log("chip_index do not exist!")
            sys.exit()
    if reference_method in ['motif', 'integration', 'both']:
        if not os.path.exists(motif_index) or motif_index == '':
            print_log("motif_index do not exist!")
            sys.exit()
    # prepare information
    information = '\n~~~~~\nWelcome to use SCRIPT. Here are your settings:\n\n'
    information += 'There are UMAP information and cluster information in your processed_adata. '
    information += 'This will borrow peaks from {cell_number_per_group} nearest cells for every cell. '.format(cell_number_per_group=cell_number_per_group)
    information += 'For each cell, peaks that confidence is less than {peak_confidence} will be deleted.\n'.format(peak_confidence=peak_confidence)
    information += 'SCRIPT randomly aggregates {cell_number_per_group} cells peaks '.format(cell_number_per_group=cell_number_per_group)
    information += 'and deletes the peaks that confidence less than {peak_confidence}. '.format(peak_confidence=peak_confidence)
    information += 'This will iterate {bg_iteration} times.\n'.format(bg_iteration=bg_iteration)
    information += 'The background peaks will be stored in {project}, '.format(project=os.path.join(project, 'enrichment', 'bg_files', 'bg_bed'))
    information += 'foreground peaks will be stored in {project}.\n'.format(project=os.path.join(project, 'enrichment', 'fg_files', 'fg_bed'))
    if reference_method == 'integration' or reference_method == 'both':
        information += 'SCRIPT will search ChIP-seq and motif index both.\n'
        if reference_method == 'integration':
            information += 'After that, SCRIPT will integrate into one result based on variability from each result.\n'
    elif reference_method == 'chip':
        information += 'SCRIPT will search ChIP-seq index.\n'
    elif reference_method == 'motif':
        information += 'SCRIPT will search motif index.\n'
    if reference_method != 'motif':
        information += 'ChIP-seq index locates at {chip_index}.\n'.format(chip_index=chip_index)
        information += 'ChIP-seq background result will be stored in {project}, '.format(project=os.path.join(project, 'enrichment', 'bg_files', 'bg_chip_result'))
        information += 'foreground result will be stored in {project}.\n'.format(project=os.path.join(project, 'enrichment', 'fg_files', 'fg_chip_result'))
    if reference_method != 'chip':
        information += 'Motif index locates at {motif_index}.\n'.format(motif_index=motif_index)
        information += 'Motif background result will be stored in {project}, '.format(project=os.path.join(project, 'enrichment', 'bg_files', 'bg_motif_result'))
        information += 'foreground result will be stored in {project}.\n'.format(project=os.path.join(project, 'enrichment', 'fg_files', 'fg_motif_result'))
    if result_adata_path != '':
        information += 'Computed result will be stored in {result_adata_path}.\n'.format(result_adata_path=result_adata_path)
    else:
        information += 'Computed result will not be stored.\n'
    information += 'All folders will{not_string} be removed after processing. '.format(not_string='' if clean == True else ' not')
    information += 'All processes will use {n_cores} cores.\n~~~~~\n'.format(n_cores=n_cores)
    print(information)
