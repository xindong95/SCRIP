#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   enrich.py
@Time    :   2021/04/16 12:34:09
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import os
import sys
import time
import shutil
import random
import scanpy as sc
import pandas as pd
from SCRIP.enrichment.bed_generation import generate_beds_by_matrix
from SCRIP.enrichment.validation import check_para
from SCRIP.enrichment.utils import EnrichRunInfo, time_estimate
from SCRIP.enrichment.calculation import score_normalization, cal_score, get_factor_source
from SCRIP.enrichment.search import search_ref_batch, read_search_result_batch
from SCRIP.utilities.utils import read_config, print_log, safe_makedirs
# from SCRIP.enhancement.enhance import determine_number_of_cells_per_group
# from SCRIP.Constants import *


def search_and_read_result(run_info, beds_path, result_path, index, n_cores):
    folder_prefix = run_info.info['project_folder']
    qpeak_length_path = os.path.join(folder_prefix, 'qpeaks_length.txt')
    qpeak_length = pd.read_csv(qpeak_length_path, sep='\t', header=None, index_col=0)/1e8
    peaks_number_path = os.path.join(index, 'peaks_number.txt')
    peaks_number = pd.read_csv(peaks_number_path, sep='\t', header=None, index_col=0)

    run_info.safe_run(search_ref_batch, [beds_path, result_path, index, n_cores], 'bed_search')
    dataset_overlap_df = run_info.safe_run_and_store(
        read_search_result_batch, [result_path, n_cores],
        os.path.join(folder_prefix, 'dataset_overlap_df.pk'),
        'dataset_overlap_df_store')
    dataset_cell_norm_df = run_info.safe_run_and_store(
        cal_score, [dataset_overlap_df, peaks_number, qpeak_length],
        os.path.join(folder_prefix, 'dataset_cell_norm_df.pk'),
        'dataset_cell_norm_df_store')
    dataset_score_source_df = run_info.safe_run_and_store(
        get_factor_source, [dataset_cell_norm_df],
        os.path.join(folder_prefix, 'dataset_score_source_df.pk'),
        'dataset_score_source_df_store')
    tf_cell_score_df = run_info.safe_run_and_store(
        score_normalization, [dataset_cell_norm_df],
        os.path.join(folder_prefix, 'tf_cell_score_df.pk'),
        'tf_cell_score_df_store')
    # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
    cell_tf_score_df = tf_cell_score_df.T
    return cell_tf_score_df, run_info, dataset_score_source_df


def enrich(cell_feature_adata, species='NA',
           project='',
           chip_index='',
           store_result_adata=True,
           n_cores=8,
           yes=False,
           clean=True):

    if project == '':
        tmp_chr_list = [chr(i) for i in range(ord("A"), ord("Z") + 1)] + [chr(i) for i in range(ord("a"), ord("z") + 1)] + [chr(i) for i in range(ord("0"), ord("9") + 1)]
        random.seed(time.time())
        tmp_prefix = str(time.time())[6:13].replace('.','') + '_' + ''.join(random.sample(tmp_chr_list, 4))
        project = 'SCRIP_' + tmp_prefix

    safe_makedirs(project)
    project_abs_path = os.path.abspath(project)
    enrich_abs_path = os.path.join(project_abs_path, 'enrichment')
    safe_makedirs(enrich_abs_path)
    ##################################
    ### check running information
    ##################################
    params = {
        'project':project_abs_path,
        'chip_index': chip_index,
        'store_result_adata': store_result_adata,
        'clean': clean,
        'species': species
    }
    run_info = EnrichRunInfo(os.path.join(enrich_abs_path, 'run_info.json'), params)
    if not run_info.check_consist(params):
        print('Not consist! Remove the folder and continue? (Y/N)')
        while True:
            confirm_info = input()
            if confirm_info == 'y' or confirm_info == 'Y':
                shutil.rmtree(enrich_abs_path)
                safe_makedirs(enrich_abs_path)
                run_info = EnrichRunInfo(os.path.join(enrich_abs_path, 'run_info.json'), params)
                break
            elif confirm_info == 'N' or confirm_info == 'n':
                print_log('Good Bye!')
                sys.exit(0)
            else:
                print('Please type Y / N.')

    beds_path = os.path.join(project, 'enrichment', 'beds')
    # total_peaks_path = os.path.join(project, 'enrichment', 'all_peaks.bed')
    qpeaks_length_path = os.path.join(project, 'enrichment', 'qpeaks_length.txt')
    chip_result_path = os.path.join(project, 'enrichment', 'ChIP_result')
    result_store_path = os.path.join(project, 'enrichment', 'SCRIP_enrichment.txt')
    safe_makedirs(beds_path)

    if species == 'hs':
        # ccre_number = 926535
        chip_factor_number = 4499
    if species == 'mm':
        # ccre_number = 339815
        chip_factor_number = 2185
    ##################################
    ### pre-check
    ##################################
    print_log('Checking parameters ...')
    check_para(cell_feature_adata, project,
               chip_index,
               beds_path, chip_result_path, result_store_path,
               yes, clean, n_cores)
    print_log('Estimating running time ...')
    elapse, future_time = time_estimate(cell_number = cell_feature_adata.shape[0],
                                        core=n_cores, chip_factor_number=chip_factor_number)
    print_log(f"It will take about {elapse} to process and finish at {future_time}.\n")

    if yes == False:
        print('Type "Y" to continue processing, "N" to abort.')
        while True:
            confirm_info = input()
            if confirm_info == 'y' or confirm_info == 'Y':
                break
            elif confirm_info == 'N' or confirm_info == 'n':
                print_log('Good Bye!')
                sys.exit(0)
            else:
                print('Please type Y / N.')
    ##################################
    ### bed generation
    ##################################
    if run_info.info['progress']['bed_generation'] == 'No':
        if os.path.exists(beds_path):
            shutil.rmtree(beds_path)
            generate_beds_by_matrix(cell_feature_adata, beds_path, qpeaks_length_path, n_cores)
        run_info.finish_stage('bed_generation')
    # affinity = get_affinity(cell_feature_adata, total_peaks_path, chip_index, ccre_number)
    ##################################
    ### Search giggle and compute enrich score
    ##################################
    cell_tf_score_df, run_info, _ = search_and_read_result(run_info, beds_path, chip_result_path, chip_index, n_cores)
    ##################################
    ### Summary results
    ##################################
    run_info.finish_stage('result_store')
    if store_result_adata == True:
        cell_tf_score_df.to_csv(result_store_path, sep='\t')
    ##################################
    ### Clean files
    ##################################
    if clean == True:
        try:
            shutil.rmtree(beds_path)
            shutil.rmtree(chip_result_path)
        except:
            pass
    return


def run_enrich(args):
    feature_matrix_path = args.feature_matrix
    species = args.species
    project = args.project
    min_cells = args.min_cells  # for removing few features cells
    min_peaks = args.min_peaks  # for removing few cells features
    max_peaks = args.max_peaks  # for removing doublet cells
    yes = args.yes
    clean = args.clean
    n_cores = args.n_cores

    CONFIG, _ = read_config()
    if species == 'hs':
        chip_index = CONFIG['index']['human_tf_index']
    elif species == 'mm':
        chip_index = CONFIG['index']['mouse_tf_index']
    else:
        pass

    print_log('Reading Files, Please wait ...')
    feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)

    sc.pp.calculate_qc_metrics(feature_matrix, percent_top=None, log1p=False, inplace=True)
    feature_mean = feature_matrix.obs.n_genes_by_counts.mean()
    feature_std = feature_matrix.obs.n_genes_by_counts.std()

    if min_cells == 'auto':
        min_cells = int(0.005 * feature_matrix.n_obs)
    else:
        min_cells = int(min_cells)
    sc.pp.filter_genes(feature_matrix, min_cells=min_cells)

    if min_peaks == 'auto':
        min_peaks = 500 if feature_mean-3*feature_std < 500 else feature_mean-3*feature_std
    else:
        min_peaks = int(min_peaks)
    sc.pp.filter_cells(feature_matrix, min_genes=min_peaks)  # filter by n_genes_by_counts

    if max_peaks == 'auto':
        max_peaks = feature_mean+3*feature_std
        sc.pp.filter_cells(feature_matrix, max_genes = max_peaks)
        # feature_matrix = feature_matrix[feature_matrix.obs.n_genes_by_counts < max_peaks, :]
    elif max_peaks == 'none':
        print_log('Skip filter by max peaks number.')
    else:
        max_peaks = int(max_peaks)
        sc.pp.filter_cells(feature_matrix, max_genes=max_peaks)
        # feature_matrix = feature_matrix[feature_matrix.obs.n_genes_by_counts < max_peaks, :]


    enrich(feature_matrix,
           species=species,
           project=project,
           chip_index=chip_index,
           store_result_adata=True,
           n_cores=n_cores,
           yes=yes,
           clean=clean
           )

    return
