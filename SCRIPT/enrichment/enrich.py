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
import pickle
import shutil
import random
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from SCRIPT.enrichment.bed_generation import generate_beds_by_matrix
from SCRIPT.enrichment.validation import check_para
from SCRIPT.enrichment.utils import EnrichRunInfo, time_estimate
from SCRIPT.enrichment.post_processing import merge_giggle_adata
from SCRIPT.enrichment.calculation import score_normalization, cal_peak_norm, cal_score, get_factor_source
from SCRIPT.enrichment.search_seqpare import search_seqpare_batch, read_seqpare_result_batch
from SCRIPT.utilities.utils import read_config, read_SingleCellExperiment_rds, print_log, store_to_pickle, read_pickle, safe_makedirs
from SCRIPT.enhancement.enhance import determine_number_of_cells_per_group
# from SCRIPT.Constants import *


def search_and_read_seqpare(run_info, beds_path, result_path, peaks_number_path, index, ccre_number, n_cores):
    folder_prefix = run_info.info['project_folder']
    # tp(type) is 'ChIP-seq' or 'motif'

    run_info.safe_run(search_seqpare_batch, [beds_path, result_path, index, n_cores], 'bed_search')
    dataset_mbm_overlap_df = run_info.safe_run_and_store(
        read_seqpare_result_batch, [result_path, n_cores],
        os.path.join(folder_prefix, 'dataset_mbm_overlap_df.pk'),
        'dataset_mbm_overlap_df_store')
    dataset_bg_peak_norm_df = run_info.safe_run_and_store(
        cal_peak_norm, [os.path.join(index, 'peaks_number.txt'), peaks_number_path, ccre_number],
        os.path.join(folder_prefix, 'dataset_bg_peak_norm_df.pk'),
        'dataset_bg_peak_norm_df_store')
    dataset_raw_score_df = run_info.safe_run_and_store(
        cal_score, [dataset_mbm_overlap_df, dataset_bg_peak_norm_df],
        os.path.join(folder_prefix, 'dataset_raw_score_df.pk'),
        'dataset_raw_score_df_store')
    dataset_score_resource_df = run_info.safe_run_and_store(
        get_factor_source, [dataset_raw_score_df],
        os.path.join(folder_prefix, 'dataset_score_resource_df.pk'),
        'dataset_score_resource_df_store')
    dataset_score_df = run_info.safe_run_and_store(
        score_normalization, [dataset_raw_score_df],
        os.path.join(folder_prefix, 'dataset_score_df.pk'),
        'dataset_score_df_store')
    # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
    fg_cell_dataset_score_df = dataset_score_df.T
    # else:
    #     run_info.safe_run(search_seqpare_batch, [beds_path, result_path, index, n_cores, tp], 'bed_motif_search')
    #     dataset_mbm_overlap_motif_df = run_info.safe_run_and_store(
    #         read_seqpare_result_batch, [result_path, n_cores, tp],
    #         os.path.join(folder_prefix, 'enrichment', 'dataset_mbm_overlap_motif_df.pk'),
    #         'dataset_mbm_overlap_motif_df_store')
    #     dataset_bg_peak_norm_motif_df = run_info.safe_run_and_store(
    #         cal_peak_norm, [os.path.join(index, 'peaks_number.txt'), peaks_number_path, ccre_number],
    #         os.path.join(folder_prefix, 'enrichment', 'dataset_bg_peak_norm_motif_df.pk'),
    #         'dataset_bg_peak_norm_motif_df_store')
    #     dataset_raw_score_motif_df = run_info.safe_run_and_store(
    #         cal_score, [dataset_mbm_overlap_motif_df, dataset_bg_peak_norm_motif_df],
    #         os.path.join(folder_prefix, 'enrichment', 'dataset_raw_score_motif_df.pk'),
    #         'dataset_raw_score_motif_df_store')
    #     dataset_zscore_motif_df = run_info.safe_run_and_store(
    #         score_normalization, [dataset_raw_score_motif_df, 'motif'],
    #         os.path.join(folder_prefix, 'enrichment', 'dataset_score_motif_df.pk'),
    #         'dataset_score_motif_df_store')
    #     # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
    #     fg_cell_dataset_score_df = dataset_zscore_motif_df.T
    return fg_cell_dataset_score_df, run_info


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
        project = 'SCRIPT_' + tmp_prefix

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
    peaks_number_path = os.path.join(project, 'enrichment', 'peaks_number.txt')
    chip_result_path = os.path.join(project, 'enrichment', 'ChIP_result')
    result_store_path = os.path.join(project, 'enrichment', 'SCRIPT_enrichment.txt')
    safe_makedirs(beds_path)

    if species == 'hs':
        ccre_number = 926535
        chip_factor_number = 4499
        # motif_factor_number = 1096
    if species == 'mm':
        ccre_number = 339815
        chip_factor_number = 2185
        # motif_factor_number = 861
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
    print_log("It will take about {elapse} to process and finish at {future_time}.\n".format(elapse = elapse, future_time = future_time))

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
            generate_beds_by_matrix(cell_feature_adata,beds_path, peaks_number_path, n_cores)
        # map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, 
        #                                 beds_path, imputed_beds_dict_path, imputed_peaks_number_path, 
        #                                 cell_number_impute, peak_confidence_impute, n_cores)
        run_info.finish_stage('bed_generation')
    # else:
    #     with open(imputed_beds_dict_path, "rb") as map_dict_file:
    #         map_dict = pickle.load(map_dict_file)

    ##################################
    ### Search giggle and compute enrich score
    ##################################
    # if reference_method in ['integration', 'both', 'chip']: # need to search chip
    #     search_chip = True
    fg_cell_dataset_score_df, run_info = search_and_read_seqpare(run_info, beds_path, chip_result_path,
                                                                      peaks_number_path, chip_index, ccre_number, n_cores)
    # if reference_method in ['integration', 'both', 'motif']: # need to search motif
    #     search_motif = True
    #     fg_cell_dataset_score_df_motif, run_info = search_and_read_seqpare(run_info, 'motif', beds_path, motif_result_path,
    #                                                                        imputed_peaks_number_path, motif_index, ccre_number, n_cores)
    ##################################
    ### Summary results
    ##################################
    # if reference_method == 'integration':
    #     regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_ChIP, 'integration', fg_cell_dataset_score_df_motif)
    # elif reference_method == 'both':
    #     regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_ChIP, 'ChIP-seq')
    #     regulation_adata = merge_giggle_adata(regulation_adata, fg_cell_dataset_score_df_motif, 'motif')
    # elif reference_method == 'chip':
    #     regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_ChIP, 'ChIP-seq')
    # elif reference_method == 'motif':
    #     regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_motif, 'motif')
    run_info.finish_stage('result_store')
    if store_result_adata == True:
        fg_cell_dataset_score_df.to_csv(result_store_path, sep='\t')
    ##################################
    ### Clean files
    ##################################
    if clean == True:
        try:
            shutil.rmtree(beds_path)
            # if search_chip == True:
            shutil.rmtree(chip_result_path)
            # if search_motif == True:
            #     shutil.rmtree(motif_result_path)
        except:
            pass
    return 


def run( args ):
    # processed_adata_path = args.processed_experiment
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
    else:
        max_peaks = int(max_peaks)
    feature_matrix = feature_matrix[feature_matrix.obs.n_genes_by_counts < max_peaks, :]

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
