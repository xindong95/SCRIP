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
from SCRIPT.enrichment.bed_generation import generate_neighbor_bed
from SCRIPT.enrichment.validation import check_para
from SCRIPT.enrichment.utils import EnrichRunInfo, time_estimate
from SCRIPT.enrichment.post_processing import map_factor_on_ChIP, merge_giggle_adata
from SCRIPT.enrichment.calculation import zscore_normalization, cal_peak_norm, cal_score
from SCRIPT.enrichment.search_seqpare import search_seqpare_batch, read_seqpare_result_batch
from SCRIPT.utilities.utils import read_config, read_SingleCellExperiment_rds, print_log, store_to_pickle, read_pickle, safe_makedirs
from SCRIPT.imputation.impute import determine_number_of_cells_per_group
# from SCRIPT.Constants import *


def search_and_read_seqpare(run_info, tp, imputed_beds_path, result_path, peaks_number_path, index, ccre_number, n_cores):
    folder_prefix = run_info.info['project_folder']
    # tp(type) is 'ChIP-seq' or 'motif'
    if tp == 'ChIP-seq':
        run_info.safe_run(search_seqpare_batch, [imputed_beds_path, result_path, index, n_cores, tp], 'bed_ChIP_search')
        dataset_mbm_overlap_ChIP_df = run_info.safe_run_and_store(
            read_seqpare_result_batch, [result_path, n_cores, tp],
            os.path.join(folder_prefix, 'enrichment', 'dataset_mbm_overlap_ChIP_df.pk'),
            'dataset_mbm_overlap_ChIP_df_store')
        dataset_bg_peak_norm_ChIP_df = run_info.safe_run_and_store(
            cal_peak_norm, [os.path.join(index, 'peaks_number.txt'), peaks_number_path, ccre_number],
            os.path.join(folder_prefix, 'enrichment', 'dataset_bg_peak_norm_ChIP_df.pk'),
            'dataset_bg_peak_norm_ChIP_df_store')
        dataset_raw_score_ChIP_df = run_info.safe_run_and_store(
            cal_score, [dataset_mbm_overlap_ChIP_df, dataset_bg_peak_norm_ChIP_df],
            os.path.join(folder_prefix, 'enrichment', 'dataset_raw_score_ChIP_df.pk'),
            'dataset_raw_score_ChIP_df_store')
        dataset_zscore_ChIP_df = run_info.safe_run_and_store(
            zscore_normalization, [dataset_raw_score_ChIP_df, 'cell'],
            os.path.join(folder_prefix, 'enrichment', 'dataset_zscore_ChIP_df.pk'),
            'dataset_zscore_ChIP_df_store')
        # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
        fg_cell_dataset_score_df = map_factor_on_ChIP(dataset_zscore_ChIP_df).T
    else:
        run_info.safe_run(search_seqpare_batch, [imputed_beds_path, result_path, index, n_cores, tp], 'bed_motif_search')
        dataset_mbm_overlap_motif_df = run_info.safe_run_and_store(
            read_seqpare_result_batch, [result_path, n_cores, tp],
            os.path.join(folder_prefix, 'enrichment', 'dataset_mbm_overlap_motif_df.pk'),
            'dataset_mbm_overlap_motif_df_store')
        dataset_bg_peak_norm_motif_df = run_info.safe_run_and_store(
            cal_peak_norm, [os.path.join(index, 'peaks_number.txt'), peaks_number_path, ccre_number],
            os.path.join(folder_prefix, 'enrichment', 'dataset_bg_peak_norm_motif_df.pk'),
            'dataset_bg_peak_norm_motif_df_store')
        dataset_raw_score_motif_df = run_info.safe_run_and_store(
            cal_score, [dataset_mbm_overlap_motif_df, dataset_bg_peak_norm_motif_df],
            os.path.join(folder_prefix, 'enrichment', 'dataset_raw_score_motif_df.pk'),
            'dataset_raw_score_motif_df_store')
        dataset_zscore_motif_df = run_info.safe_run_and_store(
            zscore_normalization, [dataset_raw_score_motif_df, 'cell'],
            os.path.join(folder_prefix, 'enrichment', 'dataset_zscore_motif_df.pk'),
            'dataset_zscore_motif_df_store')
        # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
        fg_cell_dataset_score_df = dataset_zscore_motif_df.T
    return fg_cell_dataset_score_df, run_info


def enrich(processed_adata, cell_feature_adata, species='NA',
           project='',
           cell_number_impute=20, peak_confidence_impute=2,
           chip_index='', motif_index='',
           reference_method='integration',
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
    ##################################
    ### check running information
    ##################################
    params = {
        'project':project_abs_path,
        'cell_number_impute': cell_number_impute,
        'peak_confidence_impute': peak_confidence_impute,
        'chip_index': chip_index,
        'motif_index': motif_index,
        'reference_method': reference_method,
        'store_result_adata': store_result_adata,
        'clean': clean,
        'species': species
    }
    run_info = EnrichRunInfo(os.path.join(project,'run_info.json'), params)
    if not run_info.check_consist(params):
        print('Not consist! Remove the folder and continue? (Y/N)')
        while True:
            confirm_info = input()
            if confirm_info == 'y' or confirm_info == 'Y':
                shutil.rmtree(project)
                safe_makedirs(project)
                run_info = EnrichRunInfo(os.path.join(project,'run_info.json'), params)
                break
            elif confirm_info == 'N' or confirm_info == 'n':
                print_log('Good Bye!')
                sys.exit(0)
            else:
                print('Please type Y / N.')

    imputed_beds_path = os.path.join(project, 'enrichment', 'imputed_beds')
    imputed_beds_dict_path = os.path.join(project, 'enrichment', 'imputed_beds_barcodes.pk')
    imputed_peaks_number_path = os.path.join(project, 'enrichment', 'imputed_peaks_number.txt')
    chip_result_path = os.path.join(project, 'enrichment', 'ChIP_result')
    motif_result_path = os.path.join(project, 'enrichment', 'motif_result')
    result_store_path = os.path.join(project, 'enrichment', 'SCRIPT_enrichment.h5ad')
    safe_makedirs(imputed_beds_path)

    if species == 'hs':
        ccre_number = 926535
        chip_factor_number = 4499
        motif_factor_number = 1096
    if species == 'mm':
        ccre_number = 339815
        chip_factor_number = 2185
        motif_factor_number = 861
    ##################################
    ### pre-check
    ##################################
    print_log('Checking parameters ...')
    check_para(processed_adata, cell_feature_adata, project,
               cell_number_impute, peak_confidence_impute,
               chip_index, motif_index, reference_method, 
               imputed_beds_path, chip_result_path, motif_result_path, result_store_path,
               yes, clean, n_cores)
    print_log('Estimating running time ...')
    elapse, future_time = time_estimate(cell_number = cell_feature_adata.shape[0], reference_method=reference_method, 
                                        core=n_cores, chip_factor_number=chip_factor_number, motif_factor_number=motif_factor_number)
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
        if os.path.exists(imputed_beds_path):
            shutil.rmtree(imputed_beds_path)
        map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, 
                                        imputed_beds_path, imputed_beds_dict_path, imputed_peaks_number_path, 
                                        cell_number_impute, peak_confidence_impute, n_cores)
        run_info.finish_stage('bed_generation')
    else:
        with open(imputed_beds_dict_path, "rb") as map_dict_file:
            map_dict = pickle.load(map_dict_file)

    ##################################
    ### Search giggle and compute enrich score
    ##################################
    if reference_method in ['integration', 'both', 'chip']: # need to search chip
        search_chip = True
        fg_cell_dataset_score_df_ChIP, run_info = search_and_read_seqpare(run_info, 'ChIP-seq', imputed_beds_path, chip_result_path,
                                                                          imputed_peaks_number_path, chip_index, ccre_number, n_cores)
    if reference_method in ['integration', 'both', 'motif']: # need to search motif
        search_motif = True
        fg_cell_dataset_score_df_motif, run_info = search_and_read_seqpare(run_info, 'motif', imputed_beds_path, motif_result_path,
                                                                           imputed_peaks_number_path, motif_index, ccre_number, n_cores)
    ##################################
    ### Summary results
    ##################################
    if reference_method == 'integration':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_ChIP, 'integration', fg_cell_dataset_score_df_motif)
    elif reference_method == 'both':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_ChIP, 'ChIP-seq')
        regulation_adata = merge_giggle_adata(regulation_adata, fg_cell_dataset_score_df_motif, 'motif')
    elif reference_method == 'chip':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_ChIP, 'ChIP-seq')
    elif reference_method == 'motif':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_dataset_score_df_motif, 'motif')
    run_info.finish_stage('enrich_adata_store')
    if store_result_adata == True:
        regulation_adata.write(result_store_path)
    ##################################
    ### Clean files
    ##################################
    if clean == True:
        try:
            shutil.rmtree(imputed_beds_path)
            if search_chip == True:
                shutil.rmtree(chip_result_path)
            if search_motif == True:
                shutil.rmtree(motif_result_path)
        except:
            pass
    return regulation_adata


def run( args ):
    processed_adata_path = args.processed_experiment
    feature_matrix_path = args.feature_matrix
    species = args.species
    project = args.project
    cell_number_impute = args.cell_number_impute
    peak_confidence_impute = args.peak_confidence_impute
    reference_method = args.reference
    yes = args.yes
    clean = args.clean
    n_cores = args.n_cores

    CONFIG, _ = read_config()
    if species == 'hs':
        chip_index = CONFIG['index']['human_chip_index']
        motif_index = CONFIG['index']['human_motif_index']
    elif species == 'mm':
        chip_index = CONFIG['index']['mouse_chip_index']
        motif_index = CONFIG['index']['mouse_motif_index']
    else:
        pass
    
    print_log('Reading Files, Please wait ...')
    processed_adata = ad.read_h5ad(processed_adata_path)
    feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)

    if peak_confidence_impute == 'auto':
        if cell_number_impute == 'auto':
            print_log('Estimating the number of cells per group...')
            cell_number_impute = determine_number_of_cells_per_group(
                feature_matrix, start=1, end=70, iteration_time=30, aim_peak_number=10000, peak_confidence=None)
            print_log('Merge {cell_number} neighbor cells as one!'.format(cell_number=cell_number_impute))
        else:
            cell_number_impute = int(cell_number_impute)
        peak_confidence_impute = int(np.ceil(0.1*cell_number_impute))
    else:
        peak_confidence_impute = int(peak_confidence_impute)
        if cell_number_impute == 'auto':
            print_log('Estimating the number of cells per group...')
            cell_number_impute = determine_number_of_cells_per_group(
                feature_matrix, start=1, end=70, iteration_time=30, aim_peak_number=10000, peak_confidence=peak_confidence_impute)
            print_log('Merge {cell_number} neighbor cells as one!'.format(cell_number=cell_number_impute))
        else:
            cell_number_impute = int(cell_number_impute)
        

    enrich(processed_adata, 
           feature_matrix, 
           species=species,
           project=project,
           cell_number_impute=cell_number_impute,
           peak_confidence_impute=peak_confidence_impute,
           chip_index=chip_index, 
           motif_index=motif_index, 
           reference_method=reference_method, 
           store_result_adata=True,
           n_cores=n_cores,
           yes=yes, 
           clean=clean
           )

    return
