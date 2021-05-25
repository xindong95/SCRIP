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
from SCRIPT.enrichment.bed_generation import generate_background_bed, generate_neighbor_bed
from SCRIPT.enrichment.validation import check_para
from SCRIPT.enrichment.utils import EnrichRunInfo, time_estimate
from SCRIPT.enrichment.post_processing import map_factor_on_ChIP, merge_giggle_adata
from SCRIPT.enrichment.calculation import cal_deviation_table_batch, score_normalization
from SCRIPT.enrichment.search_giggle import search_giggle_batch, read_giggle_result_odds_batch, read_giggle_result_fisher_batch
from SCRIPT.utilities.utils import read_config, read_SingleCellExperiment_rds, print_log, store_to_pickle, read_pickle, safe_makedirs
from SCRIPT.imputation.impute import determine_number_of_cells_per_group
from SCRIPT.Constants import *

def search_and_read_giggle(run_info, tp, bg_bed_path, bg_result_path, fg_bed_path, fg_result_path, fg_peaks_number_path, index, genome_length, n_cores, fg_map_dict):
    # tp(type) is 'ChIP-seq' or 'motif'
    if tp == 'ChIP-seq':
        if run_info.info['progress']['bg_bed_chip_search'] == 'No':
            search_giggle_batch(bg_bed_path, bg_result_path, index, genome_length, n_cores, tp)
            run_info.finish_stage('bg_bed_chip_search')
        if run_info.info['progress']['fg_bed_chip_search'] == 'No':
            search_giggle_batch(fg_bed_path, fg_result_path, index, genome_length, n_cores, tp)
            run_info.finish_stage('fg_bed_chip_search')
    else:
        if run_info.info['progress']['bg_bed_motif_search'] == 'No':
            search_giggle_batch(bg_bed_path, bg_result_path, index, genome_length, n_cores, tp)
            run_info.finish_stage('bg_bed_motif_search')
        if run_info.info['progress']['fg_bed_motif_search'] == 'No':
            search_giggle_batch(fg_bed_path, fg_result_path, index, genome_length, n_cores, tp)
            run_info.finish_stage('fg_bed_motif_search')

    if tp == 'ChIP-seq':
        bg_dataset_odds_ratio_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'bg_files', 'bg_dataset_odds_ratio_df_ChIP.pk')
        if run_info.info['progress']['bg_dataset_raw_odds_ratio_chip_df_store'] == 'No':
            bg_dataset_odds_ratio_df = read_giggle_result_odds_batch(bg_result_path, n_cores, 'background {tp}'.format(tp=tp))
            store_to_pickle(bg_dataset_odds_ratio_df, bg_dataset_odds_ratio_df_path)
            run_info.finish_stage('bg_dataset_raw_odds_ratio_chip_df_store')
        else:
            bg_dataset_odds_ratio_df = read_pickle(bg_dataset_odds_ratio_df_path)

        fg_dataset_odds_ratio_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_odds_ratio_df_ChIP.pk')
        if run_info.info['progress']['fg_dataset_raw_odds_ratio_chip_df_store'] == 'No':
            fg_dataset_odds_ratio_df = read_giggle_result_odds_batch(fg_result_path, n_cores, 'foreground {tp}'.format(tp=tp))
            store_to_pickle(fg_dataset_odds_ratio_df, fg_dataset_odds_ratio_df_path)
            run_info.finish_stage('fg_dataset_raw_odds_ratio_chip_df_store')
        else:
            fg_dataset_odds_ratio_df = read_pickle(fg_dataset_odds_ratio_df_path)

        fg_dataset_fisher_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_fisher_df_ChIP.pk')
        if run_info.info['progress']['fg_dataset_fisher_chip_df_store'] == 'No':
            fg_dataset_fisher_df = read_giggle_result_fisher_batch(fg_result_path, n_cores, 'foreground {tp}'.format(tp=tp))
            store_to_pickle(fg_dataset_fisher_df, fg_dataset_fisher_df_path)
            run_info.finish_stage('fg_dataset_fisher_chip_df_store')
        else:
            fg_dataset_fisher_df = read_pickle(fg_dataset_fisher_df_path)

        fg_dataset_deviation_score_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_fisher_df_ChIP.pk')
        if run_info.info['progress']['fg_dataset_deviation_score_chip_df_store'] == 'No':
            fg_dataset_deviation_score_df = cal_deviation_table_batch(fg_dataset_odds_ratio_df, bg_dataset_odds_ratio_df, n_cores)
            store_to_pickle(fg_dataset_deviation_score_df, fg_dataset_deviation_score_df_path)
            run_info.finish_stage('fg_dataset_deviation_score_chip_df_store')
        else:
            fg_dataset_deviation_score_df = read_pickle(fg_dataset_deviation_score_df_path)

        fg_dataset_cell_score_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_score_df_ChIP.pk')
        if run_info.info['progress']['fg_dataset_score_chip_df_store'] == 'No':
            fg_dataset_cell_score_df = score_normalization(fg_dataset_deviation_score_df, fg_dataset_fisher_df, os.path.join(index, 'peaks_number.txt'), fg_peaks_number_path)
            store_to_pickle(fg_dataset_cell_score_df, fg_dataset_cell_score_df_path)
            run_info.finish_stage('fg_dataset_score_chip_df_store')
        else:
            fg_dataset_cell_score_df = read_pickle(fg_dataset_cell_score_df_path)
    else:
        bg_dataset_odds_ratio_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'bg_files', 'bg_dataset_odds_ratio_df_motif.pk')
        if run_info.info['progress']['bg_dataset_odds_ratio_motif_df_store'] == 'No':
            bg_dataset_odds_ratio_df = read_giggle_result_odds_batch(bg_result_path, n_cores, 'background {tp}'.format(tp=tp))
            store_to_pickle(bg_dataset_odds_ratio_df, bg_dataset_odds_ratio_df_path)
            run_info.finish_stage('bg_dataset_odds_ratio_motif_df_store')
        else:
            bg_dataset_odds_ratio_df = read_pickle(bg_dataset_odds_ratio_df_path)

        fg_dataset_odds_ratio_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_odds_ratio_df_motif.pk')
        if run_info.info['progress']['fg_dataset_odds_ratio_motif_df_store'] == 'No':
            fg_dataset_odds_ratio_df = read_giggle_result_odds_batch(fg_result_path, n_cores, 'foreground {tp}'.format(tp=tp))
            store_to_pickle(fg_dataset_odds_ratio_df, fg_dataset_odds_ratio_df_path)
            run_info.finish_stage('fg_dataset_odds_ratio_motif_df_store')
        else:
            fg_dataset_odds_ratio_df = read_pickle(fg_dataset_odds_ratio_df_path)

        fg_dataset_fisher_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_fisher_df_motif.pk')
        if run_info.info['progress']['fg_dataset_fisher_motif_df_store'] == 'No':
            fg_dataset_fisher_df = read_giggle_result_fisher_batch(fg_result_path, n_cores, 'foreground {tp}'.format(tp=tp))
            store_to_pickle(fg_dataset_fisher_df, fg_dataset_fisher_df_path)
            run_info.finish_stage('fg_dataset_fisher_motif_df_store')
        else:
            fg_dataset_fisher_df = read_pickle(fg_dataset_fisher_df_path)

        fg_dataset_deviation_score_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_fisher_df_motif.pk')
        if run_info.info['progress']['fg_dataset_deviation_score_motif_df_store'] == 'No':
            fg_dataset_deviation_score_df = cal_deviation_table_batch(fg_dataset_odds_ratio_df, bg_dataset_odds_ratio_df, n_cores)
            store_to_pickle(fg_dataset_deviation_score_df, fg_dataset_deviation_score_df_path)
            run_info.finish_stage('fg_dataset_deviation_score_motif_df_store')
        else:
            fg_dataset_deviation_score_df = read_pickle(fg_dataset_deviation_score_df_path)

        fg_dataset_cell_score_df_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_score_df_motif.pk')
        if run_info.info['progress']['fg_dataset_score_motif_df_store'] == 'No':
            fg_dataset_cell_score_df = score_normalization(fg_dataset_deviation_score_df, fg_dataset_fisher_df, os.path.join(index, 'peaks_number.txt'), fg_peaks_number_path)
            store_to_pickle(fg_dataset_cell_score_df, fg_dataset_cell_score_df_path)
            run_info.finish_stage('fg_dataset_score_motif_df_store')
        else:
            fg_dataset_cell_score_df = read_pickle(fg_dataset_cell_score_df_path)

    # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
    if tp == 'ChIP-seq':
        fg_cell_dataset_score_df = map_factor_on_ChIP(fg_dataset_cell_score_df).T
    else:
        # motif's dataset is same to factor
        fg_cell_dataset_score_df = fg_dataset_cell_score_df.T
    return fg_cell_dataset_score_df, run_info


def enrich(processed_adata, cell_feature_adata, project='',
           cell_number_per_group=20, peak_confidence=2, bg_iteration='auto', 
           chip_index='', motif_index='', reference_method='integration', store_result_adata=True,
           yes=False, clean=True, n_cores=8, 
           processed_adata_path='NA', cell_feature_adata_path='NA', species='NA'):

    if project == '':
        tmp_chr_list = [chr(i) for i in range(ord("A"), ord("Z") + 1)] + [chr(i) for i in range(ord("a"), ord("z") + 1)] + [chr(i) for i in range(ord("0"), ord("9") + 1)]
        random.seed(time.time())
        tmp_prefix = str(time.time())[6:13].replace('.','') + '_' + ''.join(random.sample(tmp_chr_list, 4))
        project = 'SCRIPT_' + tmp_prefix

    safe_makedirs(project)
    project_abs_path = os.path.abspath(project)
    if bg_iteration != 'auto':
        bg_iteration = int(bg_iteration)
    ##################################
    ### check running information
    ##################################
    params = {
        'project':project_abs_path,
        'cell_number_per_group':cell_number_per_group,
        'peak_confidence':peak_confidence,
        'bg_iteration':bg_iteration,
        'chip_index':chip_index,
        'motif_index':motif_index,
        'reference_method':reference_method,
        'store_result_adata':store_result_adata,
        'clean':clean,
        'processed_adata_path':os.path.abspath(processed_adata_path),
        'cell_feature_adata_path':os.path.abspath(cell_feature_adata_path),
        'species':species
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

    fg_bed_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_bed')
    fg_map_dict_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_bed.pk')
    fg_peaks_number_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_peaks_number.txt')
    fg_motif_result_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_motif_result')
    fg_chip_result_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_chip_result')
    bg_bed_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_bed')
    bg_map_dict_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_bed.pk')
    bg_chip_result_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_chip_result') 
    bg_motif_result_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_motif_result')
    result_store_path = os.path.join(project, 'enrichment', 'SCRIPT_enrichment.h5ad')
    safe_makedirs(fg_bed_path)
    safe_makedirs(bg_bed_path)

    if species == 'hs':
        genome_length = '3088269832'
    if species == 'mm':
        genome_length = '2730855475'
    ##################################
    ### pre-check
    ##################################
    print_log('Checking parameters ...')
    if bg_iteration == "auto":
        bg_iteration = int(processed_adata.shape[0] * 3 / cell_number_per_group) + 1
    check_para(processed_adata, cell_feature_adata, project,
                cell_number_per_group, peak_confidence, bg_iteration, 
                chip_index, motif_index, reference_method, result_store_path,
                yes, clean, n_cores)
    print_log('Estimating running time ...')
    elapse, future_time = time_estimate(cell_number = cell_feature_adata.shape[0], bg_iteration_number=bg_iteration,                                         
                                        reference_method=reference_method, core=n_cores)
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
    if run_info.info['progress']['bg_bed_generation'] == 'No':
        if os.path.exists(bg_bed_path):
            shutil.rmtree(bg_bed_path)
        bg_map_dict = generate_background_bed(cell_feature_adata, bg_bed_path, bg_map_dict_path, 
                                              cell_number_per_group, bg_iteration, peak_confidence, n_cores)
        run_info.finish_stage('bg_bed_generation')
    else:  # run_info.info['progress']['bg_bed_generation'] == 'Finish'
        with open(bg_map_dict_path, "rb") as map_dict_file:
            bg_map_dict = pickle.load(map_dict_file)
    if run_info.info['progress']['fg_bed_generation'] == 'No':
        if os.path.exists(fg_bed_path):
            shutil.rmtree(fg_bed_path)
        fg_map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, fg_bed_path,
                                            fg_map_dict_path, fg_peaks_number_path, cell_number_per_group, peak_confidence, n_cores)
        run_info.finish_stage('fg_bed_generation')
    else:
        with open(fg_map_dict_path, "rb") as map_dict_file:
            fg_map_dict = pickle.load(map_dict_file)

    ##################################
    ### Search giggle and compute enrich score
    ##################################
    
    if reference_method in ['integration', 'both', 'chip']: # need to search chip
        search_chip = True
        fg_cell_factor_percent_df_chip, run_info = search_and_read_giggle(run_info, 'ChIP-seq', bg_bed_path, bg_chip_result_path, fg_bed_path, fg_chip_result_path, 
                                                                          fg_peaks_number_path, chip_index, genome_length, n_cores, fg_map_dict)
    if reference_method in ['integration', 'both', 'motif']: # need to search motif
        search_motif = True
        fg_cell_factor_percent_df_motif, run_info = search_and_read_giggle(run_info, 'motif', bg_bed_path, bg_motif_result_path, fg_bed_path, fg_motif_result_path, 
                                                                           fg_peaks_number_path, motif_index, genome_length, n_cores, fg_map_dict)
    ##################################
    ### Summary results
    ##################################
    if reference_method == 'integration':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_factor_percent_df_chip, 'integration', fg_cell_factor_percent_df_motif)
    elif reference_method == 'both':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_factor_percent_df_chip, 'ChIP-seq')
        regulation_adata = merge_giggle_adata(regulation_adata, fg_cell_factor_percent_df_motif, 'motif')
    elif reference_method == 'chip':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_factor_percent_df_chip, 'ChIP-seq')
    elif reference_method == 'motif':
        regulation_adata = merge_giggle_adata(processed_adata, fg_cell_factor_percent_df_motif, 'motif')
    run_info.finish_stage('enrich_adata_store')
    if store_result_adata == True:
        regulation_adata.write(result_store_path)
    ##################################
    ### Clean files
    ##################################
    if clean == True:
        try:
            shutil.rmtree(fg_bed_path)
            shutil.rmtree(bg_bed_path)
            if search_chip == True:
                shutil.rmtree(fg_chip_result_path)
                shutil.rmtree(bg_chip_result_path)
            if search_motif == True:
                shutil.rmtree(fg_motif_result_path)
                shutil.rmtree(bg_motif_result_path)
        except:
            pass
    return regulation_adata


def run( args ):
    processed_adata_path = args.processed_experiment
    feature_matrix_path = args.feature_matrix
    species = args.species
    project = args.project
    # aggregate_peak_method = args.aggregate_peak_method
    cell_number_per_group = args.cell_number_per_group
    # cell_cutoff = args.cell_cutoff
    peak_confidence = args.peak_confidence
    bg_iteration = args.bg_iter
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

    if peak_confidence == 'auto':
        if cell_number_per_group == 'auto':
            print_log('Estimating the number of cells per group...')
            cell_number_per_group = determine_number_of_cells_per_group(feature_matrix, start=1, end=70, iteration_time=30, aim_peak_number=10000, peak_confidence=None)
            print_log('Merge {cell_number} neighbor cells as one!'.format(cell_number=cell_number_per_group))
        else:
            cell_number_per_group = int(cell_number_per_group)
        peak_confidence = np.ceil(0.1*cell_number_per_group)
    else:
        peak_confidence = int(peak_confidence)
        if cell_number_per_group == 'auto':
            print_log('Estimating the number of cells per group...')
            cell_number_per_group = determine_number_of_cells_per_group(feature_matrix, start=1, end=70, iteration_time=30, aim_peak_number=10000, peak_confidence=peak_confidence)
            print_log('Merge {cell_number} neighbor cells as one!'.format(cell_number=cell_number_per_group))
        else:
            cell_number_per_group = int(cell_number_per_group)
        

    enrich(processed_adata, 
           feature_matrix, 
           project=project,
           cell_number_per_group=cell_number_per_group, 
           peak_confidence=peak_confidence, 
           bg_iteration=bg_iteration, 
           chip_index=chip_index, 
           motif_index=motif_index, 
           reference_method=reference_method, 
           store_result_adata=True,
           yes=yes, 
           clean=clean, 
           n_cores=n_cores,
           processed_adata_path=processed_adata_path,
           cell_feature_adata_path=feature_matrix_path,
           species=species
           )

    return
