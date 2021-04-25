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
from SCRIPT.enrichment.bed_generation import generate_background_bed, generate_neighbor_bed, generate_cluster_bed
from SCRIPT.enrichment.validation import check_para
from SCRIPT.enrichment.utils import EnrichRunInfo, time_estimate
from SCRIPT.enrichment.post_processing import extract_by_cell_cluster, map_factor_on_ChIP, merge_giggle_adata
from SCRIPT.enrichment.calculation import cal_rank_table_batch, cal_p_table_batch
from SCRIPT.enrichment.search_giggle import search_giggle_batch, read_giggle_result_batch
from SCRIPT.utilities.utils import read_config, read_SingleCellExperiment_rds, print_log, store_to_pickle, read_pickle, safe_makedirs
from SCRIPT.Constants import *

def search_and_read_giggle(run_info, tp, bg_bed_path, bg_result_path, fg_bed_path, fg_result_path, index, n_cores, fg_map_dict, aggregate_peak_method):
    # tp(type) is 'ChIP-seq' or 'motif'
    if tp == 'ChIP-seq':
        if run_info.info['progress']['bg_bed_chip_search'] == 'No':
            search_giggle_batch(bg_bed_path, bg_result_path, index, n_cores, tp)
            run_info.finish_stage('bg_bed_chip_search')
        if run_info.info['progress']['fg_bed_chip_search'] == 'No':
            search_giggle_batch(fg_bed_path, fg_result_path, index, n_cores, tp)
            run_info.finish_stage('fg_bed_chip_search')
    else:
        if run_info.info['progress']['bg_bed_motif_search'] == 'No':
            search_giggle_batch(bg_bed_path, bg_result_path, index, n_cores, tp)
            run_info.finish_stage('bg_bed_motif_search')
        if run_info.info['progress']['fg_bed_motif_search'] == 'No':
            search_giggle_batch(fg_bed_path, fg_result_path, index, n_cores, tp)
            run_info.finish_stage('fg_bed_motif_search')

    bg_dataset_cell_raw_score_df = read_giggle_result_batch(bg_result_path, n_cores, 'background {tp}'.format(tp=tp))
    
    
    if tp == 'ChIP-seq':
        bg_raw_score_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'bg_files', 'bg_dataset_cell_raw_score_df_ChIP.pk')
        if run_info.info['progress']['bg_dataset_cell_raw_score_chip_df_store'] == 'No':
            store_to_pickle(bg_dataset_cell_raw_score_df, bg_raw_score_path)
            run_info.finish_stage('bg_dataset_cell_raw_score_chip_df_store')
        else:
            bg_dataset_cell_raw_score_df = read_pickle(bg_raw_score_path)

        fg_raw_score_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_cell_raw_score_df_ChIP.pk')
        fg_percent_score_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_cell_percent_df_ChIP.pk')

        if run_info.info['progress']['fg_dataset_cell_raw_score_chip_df_store'] == 'No':
            fg_dataset_cell_raw_score_df = read_giggle_result_batch(fg_result_path, n_cores, 'foreground {tp}'.format(tp=tp))
            store_to_pickle(fg_dataset_cell_raw_score_df, fg_raw_score_path)
            run_info.finish_stage('fg_dataset_cell_raw_score_chip_df_store')
        else:
            fg_dataset_cell_raw_score_df = read_pickle(fg_raw_score_path)

        if run_info.info['progress']['fg_dataset_cell_percent_chip_df_store'] == 'No':
            fg_dataset_cell_percent_df = cal_p_table_batch(fg_dataset_cell_raw_score_df, bg_dataset_cell_raw_score_df, n_cores)
            store_to_pickle(fg_dataset_cell_percent_df, fg_percent_score_path)
            run_info.finish_stage('fg_dataset_cell_percent_chip_df_store')
        else:
            fg_dataset_cell_percent_df = read_pickle(fg_percent_score_path)
    else:
        bg_raw_score_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'bg_files', 'bg_dataset_cell_raw_score_df_motif.pk')
        if run_info.info['progress']['bg_dataset_cell_raw_score_motif_df_store'] == 'No':
            store_to_pickle(bg_dataset_cell_raw_score_df, bg_raw_score_path)
            run_info.finish_stage('bg_dataset_cell_raw_score_motif_df_store')
        else:
            bg_dataset_cell_raw_score_df = read_pickle(bg_raw_score_path)

        fg_raw_score_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_cell_raw_score_df_motif.pk')
        fg_percent_score_path = os.path.join(run_info.info['project_folder'], 'enrichment', 'fg_files', 'fg_dataset_cell_percent_df_motif.pk')
        if run_info.info['progress']['fg_dataset_cell_raw_score_motif_df_store'] == 'No':
            fg_dataset_cell_raw_score_df = read_giggle_result_batch(fg_result_path, n_cores, 'foreground {tp}'.format(tp=tp))
            store_to_pickle(fg_dataset_cell_raw_score_df, fg_raw_score_path)
            run_info.finish_stage('fg_dataset_cell_raw_score_motif_df_store')
        else:
            fg_dataset_cell_raw_score_df = read_pickle(fg_raw_score_path)
            
        if run_info.info['progress']['fg_dataset_cell_percent_motif_df_store'] == 'No':
            fg_dataset_cell_percent_df = cal_p_table_batch(fg_dataset_cell_raw_score_df, bg_dataset_cell_raw_score_df, n_cores)
            store_to_pickle(fg_dataset_cell_percent_df, fg_percent_score_path)
            run_info.finish_stage('fg_dataset_cell_percent_motif_df_store')
        else:
            fg_dataset_cell_percent_df = read_pickle(fg_percent_score_path)
    

    if aggregate_peak_method == "group":
        fg_dataset_cell_percent_df = extract_by_cell_cluster(fg_dataset_cell_percent_df.copy(), fg_map_dict)
    # transpose is used to better merge table to h5ad (anndata.obs's row is cell, col is variable)
    if tp == 'ChIP-seq':
        fg_cell_factor_percent_df = map_factor_on_ChIP(fg_dataset_cell_percent_df).T
    else:
        fg_cell_factor_percent_df = fg_dataset_cell_percent_df.T # motif's dataset is same to factor
    return fg_cell_factor_percent_df, run_info


def enrich(processed_adata, cell_feature_adata, project='',
           aggregate_peak_method='group', cell_number_per_group=50, cell_cutoff=20, peak_confidence=5, bg_iteration='auto', 
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
        'aggregate_peak_method':aggregate_peak_method,
        'cell_number_per_group':cell_number_per_group,
        'cell_cutoff':cell_cutoff,
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
    fg_motif_result_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_motif_result')
    fg_chip_result_path = os.path.join(project, 'enrichment', 'fg_files', 'fg_chip_result')
    bg_bed_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_bed')
    bg_map_dict_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_bed.pk')
    bg_chip_result_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_chip_result') 
    bg_motif_result_path = os.path.join(project, 'enrichment', 'bg_files', 'bg_motif_result')
    result_store_path = os.path.join(project, 'enrichment', 'SCRIPT_enrichment.h5ad')
    safe_makedirs(fg_bed_path)
    safe_makedirs(bg_bed_path)
    ##################################
    ### pre-check
    ##################################
    print_log('Checking parameters ...')
    if bg_iteration == "auto":
        bg_iteration = int(processed_adata.shape[0] * 5 / cell_number_per_group) + 1
    check_para(processed_adata, cell_feature_adata, project,
                aggregate_peak_method, cell_number_per_group, cell_cutoff, peak_confidence, bg_iteration, 
                chip_index, motif_index, reference_method, result_store_path,
                yes, clean, n_cores)
    print_log('Estimating running time ...')
    elapse, future_time = time_estimate(cell_number = cell_feature_adata.shape[0], bg_iteration_number=bg_iteration, 
                                        peak_methods=aggregate_peak_method, cell_number_per_group=cell_number_per_group, 
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
    # generate background peak, if length same as iteration, we consider it has estimated, skip generation.
    # if user generate same length background, but diff depth, may report unaccurate result.
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
        if aggregate_peak_method == "group":
            fg_map_dict = generate_cluster_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, cell_number_per_group, cell_cutoff, peak_confidence)
            run_info.finish_stage('fg_bed_generation')
        if aggregate_peak_method == "nearest":
            fg_map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, cell_number_per_group, peak_confidence, n_cores)
            run_info.finish_stage('fg_bed_generation')
    else:
        with open(fg_map_dict_path, "rb") as map_dict_file:
            fg_map_dict = pickle.load(map_dict_file)

    ##################################
    ### Search giggle and compute enrich score
    ##################################
    if reference_method in ['integration', 'both', 'chip']: # need to search chip
        fg_cell_factor_percent_df_chip, run_info = search_and_read_giggle(run_info, 'ChIP-seq', bg_bed_path, bg_chip_result_path, fg_bed_path, fg_chip_result_path, 
                                                                          chip_index, n_cores, fg_map_dict, aggregate_peak_method)
    if reference_method in ['integration', 'both', 'motif']: # need to search motif
        fg_cell_factor_percent_df_motif, run_info = search_and_read_giggle(run_info, 'motif', bg_bed_path, bg_motif_result_path, fg_bed_path, fg_motif_result_path, 
                                                                           motif_index, n_cores, fg_map_dict, aggregate_peak_method)
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
    # if search_chip == True and search_motif == True:
    #     if integration == True:
    #         regulation_adata = merge_giggle_adata(processed_adata, chip_result_p, 'integration', motif_result_p)
    #     else:
    #         regulation_adata = merge_giggle_adata(processed_adata, chip_result_p, 'ChIP-seq')
    #         regulation_adata = merge_giggle_adata(regulation_adata, motif_result_p, 'motif')
    # elif search_chip == True and search_motif != True:
    #     regulation_adata = merge_giggle_adata(processed_adata, chip_result_p, 'ChIP-seq')
    # else:
    #     regulation_adata = merge_giggle_adata(processed_adata, motif_result_p, 'motif')
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
    aggregate_peak_method = args.aggregate_peak_method
    cell_number_per_group = args.cell_number_per_group
    cell_cutoff = args.cell_cutoff
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

    processed_adata = ad.read_h5ad(processed_adata_path)
    feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)

    enrich(processed_adata, 
           feature_matrix, 
           project=project,
           aggregate_peak_method=aggregate_peak_method, 
           cell_number_per_group=cell_number_per_group, 
           cell_cutoff=cell_cutoff, 
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