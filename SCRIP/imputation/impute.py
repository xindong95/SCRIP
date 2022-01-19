#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   impute.py
@Time    :   2021/08/25 19:09:37
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''
import os
import subprocess
import random
import time
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED, as_completed
from multiprocessing import Pool
import pandas as pd
import scipy
import pybedtools
import scanpy as sc
from SCRIP.enrichment.search import read_search_result_batch
from SCRIP.enrichment.bed_generation import generate_beds_by_matrix
from SCRIP.utilities.utils import print_log, safe_makedirs, excute_info, write_to_mtx, store_to_pickle, read_config




def search_ref_factor(bed_path, result_path, index_path, factor):
    cmd = f'giggle search -i "{index_path}" -q "{bed_path}" -s -f {factor}_ > "{result_path}"\n'
    # cmd = f'igd search {index_path}/ref.igd -q {bed_path} | head -n -1 | cut -f 2,3,4 > {result_path}'
    # cmd = f'seqpare "{index_path}/*.bed.gz" "{bed_path}" -m 1 -o {result_path}\n'
    subprocess.run(cmd, shell=True, check=True)

def search_ref_factor_batch(bed_folder, result_folder, index_path, factor, n_cores=8, tp=''):
    print_log(f'Start searching beds from {tp} index ...')
    safe_makedirs(result_folder)
    beds = os.listdir(bed_folder)
    args = []
    for bed in beds:
        barcodes = bed[:-7]  # remove suffix '.bed.gz'
        args.append((os.path.join(bed_folder, bed),
                     os.path.join(result_folder, barcodes + '.txt'),
                     index_path,
                     factor))
    with Pool(n_cores) as p:
        p.starmap(search_ref_factor, args)
    print_log(f'Finished searching beds from {tp} index ...')

# def get_factor_affinity(input_mat, bed_file_path, reference, factor, ccre_number):
#     peaks = input_mat.var_names.to_list()
#     peaks_number = peaks.__len__()
#     ref_number = pd.read_csv(reference + '/peaks_number.txt', sep='\t', index_col=0, header=None)
#     factor_idx = [i for i in ref_number.index if i.startswith(factor)]
#     ref_number = ref_number.loc[factor_idx,:]
#     peaks = pd.DataFrame([p.rsplit("_", 2) for p in peaks])
#     peaks.to_csv(bed_file_path, sep="\t", header=None, index=None)
#     cmd = 'sort --buffer-size 2G -k1,1 -k2,2n -k3,3n {bed_path} | bgzip -c > {bed_path}.gz\n'.format(bed_path=bed_file_path)
#     cmd += 'rm {bed_path}'.format(bed_path=bed_file_path)
#     subprocess.run(cmd, shell=True, check=True)

#     search_seqpare_factor(bed_file_path + '.gz', bed_file_path[0:-4] + '.txt', reference, factor)
#     all_peak_result = read_seqpare_result([bed_file_path[0:-4] + '.txt'])
#     affinity = all_peak_result.iloc[:,0]/(ref_number[1]*peaks_number/ccre_number)
#     return affinity


@excute_info('Getting the best reference for each cell.')
def get_factor_source(table):
    ret_table = table.copy()
    # map factor by id "_"
    factor_index_list = []
    for i in ret_table.index:
        factor_name = i.split("_")
        factor_index_list.append(factor_name[0])
    ret_table.loc[:, "Factor"] = factor_index_list
    max_index = ret_table.groupby("Factor").idxmax()
    return max_index

def impute(input_mat_adata, impute_factor, ref_path, bed_check=True, search_check=True, path='SCRIP/imputation/', write_format='', ref_baseline=500, remove_others_source=False, n_cores=8):
    '''
    impute tf ChIP data from ATAC data
    '''

    safe_makedirs(path)
    print(input_mat_adata.X.shape)
    if bed_check == True:
        if not os.path.exists(f'{path}/imputed_beds/'):
            print_log('Generating beds...')
            generate_beds_by_matrix(input_mat_adata, f'{path}/imputed_beds/', f'{path}/imputed_beds_peaks_number.txt', n_cores)
        else:
            print_log('Skip generate beds...')
    else:
        print_log('Generating beds...')
        generate_beds_by_matrix(input_mat_adata, f'{path}/imputed_beds/', f'{path}/imputed_beds_peaks_number.txt', n_cores)

    if search_check == True:
        if not os.path.exists(f'{path}/imputed_results_{impute_factor}/'):
            search_ref_factor_batch(f'{path}/imputed_beds/', f'{path}/imputed_results_{impute_factor}/', ref_path, impute_factor, n_cores)
        else:
            print_log('Skip searching beds...')
    else:
        search_ref_factor_batch(f'{path}/imputed_beds/', f'{path}/imputed_results_{impute_factor}/', ref_path, impute_factor, n_cores)

    print_log('Calculating score...')
    factor_enrich = read_search_result_batch(f'{path}/imputed_results_{impute_factor}/', n_cores)

    peaks_number = pd.read_csv(os.path.join(ref_path, 'peaks_number.txt'), sep='\t', header=None, index_col=0)
    peaks_number_factor = peaks_number.loc[[i for i in peaks_number.index if i.startswith(f'{impute_factor}_')], :].copy()
    peaks_number_baseline_index = peaks_number_factor.index[peaks_number_factor[1] > ref_baseline]
    # peaks_number_factor = peaks_number_factor.loc[peaks_number_baseline_index, :]

    factor_enrich = factor_enrich.loc[peaks_number_baseline_index, :].copy()
    factor_score = (factor_enrich.T/peaks_number_factor.loc[factor_enrich.index, 1]).T

    factor_source = get_factor_source(factor_score)
    store_to_pickle(factor_source, f'{path}/{impute_factor}_dataset_source.pk')

    chip_bed_list = [pybedtools.BedTool(os.path.join(ref_path, 'raw_beds', i + '.bed.gz')) for i in factor_source.iloc[0, :].unique()]
    if len(chip_bed_list) == 1:
        chip_bed = chip_bed_list[0]
    else:
        chip_bed = chip_bed_list[0].cat(*chip_bed_list[1:])
    data_bed = pybedtools.BedTool('\n'.join(['\t'.join(p.rsplit('_', maxsplit=2)) for p in input_mat_adata.var_names]), from_string=True)
    intersect_bed = data_bed.intersect(chip_bed, u=True)
    imputed_chip_peak = str(intersect_bed).replace('\t', '_').split('\n')[0:-1]

    chip_cell_peak = input_mat_adata[:, imputed_chip_peak].copy()
    chip_cell_peak_df = chip_cell_peak.to_df()
    if remove_others_source == True:
        for i in factor_source.iloc[0, :].unique():
            cellbc = factor_source.columns[factor_source.iloc[0, :] == i]
            tmp_dataset_bed = pybedtools.BedTool(os.path.join(ref_path, 'raw_beds', i + '.bed.gz'))
            exclude_chip_peak = str(intersect_bed.intersect(tmp_dataset_bed, v=True)).replace('\t', '_').split('\n')[0:-1]
            chip_cell_peak_df.loc[cellbc, exclude_chip_peak] = 0
    chip_cell_peak = sc.AnnData(chip_cell_peak_df)
    chip_cell_peak.X = scipy.sparse.csr.csr_matrix(chip_cell_peak.X)
    if write_format != '':
        print_log('Writing results...')
        if write_format == 'h5ad':
            chip_cell_peak.write_h5ad(f'{path}/imputed_{impute_factor}.h5ad')
        elif write_format == 'mtx':
            write_to_mtx(chip_cell_peak, f'{path}/imputed_{impute_factor}_mtx/')
    print_log('Finished!')
    return chip_cell_peak


def run_impute(args):
    feature_matrix_path = args.feature_matrix
    species = args.species
    factor = args.factor
    project = args.project
    write_format = args.format
    min_cells = args.min_cells  # for removing few features cells
    min_peaks = args.min_peaks  # for removing few cells features
    max_peaks = args.max_peaks  # for removing doublet cells
    remove_others = args.remove_others
    ref_baseline = args.ref_baseline
    n_cores = args.n_cores

    if project == '':
        tmp_chr_list = [chr(i) for i in range(ord("A"), ord("Z") + 1)] + [chr(i) for i in range(ord("a"), ord("z") + 1)] + [chr(i) for i in range(ord("0"), ord("9") + 1)]
        random.seed(time.time())
        tmp_prefix = str(time.time())[6:13].replace('.', '') + '_' + ''.join(random.sample(tmp_chr_list, 4))
        project = 'SCRIP_' + tmp_prefix

    if factor in ['H3K4me3', 'H3K4me2', 'H3K27ac', 'H3K9ac', 'H3K4me1']:
        factor_type='histone'
    else:
        factor_type='TR'

    CONFIG, _ = read_config()
    if factor_type== 'TR':
        if species == 'hs':
            chip_index = CONFIG['index']['human_tf_index']
        elif species == 'mm':
            chip_index = CONFIG['index']['mouse_tf_index']
        else:
            pass
    else:
        if species == 'hs':
            chip_index = CONFIG['index']['human_hm_index']
        elif species == 'mm':
            chip_index = CONFIG['index']['mouse_hm_index']
        else:
            pass

    if feature_matrix_path.endswith('.h5'):
        feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)
    elif feature_matrix_path.endswith('.h5ad'):
        feature_matrix = sc.read_h5ad(feature_matrix_path)

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
        sc.pp.filter_cells(feature_matrix, max_genes=max_peaks)
        # feature_matrix = feature_matrix[feature_matrix.obs.n_genes_by_counts < max_peaks, :]
    elif max_peaks == 'none':
        print_log('Skip filter by max peaks number.')
    else:
        max_peaks = int(max_peaks)
        sc.pp.filter_cells(feature_matrix, max_genes=max_peaks)
        # feature_matrix = feature_matrix[feature_matrix.obs.n_genes_by_counts < max_peaks, :]

    impute(input_mat_adata=feature_matrix,
            impute_factor=factor,
            ref_path=chip_index,
            bed_check=True,
            search_check=True,
            path=f'{project}/imputation/',
            write_format=write_format,
            ref_baseline=ref_baseline,
            remove_others_source=remove_others,
            n_cores=n_cores)
    return
