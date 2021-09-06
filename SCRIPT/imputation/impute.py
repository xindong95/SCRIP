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
import pandas as pd
import scanpy as sc
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED, as_completed
from multiprocessing import Process, Pool
from SCRIPT.enrichment.search_seqpare import read_seqpare_result
from SCRIPT.utilities.utils import print_log, safe_makedirs, excute_info

def search_seqpare_factor(bed_path, result_path, index_path, factor):
    cmd = 'seqpare "{index_path}/{factor}*.bed.gz" "{bed_path}" -m 1 -o {result_path}\n'.format(
        index_path=index_path, result_path=result_path, bed_path=bed_path, factor=factor)
    subprocess.run(cmd, shell=True, check=True)


def search_seqpare_factor_batch(bed_folder, result_folder, index_path, factor, n_cores=8):
    print_log('Start searching beds from {factor} index ...'.format(factor=factor))
    safe_makedirs(result_folder)
    beds = os.listdir(bed_folder)
    args = []
    for bed in beds:
        barcodes = bed[:-7]  # remove suffix '.bed.gz'
        args.append((os.path.join(bed_folder, bed),
                     os.path.join(result_folder, barcodes + '.txt'),
                     index_path, factor))
    with Pool(n_cores) as p:
        p.starmap(search_seqpare_factor, args)
    print_log('Finished searching beds from {factor} index ...'.format(factor=factor))


def generate_peak_list(cells, input_mat, peak_confidence=1):
    cell_above_cutoff_index = sc.pp.filter_genes(
        input_mat[cells, :], min_cells=peak_confidence, inplace=False)[0]
    peaks = input_mat.var_names[cell_above_cutoff_index].to_list()
    return peaks


def generate_beds(file_path, cells, input_mat, peak_confidence=1):
    peaks = generate_peak_list(cells, input_mat, peak_confidence)
    cell_barcode = os.path.basename(file_path)[:-4]  # remove .bed
    if peaks.__len__() == 0:
        print_log('Warning: No peaks in {bed_path}, skip generation'.format(bed_path=file_path[:-4]))
    else:
        peaks = pd.DataFrame([p.rsplit("_", 2) for p in peaks])
        peaks.to_csv(file_path, sep="\t", header=None, index=None)
        cmd = 'sort --buffer-size 2G -k1,1 -k2,2n -k3,3n {bed_path} | bgzip -c > {bed_path}.gz\n'.format(bed_path=file_path)
        cmd += 'rm {bed_path}'.format(bed_path=file_path)
        subprocess.run(cmd, shell=True, check=True)
    return [cell_barcode, peaks.__len__()]


@excute_info('Start generating cells beds ...', 'Finished generating cells beds!')
def generate_beds_by_matrix(cell_feature_adata, beds_path, peaks_number_path, n_cores):
    safe_makedirs(beds_path)
    # total_cnt = adata.obs.index.__len__()
    executor = ThreadPoolExecutor(max_workers=n_cores)
    all_task = []
    for cell in cell_feature_adata.obs.index:
        # neighbor_cells = find_nearest_cells(cell, coor_table, n_neighbor, step)
        # map_dict[cell] = neighbor_cells
        all_task.append(executor.submit(generate_beds, beds_path + "/" + str(cell) + ".bed", cell, cell_feature_adata))
    wait(all_task, return_when=ALL_COMPLETED)
    pd.DataFrame([_.result() for _ in as_completed(all_task)]).to_csv(peaks_number_path, header=None, index=None, sep='\t')
    return


def get_factor_affinity(input_mat, bed_file_path, reference, factor, ccre_number):
    peaks = input_mat.var_names.to_list()
    peaks_number = peaks.__len__()
    ref_number = pd.read_csv(reference + '/peaks_number.txt', sep='\t', index_col=0, header=None)
    factor_idx = [i for i in ref_number.index if i.startswith(factor)]
    ref_number = ref_number.loc[factor_idx,:]
    
    peaks = pd.DataFrame([p.rsplit("_", 2) for p in peaks])
    peaks.to_csv(bed_file_path, sep="\t", header=None, index=None)
    cmd = 'sort --buffer-size 2G -k1,1 -k2,2n -k3,3n {bed_path} | bgzip -c > {bed_path}.gz\n'.format(bed_path=bed_file_path)
    cmd += 'rm {bed_path}'.format(bed_path=bed_file_path)
    subprocess.run(cmd, shell=True, check=True)

    search_seqpare_factor(bed_file_path + '.gz', bed_file_path[0:-4] + '.txt', reference, factor)
    all_peak_result = read_seqpare_result([bed_file_path[0:-4] + '.txt'])
    affinity = all_peak_result.iloc[:,0]/(ref_number[1]*peaks_number/ccre_number)
    return affinity

def impute():
    return

def run(args):
    return 
