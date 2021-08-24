#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   bed_generation.py
@Time    :   2021/04/16 12:34:55
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import scanpy as sc
import subprocess
import os
import sys
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED, as_completed
import pandas as pd
import random
import pickle
import anndata as ad
from SCRIPT.utilities.utils import print_log, excute_info, safe_makedirs


def generate_peak_list(cells, input_mat, peak_confidence=1):
    cell_above_cutoff_index = sc.pp.filter_genes(
        input_mat[cells, :], min_cells=peak_confidence, inplace=False)[0]
    peaks = input_mat.var_names[cell_above_cutoff_index].to_list()
    return peaks

def generate_beds(file_path, cells, input_mat, peak_confidence=1):
    peaks = generate_peak_list(cells, input_mat, peak_confidence)
    cell_barcode = os.path.basename(file_path)[:-4]# remove .bed
    if peaks.__len__() == 0:
        print_log('Warning: No peaks in {bed_path}, skip generation'.format(bed_path = file_path[:-4]))
    else:
        peaks = pd.DataFrame([p.rsplit("_", 2) for p in peaks])
        peaks.to_csv(file_path, sep="\t", header= None, index=None)
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


