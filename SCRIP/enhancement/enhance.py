#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   enhance.py
@Time    :   2021/08/24 10:40:23
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''


from sklearn.neighbors import BallTree
from sklearn.neighbors import KDTree
import sys
from SCRIP.utilities.utils import print_log, excute_info, safe_makedirs, store_to_pickle
# from SCRIP.enrichment.bed_generation import generate_beds
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED, as_completed
import scanpy as sc
import anndata as ad
import random
import scipy
import numpy as np
import pandas as pd
from SCRIP.enrichment.bed_generation import generate_peak_list
from multiprocessing import Pool


def determine_number_of_cells_per_group(input_mat, start=1, end=70, iteration_time=30, aim_peak_number=10000, peak_confidence=None):
    peak_number = 0
    cell_pool = input_mat.obs_names.to_list()
    cell_number = start - 1
    while peak_number <= aim_peak_number:
        cell_number += 1
        if peak_confidence == None:
            peak_confidence = np.ceil(0.2*cell_number)
        tmp_peak_number = []
        for i in range(iteration_time):
            random.seed(i)
            cells = random.sample(cell_pool, cell_number)
            tmp_peak_number.append(generate_peak_list(cells, input_mat, peak_confidence).__len__())
        peak_number = sum(tmp_peak_number)/iteration_time
        if cell_number > end:
            print_log('Can not find enough cells!')
            sys.exit()
    return cell_number


def sub_coor_table_in_small_square(cell, coor_table, step, t):
    tmp = coor_table[(coor_table['X'] < coor_table.loc[cell, 'X'] + step * t) &
                     (coor_table['X'] > coor_table.loc[cell, 'X'] - step * t) &
                     (coor_table['Y'] < coor_table.loc[cell, 'Y'] + step * t) &
                     (coor_table['Y'] > coor_table.loc[cell, 'Y'] - step * t)]
    return tmp.copy()


def find_nearest_cells_umap(cell, coor_table, n_neighbor=20, step=None):
    #   coor_table has two columns: X, Y, the index of coortable is cell barcodes
    if step == None:
        up_limit = coor_table.max()
        down_limit = coor_table.min()
        width = up_limit.X - down_limit.X
        height = up_limit.Y - down_limit.Y
        step = min(width, height)/200
    t = 1
    tmp = pd.DataFrame()
    while tmp.shape[0] < n_neighbor:
        tmp = sub_coor_table_in_small_square(cell, coor_table, step, t)
        t += 1
    squ_distance = [(tmp.loc[i, "X"] - tmp.loc[cell, "X"])**2 + (tmp.loc[i, "Y"] - tmp.loc[cell, "Y"])**2 for i in tmp.index]
    tmp['distance'] = squ_distance
    tmp = tmp.sort_values(by='distance', ascending=True)
    # since we use square to check the distance, check whether the point is out the circle
    while tmp.iloc[n_neighbor-1, 2] > ((step*t)**2):
        t += 1
        tmp = sub_coor_table_in_small_square(cell, coor_table, step, t)
        squ_distance = [(tmp.loc[i, "X"] - tmp.loc[cell, "X"])**2 + (tmp.loc[i, "Y"] - tmp.loc[cell, "Y"])**2 for i in tmp.index]
        tmp['distance'] = squ_distance
        tmp = tmp.sort_values(by='distance', ascending=True)
    neighbor_bcs = tmp.index[0:n_neighbor].tolist()
    return neighbor_bcs

def cal_neighbor_cell_peak_mat_umap(table, input_mat, coor_table, impute_n, i):
    # print(i)
    for cell_barcode in table.index:
        nearest_bc = find_nearest_cells_umap(cell_barcode, coor_table, n_neighbor=impute_n)
#         table.loc[cell_barcode,:] = [1 if i else 0 for i in input_mat.loc[nearest_bc,:].sum() >= 1]
        table.loc[cell_barcode, :] = input_mat.loc[nearest_bc, :].sum()
    return table

def cal_neighbor_cell_peak_mat_umap_batch(input_mat, coor_table, impute_n, n_cores=8):
    print_log("Generating neighbor cells peak matrix, divide into {n} chunks...".format(n=n_cores))
    input_table_split = np.array_split(input_mat, n_cores)
    args = [[table, input_mat, coor_table, impute_n, i] for (i, table) in enumerate(input_table_split)]
    with Pool(n_cores) as p:
        result = p.starmap(cal_neighbor_cell_peak_mat_umap, args)
    cell_peak = pd.concat([i for i in result])
    print_log('Finished!')
#     return scipy.sparse.csr_matrix(cell_peak)
    return cell_peak


def find_nearest_cells(q_point, tree, n_neighbor):
    _, ind = tree.query(q_point, k=n_neighbor+1)
    return ind[0][1:]


def cal_neighbor_cell_peak_mat(sub_mat, input_mat, tree, pc_table, impute_n, start_idx, i):
    end_index = start_idx + sub_mat.shape[0]
    k = 0
    for idx in range(start_idx, end_index):
        nearest_bc_idx = find_nearest_cells(np.reshape(pc_table[idx, :], (1, -1)), tree, n_neighbor=impute_n)
#         scipy.sparse.csr_matrix(input_mat[nearest_bc_idx,:].sum(0))
        sub_mat[k, :] = input_mat[nearest_bc_idx, :].sum(0)
        k += 1
    return sub_mat


def cal_neighbor_cell_peak_mat_batch(input_mat, impute_n=5, KD_leafsize=80, nPC=50, n_cores=8):
    '''
    input_mat:
    a csr sparse matrix, which can be get by adata.X
    
    '''
    print_log("Building KD tree...")
    pc_table = sc.tl.pca(input_mat, n_comps=50, svd_solver='arpack')
    tree = BallTree(pc_table, KD_leafsize)
    print_log("Calculating neighbors, divide into {n} chunks...".format(n=n_cores))
    cell_number = input_mat.shape[0]
    index_split = [i for i in range(0, cell_number, int(cell_number/n_cores))] + [cell_number]
#     input_table_split = np.array_split(input_mat_dense, n_cores)
    input_mat_lil = input_mat.tolil()
    input_mat_split = [input_mat_lil[index_split[i]:index_split[i+1], :] for i in range(index_split.__len__()-1)]
    args = [[sub_mat, input_mat, tree, pc_table, impute_n, index_split[i], i] for (i, sub_mat) in enumerate(input_mat_split)]
#     print(args)
    with Pool(n_cores) as p:
        result = p.starmap(cal_neighbor_cell_peak_mat, args)
    cell_peak_csr = scipy.sparse.vstack(result).tocsr()
    print_log('Finished!')
    return cell_peak_csr

# @excute_info('Start generating nearest neighbor cells beds ...', 'Finished generating nearest neighbor cells beds!')
# def generate_neighbor_bed_umap(adata, input_mat, bed_path, map_dict_store_path, peaks_number_store_path, n_neighbor=10, peak_confidence=2, n_cores=8):
#     coor_table = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index, columns=["X", "Y"])
#     width = coor_table.max().X - coor_table.min().X
#     height = coor_table.max().Y - coor_table.min().Y
#     step = min(width, height)/200
#     map_dict = {}
#     safe_makedirs(bed_path)
#     # total_cnt = adata.obs.index.__len__()
#     executor = ThreadPoolExecutor(max_workers=n_cores)
#     all_task = []
#     for cell in adata.obs.index:
#         neighbor_cells = find_nearest_cells_umap(cell, coor_table, n_neighbor, step)
#         map_dict[cell] = neighbor_cells
#         all_task.append(executor.submit(generate_beds, bed_path + "/" + str(cell) + ".bed", neighbor_cells, input_mat, peak_confidence))
#     wait(all_task, return_when=ALL_COMPLETED)
#     pd.DataFrame([_.result() for _ in as_completed(all_task)]).to_csv(peaks_number_store_path, header=None, index=None, sep='\t')
#     with open(map_dict_store_path, "wb") as map_dict_file:
#         pickle.dump(map_dict, map_dict_file)
#     return map_dict


def enhance(input_mat, impute_n=5, KD_leafsize=80, nPC=50, path='SCRIP/enhancement/', binarize=True, n_cores=8):
    '''
    input_mat:
    a csr sparse matrix
    
    '''
    safe_makedirs(path)
    imputed_csr = cal_neighbor_cell_peak_mat_batch(input_mat, impute_n=impute_n, KD_leafsize=KD_leafsize, nPC=nPC, n_cores=n_cores)
    if binarize == True:
        imputed_csr[imputed_csr > 1] = 1
    store_to_pickle(imputed_csr, path + 'imputed.csr.pk')
    return imputed_csr

def run(args):
    processed_adata_path = args.processed_experiment
    feature_matrix_path = args.feature_matrix
    project = args.project
    cell_number_impute = args.cell_number_impute
    peak_confidence_impute = args.peak_confidence_impute
    yes = args.yes
    clean = args.clean
    n_cores = args.n_cores

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
