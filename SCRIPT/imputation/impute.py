#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   impute.py
@Time    :   2021/05/20 09:49:38
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import sys
from SCRIPT.utilities.utils import print_log
# import scanpy as sc
import random
import numpy as np
import pandas as pd
from SCRIPT.enrichment.bed_generation import generate_peak_list,find_nearest_cells
from multiprocessing import Pool


def cal_neighbor_cell_peak_mat(table, input_mat, coor_table, impute_n, i):
    # print(i)
    for cell_barcode in table.index:
        nearest_bc = find_nearest_cells(cell_barcode, coor_table, n_neighbor=impute_n)
#         table.loc[cell_barcode,:] = [1 if i else 0 for i in input_mat.loc[nearest_bc,:].sum() >= 1]
        table.loc[cell_barcode, :] = input_mat.loc[nearest_bc, :].sum()
    return table

def cal_neighbor_cell_peak_mat_batch(input_mat, coor_table, impute_n, n_cores=8):
    print_log("Generating neighbor cells peak matrix, divide into {n} chunks...".format(n=n_cores))
    input_table_split = np.array_split(input_mat, n_cores)
    args = [[table, input_mat, coor_table, impute_n, i] for (i, table) in enumerate(input_table_split)]
    with Pool(n_cores) as p:
        result = p.starmap(cal_neighbor_cell_peak_mat, args)
    cell_peak = pd.concat([i for i in result])
    print_log('Finished!')
#     return scipy.sparse.csr_matrix(cell_peak)
    return cell_peak

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
