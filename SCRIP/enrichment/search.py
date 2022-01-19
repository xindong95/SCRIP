#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   search_giggle.py
@Time    :   2021/04/16 12:34:35
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''


import subprocess
import os
from multiprocessing import Pool
import warnings
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED
import pandas as pd
import numpy as np
from SCRIP.utilities.utils import print_log, safe_makedirs

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def search_ref(bed_path, result_path, index_path):
    cmd = f'giggle search -i "{index_path}" -q "{bed_path}" -s > "{result_path}"\n'
    # cmd = f'igd search {index_path}/ref.igd -q {bed_path} | head -n -1 | cut -f 2,3,4 > {result_path}'
    # cmd = f'seqpare "{index_path}/*.bed.gz" "{bed_path}" -m 1 -o {result_path}\n'
    subprocess.run(cmd, shell=True, check=True)

def search_ref_batch(bed_folder, result_folder, index_path, n_cores=8, tp=''):
    print_log(f'Start searching beds from {tp} index ...')
    safe_makedirs(result_folder)
    beds = os.listdir(bed_folder)
    args = []
    for bed in beds:
        barcodes = bed[:-7]  # remove suffix '.bed.gz'
        args.append((os.path.join(bed_folder, bed),
                     os.path.join(result_folder, barcodes + '.txt'),
                     index_path))
    with Pool(n_cores) as p:
        p.starmap(search_ref, args)
    print_log(f'Finished searching beds from {tp} index ...')


def read_search_result(files):
    for i in enumerate(files):
        result_name = os.path.basename(i[1])
        cell_bc = result_name[:-4]  # remove suffix '.txt'
        dtframe = pd.read_csv(i[1], sep="\t", index_col=0, comment='#', header=None)
        read_col = 2  # 1 file_size 2 overlaps 3 odds_ratio 4 fishers_two_tail 5 fishers_left_tail 6 fishers_right_tail 7 combo_score
        if i[0] == 0:
            dtframe = dtframe.loc[:, [read_col]].copy()
            dataset_cell_score_df = dtframe.rename(columns={read_col: cell_bc}).copy()
        else:
            dataset_cell_score_df[cell_bc] = dtframe.loc[:, read_col]
    dataset_cell_score_df.index = [n.rsplit('/', 1)[0][:-7] for n in dataset_cell_score_df.index]  # remove suffix '.bed.gz'
    return dataset_cell_score_df


def read_search_result_batch(path, n_cores=8, tp=''):
    print_log(f"Reading searching results, using {n_cores} cores...")
    file_list = os.listdir(path)
    result_split = np.array_split(file_list, n_cores)
    args = [[[os.path.join(path, j) for j in list_chunk]] for list_chunk in result_split]
    with Pool(n_cores) as p:
        result = p.starmap(read_search_result, args)
    dataset_cell_score_df = pd.concat([i for i in result], axis=1)
    print_log(f"Finished reading {tp} index search result!")
    return dataset_cell_score_df
