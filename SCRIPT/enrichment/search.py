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
from SCRIPT.utilities.utils import print_log, safe_makedirs

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def search_seqpare(bed_path, result_path, index_path):
    cmd = f'seqpare "{index_path}/*.bed.gz" "{bed_path}" -m 1 -o {result_path}\n'
    subprocess.run(cmd, shell=True, check=True)

def search_seqpare_batch(bed_folder, result_folder, index_path, n_cores=8, tp=''):
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
        p.starmap(search_seqpare, args)
    print_log(f'Finished searching beds from {tp} index ...')


def read_seqpare_result(files):
    for i in range(len(files)):
        giggle_result = os.path.basename(files[i])
        cell_bc = giggle_result[:-4]  # remove suffix '.txt'
        dtframe = pd.read_csv(files[i], sep="\t", index_col=5)
        if i == 0:
            dtframe = dtframe.loc[:, ['teo']].copy()
            dataset_cell_score_df = dtframe.rename(columns={'teo': cell_bc})
        else:
            dataset_cell_score_df[cell_bc] = dtframe.loc[:, "teo"]
    dataset_cell_score_df.index = [i.rsplit('/', 1)[1][:-7] for i in dataset_cell_score_df.index]  # remove suffix '.bed.gz'
    return dataset_cell_score_df

def read_seqpare_result_batch(path, n_cores=8, tp=''):
    print_log(f"Reading searching results, using {n_cores} cores...")
    file_list = os.listdir(path)
    result_split = np.array_split(file_list, n_cores)
    args = [[[os.path.join(path, j) for j in list_chunk]] for list_chunk in result_split]
    with Pool(n_cores) as p:
        result = p.starmap(read_seqpare_result, args)
    dataset_cell_score_df = pd.concat([i for i in result], axis=1)
    print_log(f"Finished reading {tp} index search result!")
    return dataset_cell_score_df
