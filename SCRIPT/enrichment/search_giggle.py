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
import sys
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED
from SCRIPT.utilities.utils import print_log, excute_info, safe_makedirs
from multiprocessing import Process, Pool

def search_giggle(bed_path, result_path, index_path, genome_length):
#     bed = bed_path.split("/")[-1]
#     if bed_path.endswith('gz'):
    cmd = 'giggle search -i "{index_path}" -s -q "{bed_path}" -g {genome_length}> "{result_path}"\n'.format(
        index_path=index_path, result_path=result_path, bed_path=bed_path, genome_length=genome_length)
#     print(cmd)
#     else:
#         cmd = 'sort --buffer-size 2G -k1,1 -k2,2n -k3,3n {bed_path} | bgzip -c > {bed_path}.gz\n'.format(bed_path=bed_path)
#         cmd += 'giggle search -i {index_path} -s -q {bed_path}.gz > {result_path}\n'.format(index_path=index_path, result_path=result_path, bed_path=bed_path)
#         cmd += 'rm {bed_path}'.format(bed_path=bed_path)
    subprocess.run(cmd, shell=True, check=True)


def search_giggle_batch(bed_folder, result_folder, index_path, genome_length, n_cores=8, tp=''):
    print_log('Start searching foreground beds from {tp} index ...'.format(tp=tp))
    safe_makedirs(result_folder)
    beds = os.listdir(bed_folder)
    args = []
    for bed in beds:
        barcodes = bed[:-7] # remove suffix '.bed.gz'
        args.append((os.path.join(bed_folder, bed), 
                     os.path.join(result_folder, barcodes + '.txt'), 
                     index_path, 
                     genome_length))
    with Pool(n_cores) as p:
        p.starmap(search_giggle, args)
    print_log('Finished searching foreground beds from {tp} index ...'.format(tp=tp))

# def read_giggle_result(path):
#     """For giggle stored path, return a table, col is cell cluster / cell and row is factor"""
#     file_list = os.listdir(path)
#     print_log('Reading search result in {path}, in total {number} files.'.format(path=path, number=len(file_list)))
#     for i in range(len(file_list)):
#         if i%50 == 0:
#             print_log("finished {percentage:.2f} %".format(percentage = i*100/file_list.__len__()), end="\r")
#         giggle_result = file_list[i]
#         cell_bc = giggle_result[:-4] # remove suffix '.txt'
#         dtframe = pd.read_csv(os.path.join(path, giggle_result), sep="\t", index_col=False)
#         if i == 0:
#             dtframe = dtframe.loc[:,["#file", "combo_score"]]
#             total = dtframe.rename(columns={'combo_score':cell_bc}).copy()
#         else:
#             newcol = dtframe[["combo_score"]]
#             total[cell_bc] = newcol
#     idList = [i.replace(".bed.gz","") for i in total['#file']]
#     total = total.rename(columns={"#file":"id"})
#     total["id"] = idList
#     total = total.set_index("id")
#     return total

def read_giggle_result(files, i):
#     giggle_result = os.path.basename(files[0])
#     cell_bc = giggle_result[:-4] # remove suffix '.txt'
#     dtframe = pd.read_csv(files[0], sep="\t", index_col=False)
#     dtframe = dtframe.loc[:,["#file", "combo_score"]]
#     total = dtframe.rename(columns={'combo_score':cell_bc}).copy()
    for i in range(len(files)):
        giggle_result = os.path.basename(files[i])
        cell_bc = giggle_result[:-4] # remove suffix '.txt'
        dtframe = pd.read_csv(files[i], sep="\t", index_col=False)
        if i == 0:
            dtframe = dtframe.loc[:,["#file", "combo_score"]]
            dataset_cell_score_df = dtframe.rename(columns={'combo_score':cell_bc}).copy()
        else:
            newcol = dtframe[["combo_score"]]
            dataset_cell_score_df[cell_bc] = newcol
    idList = [i[:-7] for i in dataset_cell_score_df['#file']] # remove suffix '.bed.gz'
    dataset_cell_score_df = dataset_cell_score_df.rename(columns={"#file":"id"})
    dataset_cell_score_df["id"] = idList
    dataset_cell_score_df = dataset_cell_score_df.set_index("id")
    return dataset_cell_score_df

def read_giggle_result_batch(path, n_cores=8, tp=''):
    print_log("Reading searching results, using {n} cores...".format(n=n_cores))
    file_list = os.listdir(path)
    result_split = np.array_split(file_list, n_cores)
    args = [[[os.path.join(path, j) for j in list_chunk], i] for (i, list_chunk) in enumerate(result_split)]
    with Pool(n_cores) as p:
        result = p.starmap(read_giggle_result, args)
    dataset_cell_score_df = pd.concat([i for i in result], axis=1)
    print_log("Finished reading {tp} index search result!".format(tp=tp))
    return dataset_cell_score_df

def read_giggle_result_fisher(files, i):
    for i in range(len(files)):
        giggle_result = os.path.basename(files[i])
        cell_bc = giggle_result[:-4]  # remove suffix '.txt'
        dtframe = pd.read_csv(files[i], sep="\t", index_col=False)
        if i == 0:
            dtframe = dtframe.loc[:, ["#file", "fishers_two_tail"]]
            total = dtframe.rename(
                columns={'fishers_two_tail': cell_bc}).copy()
        else:
            newcol = dtframe[["fishers_two_tail"]]
            total[cell_bc] = newcol
    idList = [i[:-7] for i in total['#file']]  # remove suffix '.bed.gz'
    total = total.rename(columns={"#file": "id"})
    total["id"] = idList
    total = total.set_index("id")
    return total

def read_giggle_result_fisher_batch(path, n_cores=8, tp=''):
    print_log(
        "Reading searching results of fisher's test, using {n} cores...".format(n=n_cores))
    file_list = os.listdir(path)
    result_split = np.array_split(file_list, n_cores)
    args = [[[os.path.join(path, j) for j in list_chunk], i]
            for (i, list_chunk) in enumerate(result_split)]
    with Pool(n_cores) as p:
        result = p.starmap(read_giggle_result_fisher, args)
    total = pd.concat([i for i in result], axis=1)
    return total

def read_giggle_result_odd(files, i):
    for i in range(len(files)):
        giggle_result = os.path.basename(files[i])
        cell_bc = giggle_result[:-4]  # remove suffix '.txt'
        dtframe = pd.read_csv(files[i], sep="\t", index_col=False)
        if i == 0:
            dtframe = dtframe.loc[:, ["#file", "odds_ratio"]]
            total = dtframe.rename(columns={'odds_ratio': cell_bc}).copy()
        else:
            newcol = dtframe[["odds_ratio"]]
            total[cell_bc] = newcol
    idList = [i[:-7] for i in total['#file']]  # remove suffix '.bed.gz'
    total = total.rename(columns={"#file": "id"})
    total["id"] = idList
    total = total.set_index("id")
    return total

def read_giggle_result_odds_batch(path, n_cores=8, tp=''):
    print_log(
        "Reading searching results of odds ratio, using {n} cores...".format(n=n_cores))
    file_list = os.listdir(path)
    result_split = np.array_split(file_list, n_cores)
    args = [[[os.path.join(path, j) for j in list_chunk], i]
            for (i, list_chunk) in enumerate(result_split)]
    with Pool(n_cores) as p:
        result = p.starmap(read_giggle_result_odd, args)
    total = pd.concat([i for i in result], axis=1)
    return total
