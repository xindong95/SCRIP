#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   calculation.py
@Time    :   2021/04/16 12:35:01
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import os
import sys
import numpy as np
import pandas as pd
# import scipy
from sklearn import preprocessing
from SCRIPT.utilities.utils import excute_info, print_log
# from multiprocessing import Process, Pool


def standardScaler(t):
    scaler = preprocessing.StandardScaler().fit(t)
    rt = pd.DataFrame(scaler.transform(t), index=t.index, columns=t.columns)
    return rt


@excute_info('Summary result from dataset level to factor level.')
def map_factor_on_ChIP(table):
    ret_table = table.copy()
    # map factor by id "_"
    factor_index_list = []
    for i in ret_table.index:
        factor_name = i.split("_")
        factor_index_list.append(factor_name[0])
    ret_table.loc[:, "Factor"] = factor_index_list
    factor_table = ret_table.groupby("Factor").max()
    return factor_table

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

def cal_peak_norm(ref_peak_number_path, peaks_number_path, ccre_number, affinity):
    ref_peak_number = pd.read_csv(ref_peak_number_path, sep='\t', header=None, index_col=0)
    data_peak_number = pd.read_csv(peaks_number_path, sep='\t', header=None, index_col=0)
    peak_norm_hyper = (ref_peak_number.dot(data_peak_number.T)/ccre_number)
    peak_norm_hyper = (peak_norm_hyper.T * affinity).T
    return peak_norm_hyper


def cal_score(dataset_mbm_overlap_df, dataset_bg_peak_norm_df):
    intersect_frame = dataset_mbm_overlap_df.copy()
    peak_hyper_frame = dataset_bg_peak_norm_df.copy().reindex(index=intersect_frame.index, columns=intersect_frame.columns)
    fg_dataset_cell_raw_score_df = np.log2((intersect_frame+1)/(peak_hyper_frame+1))
    return fg_dataset_cell_raw_score_df


# def zscore_normalization(data_cell_frame, by='cell'):
#     if by=="cell":
#         data_cell_frame_zscore = data_cell_frame.apply(scipy.stats.zscore, axis=0)
#     elif by == "factor":
#         data_cell_frame_zscore = data_cell_frame.apply(scipy.stats.zscore, axis=1)
#     else:
#         pass
#     return data_cell_frame_zscore


def score_normalization(data_cell_frame):
    tf_cell_table = map_factor_on_ChIP(data_cell_frame)
    tf_cell_frame_zscore = standardScaler(tf_cell_table)
    tf_cell_frame_score = (tf_cell_frame_zscore.T - tf_cell_frame_zscore.mean(1)).T
    # tf_cell_frame_zscore = tf_cell_table.apply(scipy.stats.zscore, axis=1)
    # tf_cell_frame_zscore = tf_cell_frame_zscore.apply(scipy.stats.zscore, axis=0)
    return tf_cell_frame_score
