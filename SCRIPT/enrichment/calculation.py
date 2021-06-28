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
import scipy
from SCRIPT.utilities.utils import print_log
from multiprocessing import Process, Pool

# def cal_hist_auc(arrays, bins = 500):
#     """Calculate the AUC of the arrays distribution
#     """
#     y_score = np.array(arrays)
#     hist = arrays.value_counts(bins=bins)
#     hist = hist.sort_index()
#     append_right = hist.index[bins-1].right
#     hist.index = hist.index.left
#     # bins is the regions, like [1100, 1200, 1300 ...]
#     bins = np.asarray(hist.index)
#     # bins[0] = bins[1] - np.diff(bins)[1]
#     bins = np.append(bins, append_right)
#     # values is the height of each bins, like [0, 0, 1, 2, 0, 1 ...]
#     values = np.asarray(hist.values)
#     area = sum(np.diff(bins)*values)
#     return area, bins, values

# def find_bin_idx_of_value(bins, value):
#     """Finds the bin which the value corresponds to."""
#     array = np.asarray(value)
#     idx = np.digitize(array, bins)
#     if idx == 0:
#         return 0
#     return idx-1

# def area_after_val(values, bins, val):
#     """Calculates the area of the hist after a certain value"""
#     left_bin_edge_index = find_bin_idx_of_value(bins, val)
# #     bin_width = np.diff(bins)[1]
# #     print(bin_width)
#     area = sum(np.diff(bins)[left_bin_edge_index:] * values[left_bin_edge_index:])
#     return area

# def cal_p(fg_value, bg_area, bg_bins, bg_values):
#     pvalue = area_after_val(bg_values, bg_bins, fg_value)/bg_area
#     if pvalue > 1:
#         pvalue = 1.0
#     return pvalue

# def cal_fc(fg_value, bg_mean):
#     if bg_mean * fg_value > 0: # if mean and value are opposite
#         fc = fg_value / bg_mean
#     elif bg_mean == 0:
#         fc = fg_value
#     else:
#         fc = bg_mean / (bg_mean - fg_value)
#     return fc

# calculate p value by area at the right of curve
# calculate fc by value / background average

# def cal_rank(val, bg):
#     rank = [i for i in bg if i >= val].__len__()
#     if rank == 0:
#         rank = 1
#     ret = rank/bg.__len__()
#     return ret

# def cal_rank_table(fg_table, bg_table, i):
#     print_log('chunk {i} calculating ...'.format(i=i))
#     result_table_p = fg_table.copy()
# #     result_table_fc = fg_table.copy()
#     for factor in fg_table.index:
#         factor_bg = bg_table.loc[factor,:].tolist()
#         bg_mean = np.mean(factor_bg)
#         bg_std = np.std(factor_bg)
#         if bg_mean != 0 and bg_std != 0:
#             result_table_p.loc[factor,:] = fg_table.loc[factor,:].apply(cal_rank, **{'bg':factor_bg})
#         else:
#             result_table_p.loc[factor,:] = 1
# #         result_table_fc.loc[factor,:] = fg_table.loc[factor,:].apply(cal_fc, **{'bg_mean': bg_mean})
#     print_log('chunk {i} finished calculation!'.format(i=i))
#     return result_table_p

# def cal_rank_table_batch(fg_table, bg_table, n_cores=8):
#     print_log("Calculating enrichment, divide into {n} chunks...".format(n=n_cores))
#     fg_table_split = np.array_split(fg_table, n_cores)
#     args = [[table, bg_table, i] for (i, table) in enumerate(fg_table_split)]
#     with Pool(n_cores) as p:
#         result = p.starmap(cal_rank_table, args)
#     print_log("Generating P value table ...")
#     result_table_p = pd.concat([i for i in result])
# #     print("INFO {time}, Generating FC value table ...".format(time=datetime.now()))
# #     result_table_fc = pd.concat([i[1] for i in result])
#     print_log('Finished calculation enrichment!')
#     return result_table_p


# def cal_z(score, mean, std):
#     return (score-mean)/std

# def cal_z_table(fg_table, bg_table, i):
#     print_log('chunk {i} calculating ...'.format(i=i))
#     result_table_z = fg_table.copy()
# #     result_table_fc = fg_table.copy()
#     for factor in fg_table.index:
#         factor_bg = bg_table.loc[factor,:]
#         bg_mean = np.mean(factor_bg)
#         bg_std = np.std(factor_bg)
#         if bg_mean != 0 and bg_std != 0:
#             result_table_z.loc[factor,:] = fg_table.loc[factor,:].apply(cal_z, args=(bg_mean, bg_std))
#         else:
#             result_table_z.loc[factor,:] = 0
# #         result_table_fc.loc[factor,:] = fg_table.loc[factor,:].apply(cal_fc, **{'bg_mean': bg_mean})
#     print_log('chunk {i} finished calculation!'.format(i=i))
#     return result_table_z

# def cal_p_table(fg_table, bg_table, i):
#     print_log('chunk {i} calculating ...'.format(i=i))
#     result_table_p = fg_table.copy()
# #     result_table_fc = fg_table.copy()
#     for factor in fg_table.index:
#         factor_bg = bg_table.loc[factor,:]
#         bg_mean = np.mean(factor_bg)
#         bg_std = np.std(factor_bg)
#         if bg_mean != 0 and bg_std != 0:
#             result_table_p.loc[factor,:] = fg_table.loc[factor,:].apply(scipy.stats.norm.sf, args=(bg_mean, bg_std)) # sf is more accurate than cdf, return 1-cdf
#         else:
#             result_table_p.loc[factor,:] = 0
# #         result_table_fc.loc[factor,:] = fg_table.loc[factor,:].apply(cal_fc, **{'bg_mean': bg_mean})
#     print_log('chunk {i} finished calculation!'.format(i=i))
#     return result_table_p

# def cal_p_table_batch(fg_table, bg_table, n_cores=8):
#     print_log("Calculating enrichment, divide into {n} chunks...".format(n=n_cores))
#     fg_table_split = np.array_split(fg_table, n_cores)
#     args = [[table, bg_table, i] for (i, table) in enumerate(fg_table_split)]
#     with Pool(n_cores) as p:
#         result = p.starmap(cal_p_table, args)
#     print_log("Generating P value table ...")
#     result_table_p = pd.concat([i for i in result])
# #     print("INFO {time}, Generating FC value table ...".format(time=datetime.now()))
# #     result_table_fc = pd.concat([i[1] for i in result])
#     print_log('Finished calculation enrichment!')
#     return result_table_p

# def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):                
#     """                                                                                                   
#     consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
#     modified from https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python/7453313
#     """
#     pvalues = np.array(pvalues) 
#     n = pvalues.shape[0]     
#     new_pvalues = np.empty(n)
#     if correction_type == "Bonferroni":   
#         new_pvalues = n * pvalues
#     elif correction_type == "Bonferroni-Holm": 
#         values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]                                      
#         values.sort()
#         for rank, vals in enumerate(values):                                                              
#             pvalue, i = vals
#             new_pvalues[i] = (n-rank) * pvalue                                                            
#     elif correction_type == "Benjamini-Hochberg":                                                         
#         values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]                                      
#         values.sort()
#         values.reverse()
#         new_values = []
#         for i, vals in enumerate(values): 
#             rank = n - i
#             pvalue, index = vals 
#             new_values.append((n/rank) * pvalue) 
#         for i in range(0, int(n)-1): 
#             if new_values[i] < new_values[i+1]:                                                           
#                 new_values[i+1] = new_values[i]                                                           
#         for i, vals in enumerate(values):
#             pvalue, index = vals
#             new_pvalues[index] = new_values[i]                                                                                                                  
#     return new_pvalues


def cal_deviation(val, mean, std):
    return (val-mean)/std


def cal_deviation_table(fg_table, bg_table, i):
    print_log('chunk {i} calculating ...'.format(i=i))
    dts_cell_result_table_deviation = fg_table.copy()
    for factor in fg_table.index:
        factor_bg = bg_table.loc[factor, :].tolist()
        bg_mean = np.mean(factor_bg)
        bg_std = np.std(factor_bg)
        if bg_mean != 0 and bg_std != 0:
            dts_cell_result_table_deviation.loc[factor, :] = fg_table.loc[factor, :].apply(
                cal_deviation, **{'mean': bg_mean, 'std': bg_std})
        else:
            dts_cell_result_table_deviation.loc[factor, :] = 0
    print_log('chunk {i} finished calculation!'.format(i=i))
    return dts_cell_result_table_deviation


def cal_deviation_table_batch(fg_table, bg_table, n_cores=8):
    print_log("Calculating enrichment, divide into {n} chunks...".format(n=n_cores))
    fg_table_split = np.array_split(fg_table, n_cores)
    args = [[table, bg_table, i] for (i, table) in enumerate(fg_table_split)]
    with Pool(n_cores) as p:
        result = p.starmap(cal_deviation_table, args)
    dts_cell_result_table_deviation = pd.concat([i for i in result])
    print_log('Finished calculation enrichment!')
    return dts_cell_result_table_deviation


def cal_peak_norm_matrix(index_peak_number_path, peaks_number_path):
    index_peak_number = pd.read_csv(index_peak_number_path, sep='\t', header=None, index_col=0)
    data_peak_number = pd.read_csv(peaks_number_path, sep='\t', header=None, index_col=0)
    data_peak_number.index = data_peak_number.index.astype(str)
    peak_cell_index_norm_table = pd.DataFrame(np.zeros([index_peak_number.index.__len__(), data_peak_number.index.__len__()]),
                                              index=index_peak_number.index, columns=data_peak_number.index)
    for dts in peak_cell_index_norm_table.index:
        dts_number = index_peak_number.loc[dts, 1]
        peak_cell_index_norm_table.loc[dts, :] = [i/dts_number if i <= dts_number else dts_number/i for i in data_peak_number[1]]
    peak_cell_index_norm_table = peak_cell_index_norm_table.reindex(index=index_peak_number.index, columns=data_peak_number.index)

    peak_number_dts_mean_value = data_peak_number[1].mean()
    peak_number_dts_norm_series = data_peak_number[1].apply(lambda x: 2**-(abs(x-peak_number_dts_mean_value)/x))
    peak_number_dts_norm_series = peak_number_dts_norm_series.reindex(index=data_peak_number.index)

    peak_number_norm_coef = peak_cell_index_norm_table * peak_number_dts_norm_series

    return peak_number_norm_coef


def score_normalization(dataset_odds_ratio_df, dataset_fisher_df, peak_number_norm_coef):
    # fisher normalize matrix
    fisher_log_foreground = -np.log10(dataset_fisher_df)
    fisher_log_foreground_true_table = (fisher_log_foreground.T/fisher_log_foreground.max(1)).T
    fisher_log_foreground_true_table = fisher_log_foreground_true_table.reindex(index=dataset_odds_ratio_df.index, columns=dataset_odds_ratio_df.columns)

    peak_number_norm_coef_table = peak_number_norm_coef.reindex(
        index=dataset_odds_ratio_df.index, columns=dataset_odds_ratio_df.columns)

    # key multiply
    fg_dataset_cell_score_df = dataset_odds_ratio_df * fisher_log_foreground_true_table * peak_number_norm_coef_table
    return fg_dataset_cell_score_df
