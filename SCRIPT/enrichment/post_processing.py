#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   post_processing.py
@Time    :   2021/04/16 12:34:26
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''


# import subprocess
import scipy as sp
import numpy as np
import pandas as pd
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED
from SCRIPT.utilities.utils import print_log, excute_info


@excute_info('Summary result from dataset level to factor level.') 
def map_factor_on_ChIP(table):
    # map factor by id "_"
    factor_index_list = []
    for i in table.index:
        factor_name = i.split("_")
        factor_index_list.append(factor_name[1])
    table.loc[:,"Factor"] = factor_index_list
    return table.groupby("Factor").max()

# @excute_info('Extract cell by map dictionary ...', 'Finished all clusters!')
# def extract_by_cell_cluster(result_table, map_dict):
#     r_table = pd.DataFrame()
#     for key in result_table.columns:
#         key_cell_number = len(map_dict[key])
#         table_index = result_table.index
#         cell_bc = map_dict[key]
#         new_table_np = np.tile(result_table[key].to_numpy(),[key_cell_number,1]).T
#         tmp_table = pd.DataFrame(new_table_np, index = table_index, columns=cell_bc)
#         r_table = pd.concat([r_table, tmp_table], axis=1)
#     return r_table

# def p_to_z_transform(p_table):
#     z_table = -np.log10(p_table).T.apply(sp.stats.zscore, axis=0).T.copy()
#     return z_table

@excute_info('Generating merged anndata ... ', 'Finished Generating merged anndata!')
def merge_giggle_adata(adata, table, data_type, table2=''):
    # table : each row is a cell, each column is a TF/Gene
    new_adata = adata.copy()
    if data_type == 'integration':
        d_table = table.copy()
        d_table2 = table2.copy()
        new_adata.uns['ChIP'] = d_table
        new_adata.uns['motif'] = d_table2
        # extract factor
        chip_tfs = d_table.columns.tolist()
        motif_tfs = d_table2.columns.tolist()
        overlap_tf_list = list(set(chip_tfs).intersection(set(motif_tfs)))
        chip_other_tf_list = list(set(chip_tfs) - set(overlap_tf_list))
        motif_other_tf_list = list(set(motif_tfs) - set(overlap_tf_list))
        ovlp_chip_table = d_table[overlap_tf_list]
        ovlp_motif_table = d_table2[overlap_tf_list]
        unique_chip_table = d_table[chip_other_tf_list]
        unique_motif_table = d_table2[motif_other_tf_list]
        ovlp_final_table = pd.DataFrame()
        for ovlp_factor in overlap_tf_list:
            if np.std(ovlp_chip_table[ovlp_factor]) >= np.std(ovlp_motif_table[ovlp_factor]):
                ovlp_final_table[ovlp_factor] = ovlp_chip_table[ovlp_factor]
            else:
                ovlp_final_table[ovlp_factor] = ovlp_motif_table[ovlp_factor]
                
        final_table = pd.concat([ovlp_final_table,unique_chip_table,unique_motif_table], axis=1)
        final_table.columns = ['I_' + tf for tf in final_table.columns.tolist()]
        new_adata.uns['integrated'] = final_table
        new_adata.obs = pd.concat([new_adata.obs, final_table.reindex(new_adata.obs.index.tolist())], axis=1)
    else:
        d_table = table.copy()
        if data_type == 'ChIP-seq':
            new_adata.uns['ChIP'] = d_table
            d_table.columns = ['C_' + tf for tf in d_table.columns.tolist()]
        elif data_type == 'motif':
            new_adata.uns['motif'] = d_table
            d_table.columns = ['M_' + tf for tf in d_table.columns.tolist()]


        new_adata.obs = pd.concat([new_adata.obs, d_table.reindex(new_adata.obs.index.tolist())], axis=1)

    return new_adata
