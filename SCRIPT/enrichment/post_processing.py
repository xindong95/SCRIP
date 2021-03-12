import subprocess
import scipy as sp
import numpy as np
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED
from SCRIPT.utilities.utils import print_log, excute_info


@excute_info('Summary result from dataset level to factor level.') 
def map_factor_on_ChIP(table):
    # map factor by id "_"
    factor_index_list = []
    for i in table.index:
        factor_name = i.split("_")
        if len(factor_name) != 2:
            try:
                factor_index_list.append(factor_dict[int(factor_name[0])])
            except:
                print(factor_name)
                factor_index_list.append("None")
        else:
            factor_index_list.append(factor_name[1])
    table["Factor"] = factor_index_list
    return table.groupby("Factor").min()

@excute_info('Extract cell by map dictionary ...', 'Finished all clusters!')
def extract_by_cell_cluster(result_table, map_dict):
    r_table = pd.DataFrame()
    for key in result_table.columns:
        key_cell_number = len(map_dict[key])
        table_index = result_table.index
        cell_bc = map_dict[key]
        new_table_np = np.tile(result_table[key].to_numpy(),[key_cell_number,1]).T
        tmp_table = pd.DataFrame(new_table_np, index = table_index, columns=cell_bc)
        r_table = pd.concat([r_table, tmp_table], axis=1)
    return r_table

def p_to_z_transform(p_table):
    z_table = -np.log10(p_table).T.apply(sp.stats.zscore, axis=0).T.copy()
    return z_table

@excute_info('Generating merged anndata ... ', 'Finished Generating merged anndata!')
def merge_giggle_adata(adata, table, data_type, table2=''):
    # table : each row is a cell, each column is a TF/Gene
    new_adata = adata.copy()
    if data_type == 'integration':
        p_table = table.copy()
        p_table2 = table2.copy()
        new_adata.uns['ChIP_p'] = p_table
        new_adata.uns['motif_p'] = p_table2
        # table1, usually chip_table
        z_table = p_to_z_transform(p_table)
        new_adata.uns['ChIP_z'] = z_table
        # table2, usually motif table
        p_table2 = p_table2.reindex(index = p_table.index)
        z_table2 = p_to_z_transform(p_table2)
        new_adata.uns['motif_z'] = z_table2

        # extract factor
        chip_tfs = z_table.columns.tolist()
        motif_tfs = z_table2.columns.tolist()
        overlap_tf_list = list(set(chip_tfs).intersection(set(motif_tfs)))
        chip_other_tf_list = list(set(chip_tfs) - set(overlap_tf_list))
        motif_other_tf_list = list(set(motif_tfs) - set(overlap_tf_list))
        ovlp_chip_table = z_table[overlap_tf_list]
        ovlp_motif_table = z_table2[overlap_tf_list]
        unique_chip_table = z_table[chip_other_tf_list]
        unique_motif_table = z_table2[motif_other_tf_list]
        ovlp_final_table = pd.DataFrame()
        for ovlp_factor in overlap_tf_list:
            if np.std(ovlp_chip_table[ovlp_factor]) >= np.std(ovlp_motif_table[ovlp_factor]):
                ovlp_final_table[ovlp_factor] = ovlp_chip_table[ovlp_factor]
            else:
                ovlp_final_table[ovlp_factor] = ovlp_motif_table[ovlp_factor]
                
        final_table = pd.concat([ovlp_final_table,unique_chip_table,unique_motif_table], axis=1)
        final_table.columns = ['I_' + tf for tf in final_table.columns.tolist()]
        new_adata.uns['integrated_z'] = final_table
        new_adata.obs = pd.concat([new_adata.obs, final_table.reindex(new_adata.obs.index.tolist())], axis=1)
#         new_adata.uns['integrated_quantile'] = pd.DataFrame(quantile_transform(final_table.T, axis=0, output_distribution='uniform'), 
#                                                             index=final_table.index, columns=final_table.columns).T
    else:
        p_table = table.copy()
        if data_type == 'ChIP-seq':
            new_adata.uns['ChIP_p'] = p_table
            p_table.columns = ['C_' + tf for tf in p_table.columns.tolist()]
            z_table = p_to_z_transform(p_table)
            new_adata.uns['ChIP_z'] = z_table
        elif data_type == 'motif':
            new_adata.uns['motif_p'] = p_table
            p_table.columns = ['M_' + tf for tf in p_table.columns.tolist()]
            z_table = p_to_z_transform(p_table)
            new_adata.uns['ChIP_z'] = z_table

        new_adata.obs = pd.concat([new_adata.obs, z_table.reindex(new_adata.obs.index.tolist())], axis=1)
#         new_adata.uns['TF_quantile'] = pd.DataFrame(quantile_transform(z_table, axis=0, output_distribution='normal'), 
#                                                                index=z_table.index, columns=z_table.columns)
    return new_adata