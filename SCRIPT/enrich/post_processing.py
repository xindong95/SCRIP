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
    z_table = -np.log10(np.sqrt(p_table+1e-8)).T.apply(sp.stats.zscore, axis=0).T.copy()
    return z_table

@excute_info('Generating merged anndata ... ', 'Finished Generating merged anndata!')
def merge_giggle_singlecell_experiment(adata, table, data_type, table2=''):
    # table : each row is a cell, each column is a TF/Gene
    new_adata = adata.copy()
    if data_type == 'integration':
        p_table = table.copy()
        p_table2 = table2.copy()
        new_adata.uns['ChIP_p'] = p_table
        new_adata.uns['motif_p'] = p_table2
        # tmp table, usually chip_table
#         tmp_table = tmp_table.T.apply(correct_pvalues_for_multiple_testing, axis=0)
#         tmp_table = -np.log10(tmp_table + 1e-10).apply(sp.stats.zscore, axis=0)
        new_adata.uns['ChIP_z'] = p_to_z_transform(p_table)
        # tmp table, usually motif table
        p_table2 = p_table2.reindex(index = p_table.index)
#         tmp_table2 = tmp_table2.T.apply(correct_pvalues_for_multiple_testing, axis=0)
#         tmp_table2 = -np.log10(tmp_table2 + 1e-10).apply(sp.stats.zscore, axis=0)
        new_adata.uns['motif_z'] = p_to_z_transform(p_table2)
        # extract factor
        chip_tfs = tmp_table.columns.tolist()
        motif_tfs = tmp_table2.columns.tolist()
        overlap_tf_list = list(set(chip_tfs).intersection(set(motif_tfs)))
        chip_other_tf_list = list(set(chip_tfs) - set(overlap_tf_list))
        motif_other_tf_list = list(set(motif_tfs) - set(overlap_tf_list))
        ovlp_chip_table = tmp_table[overlap_tf_list]
        ovlp_motif_table = tmp_table2[overlap_tf_list]
        unique_chip_table = tmp_table[chip_other_tf_list]
        unique_motif_table = tmp_table2[motif_other_tf_list]
        ovlp_final_table = pd.DataFrame()
        for ovlp_factor in overlap_tf_list:
            if np.std(ovlp_chip_table[ovlp_factor]) >= np.std(ovlp_motif_table[ovlp_factor]):
                ovlp_final_table[ovlp_factor] = ovlp_chip_table[ovlp_factor]
            else:
                ovlp_final_table[ovlp_factor] = ovlp_motif_table[ovlp_factor]
        final_table = pd.concat([ovlp_final_table,unique_chip_table,unique_motif_table], axis=1)
        final_table.columns = ['I_' + tf for tf in final_table.columns.tolist()]
        new_adata.uns['integrated_TF_z'] = final_table
        new_adata.obs = pd.concat([new_adata.obs, final_table.reindex(new_adata.obs.index.tolist())], axis=1)
        new_adata.uns['integrated_TF_quantile'] = pd.DataFrame(quantile_transform(final_table, axis=0, output_distribution='normal'), index=final_table.index, columns=final_table.columns)
    else:
        tmp_table = table.copy()
        if data_type == 'ChIP-seq':
            new_adata.uns['ChIP_p'] = tmp_table
            tmp_table.columns = ['C_' + tf for tf in tmp_table.columns.tolist()]
            tmp_table = -np.log10(np.sqrt(tmp_table+1e-8)).T.apply(sp.stats.zscore, axis=0).T.copy()
            new_adata.uns['ChIP_z'] = tmp_table
        elif data_type == 'motif':
            new_adata.uns['motif_p'] = tmp_table
            tmp_table.columns = ['M_' + tf for tf in tmp_table.columns.tolist()]
            tmp_table = -np.log10(np.sqrt(tmp_table+1e-8)).T.apply(sp.stats.zscore, axis=0).T.copy()
            new_adata.uns['ChIP_z'] = tmp_table
        
        new_adata.obs = pd.concat([new_adata.obs, tmp_table.reindex(new_adata.obs.index.tolist())], axis=1)
        new_adata.uns['integrated_TF_quantile'] = pd.DataFrame(quantile_transform(final_table, axis=0, output_distribution='normal'), index=final_table.index, columns=final_table.columns)
    return new_adata


