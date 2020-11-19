def map_factor_on_ChIP(table):
    print_log('Summary result from dataset level to factor level.')
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

def extract_by_cell_cluster(result_table, map_dict):
    print_log('Extract cell by map dictionary ...')
    r_table = pd.DataFrame()
    for key in result_table.columns:
        key_cell_number = len(map_dict[key])
        table_index = result_table.index
        cell_bc = map_dict[key]
        new_table_np = np.tile(result_table[key].to_numpy(),[key_cell_number,1]).T
        tmp_table = pd.DataFrame(new_table_np, index = table_index, columns=cell_bc)
        r_table = pd.concat([r_table, tmp_table], axis=1)
    print_log('Finished all clusters!')
    return r_table

def merge_giggle_singlecell_experiment(adata, table, data_type, table2=''):
    print('Generating merged anndata ... ')
    # table : each row is a cell, each column is a TF/Gene
    new_adata = adata.copy()
    tmp_table = table.copy()
    if data_type == 'integration':
        tmp_table2 = table2.copy()
        new_adata.uns['table'] = tmp_table
        new_adata.uns['table2'] = tmp_table2
        # tmp table, usually chip_table
#         tmp_table = tmp_table.T.apply(correct_pvalues_for_multiple_testing, axis=0)
        tmp_table = -np.log10(tmp_table.T + 0.000001).apply(sp.stats.zscore, axis=0).T
        new_adata.uns['table_z'] = tmp_table
        # tmp table, usually motif table
        tmp_table2 = tmp_table2.reindex(index = tmp_table.index)
#         tmp_table2 = tmp_table2.T.apply(correct_pvalues_for_multiple_testing, axis=0)
        tmp_table2 = -np.log10(tmp_table2.T + 0.000001).apply(sp.stats.zscore, axis=0).T
        new_adata.uns['table2_z'] = tmp_table2
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
        new_adata.uns['integrated_TF'] = final_table
        new_adata.obs = pd.concat([new_adata.obs, final_table.reindex(new_adata.obs.index.tolist())], axis=1)
    else:
        if data_type == 'ChIP-seq':
            new_adata.uns['ChIP_p'] = tmp_table
            tmp_table.columns = ['C_' + tf for tf in tmp_table.columns.tolist()]
        elif data_type == 'motif':
            new_adata.uns['motif_p'] = tmp_table
            tmp_table.columns = ['M_' + tf for tf in tmp_table.columns.tolist()]
#         tmp_table = tmp_table.T.apply(correct_pvalues_for_multiple_testing, axis=0)
        tmp_table = -np.log10(tmp_table.T + 0.000001).apply(sp.stats.zscore, axis=0).T
        new_adata.obs = pd.concat([new_adata.obs, tmp_table.reindex(new_adata.obs.index.tolist())], axis=1)
    print_log('Finished Generating merged anndata!')
    return new_adata