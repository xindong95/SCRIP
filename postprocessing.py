def map_factor_on_ChIP(table, data_type):
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
    if data_type == "p":
        return table.groupby("Factor").min()
    elif data_type == "fc":
        return table.groupby("Factor").max()
    else:
        return None

def extract_by_cell_cluster(result_table, map_dict):
    r_table = pd.DataFrame()
    for key in result_table.columns:
        key_cell_number = len(map_dict[key])
        table_index = result_table.index
        cell_bc = map_dict[key]
        new_table_np = np.tile(result_table[key].to_numpy(),[key_cell_number,1]).T
        tmp_table = pd.DataFrame(new_table_np, index = table_index, columns=cell_bc)
        r_table = pd.concat([r_table, tmp_table], axis=1)
        sys.stdout.write('Finished Cluster ' + key + ' ... \t')
    print('Finished all!')
    return r_table

def merge_giggle_singlecell_experiment(adata, table, data_method):
    # giggle frame : each row is a cell, each column is a TF/Gene
    new_adata = adata.copy()
    tmp_table = table.copy()
    if data_method == 'ChIP-seq':
        new_adata.uns['ChIP_p'] = tmp_table
        tmp_table.columns = ['C_' + tf for tf in tmp_table.columns.tolist()]
    else:
        new_adata.uns['motif_p'] = tmp_table
        tmp_table.columns = ['M_' + tf for tf in tmp_table.columns.tolist()]
    tmp_table = -np.log10(tmp_table + 0.000001).apply(sp.stats.zscore)
    new_adata.obs = pd.concat([new_adata.obs, tmp_table.reindex(new_adata.obs.index.tolist())], axis=1)
    return new_adata