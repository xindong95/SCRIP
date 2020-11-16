def generate_beds(file_path, cells, input_mat, peak_confidence):
    peak_sub = input_mat.loc[:,cells].copy()
    peak_sub_matrix = peak_sub.loc[peak_sub.sum(1) > peak_confidence,]
    peaks = peak_sub_matrix.index.to_list()
    peaks = pd.DataFrame([p.split("_") for p in peaks])
    peaks.to_csv(file_path, sep="\t", header= None, index=None)

def generate_background_bed(input_mat, bg_bed_path, map_dict_store_path, step=50, iteration=1000, peak_confidence=5, n_cores=8):
    map_dict = {}
    if not os.path.exists(bg_bed_path):
        os.makedirs(bg_bed_path)
    cl_name = input_mat.columns.to_list()
    cnt = 0
    total_cnt = iteration
    executor = ThreadPoolExecutor(max_workers=n_cores)
    for i in range(0,iteration): 
        map_dict[i] = random.sample(cl_name, step)
        executor.submit(generate_beds, bg_bed_path + "/" + str(i) + ".bed", map_dict[i], input_mat, peak_confidence)
        if cnt%50 == 0:
            print_log("Finished {percentage:.2f} %".format(percentage = cnt*100/total_cnt), end="\r")
        cnt += 1    
    with open(map_dict_store_path, "wb") as map_dict_file:
        pickle.dump(map_dict, map_dict_file)
    return map_dict

def generate_cluster_bed(adata, input_mat, bed_path, map_dict_store_path, step=50, cell_cutoff=20, peak_confidence=5):
    metadata = adata.obs
    map_dict = {}
    if not os.path.exists(bed_path):
        os.makedirs(bed_path)
    cluster_info = metadata["seurat_clusters"].unique().to_list()
    for i in cluster_info:
        cluster_cell_name = metadata[metadata['seurat_clusters'] == i].index.to_list()
        start = range(0, len(cluster_cell_name), step)
        j = 0
        for s in start:
            cl_name = cluster_cell_name[s:(s + step)]
            key = str(i) + "_" + str(j)
            map_dict[key] = cl_name
            if len(cl_name) < cell_cutoff:
                print_log(key + " only have " + str(len(cl_name)), " cells, skip generate bed file.\n")
                continue
            else:
                generate_beds(bed_path + "/" + str(key) + ".bed", map_dict[key], input_mat, peak_confidence)
                j = j+1
    with open(map_dict_store_path, "wb") as map_dict_file:
        pickle.dump(map_dict, map_dict_file)
    return map_dict

def sub_coor_table_in_small_square(cell, coor_table, step, t):
    tmp = coor_table[(coor_table['X'] < coor_table.loc[cell, 'X'] + step * t) & 
                     (coor_table['X'] > coor_table.loc[cell, 'X'] - step * t) &  
                     (coor_table['Y'] < coor_table.loc[cell, 'Y'] + step * t) & 
                     (coor_table['Y'] > coor_table.loc[cell, 'Y'] - step * t)]
    return tmp.copy()

def find_nearest_cells(cell, coor_table, n_neighbor=20, step=None):
#   coor_table has two columns: X, Y, the index of coortable is cell barcodes
    if step == None:
        up_limit = coor_table.max()
        down_limit = coor_table.min()
        width = up_limit.X - down_limit.X
        height = up_limit.Y - down_limit.Y
        step = min(width, height)/200
    t = 1
    tmp = pd.DataFrame()
    while tmp.shape[0] < n_neighbor:
        tmp = sub_coor_table_in_small_square(cell, coor_table, step, t)
        t += 1
    squ_distance = [(tmp.loc[i,"X"] - tmp.loc[cell,"X"])**2 + (tmp.loc[i,"Y"] - tmp.loc[cell,"Y"])**2 for i in tmp.index]
    tmp['distance'] = squ_distance
    tmp = tmp.sort_values(by='distance', ascending=True)
    # since we use square to check the distance, check whether the point is out the circle
    while tmp.iloc[n_neighbor-1, 2] > ((step*t)**2):
        t += 1
        tmp = sub_coor_table_in_small_square(cell, coor_table, step, t)
        squ_distance = [(tmp.loc[i,"X"] - tmp.loc[cell,"X"])**2 + (tmp.loc[i,"Y"] - tmp.loc[cell,"Y"])**2 for i in tmp.index]
        tmp['distance'] = squ_distance
        tmp = tmp.sort_values(by='distance', ascending=True) 
    neighbor_bcs = tmp.index[0:n_neighbor].tolist()
    return neighbor_bcs

def generate_neighbor_bed(adata, input_mat, bed_path, map_dict_store_path, n_neighbor=50, peak_confidence=5, n_cores=8):
    coor_table = pd.DataFrame(adata.obsm['X_umap'], index = adata.obs.index, columns=["X","Y"])
    width = coor_table.max().X - coor_table.min().X
    height = coor_table.max().Y - coor_table.min().Y
    step = min(width, height)/200
    map_dict = {}
    if not os.path.exists(bed_path):
        os.makedirs(bed_path)
    cnt = 0
    total_cnt = adata.obs.index.__len__()
    executor = ThreadPoolExecutor(max_workers=n_cores)
    for cell in adata.obs.index:
        neighbor_cells = find_nearest_cells(cell, coor_table, n_neighbor)
        map_dict[cell] = neighbor_cells
        executor.submit(generate_beds, bed_path + "/" + str(cell) + ".bed", neighbor_cells, input_mat, peak_confidence)
        if cnt%50 == 0:
            print_log("Finished {percentage:.2f} %".format(percentage = cnt*100/total_cnt), end="\r")
        cnt += 1    
    with open(map_dict_store_path, "wb") as map_dict_file:
        pickle.dump(map_dict, map_dict_file)
    return map_dict
