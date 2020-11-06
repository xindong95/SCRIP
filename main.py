def time_estimate(cell_number, bg_iteration_number, peak_methods, cell_number_per_group, chip_process, motif_process, core):
#   generate bed: 0.5 second per cell,  search index: 5 seconds per cell, cal p: 2 seconds per factor
    if peak_methods == "group":
        bed_number = int(cell_number/cell_number_per_group) + bg_iteration_number
        seconds = 1.5 * bed_number
        if chip_process == True:
            seconds += 5 * bed_number
        if motif_process == True:
            seconds += 5 * bed_number
    elif peak_methods == "nearest":
        bed_number = cell_number + bg_iteration_number
        seconds = 1.5 * bed_number
        if chip_process == True:
            seconds += 5 * bed_number
        if motif_process == True:
            seconds += 5 * bed_number
    chip_factor_number = 5084
    motif_factor_number = 1147
    if chip_process == True:
        seconds += 2 * chip_factor_number
    if motif_process == True:
        seconds += 2 * motif_factor_number
    seconds = seconds/core
    now_time = datetime.now()
    future_time = (datetime.now() + timedelta(seconds=seconds)).strftime("%Y-%m-%d %H:%M:%S")
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    elapse = "{hour:0>2d}:{minute:0>2d}:{second:0>2d}".format(hour=int(h), minute=int(m), second=int(s))
    return elapse, future_time    

def check_para(processed_adata, feature_matrix, 
               aggregate_peak_method, cell_number_per_group, cell_cutoff, peak_confidence, 
               fg_bed_path, fg_map_dict_path, 
               bg_bed_path, bg_map_dict_path, bg_iteration, 
               search_chip, chip_index, chip_index_prefix, fg_chip_result_path, bg_chip_result_path,
               search_motif, motif_index, motif_index_prefix, fg_motif_result_path, bg_motif_result_path, 
               integration, confirm, clean, n_cores):
    if processed_adata.shape[0] > feature_matrix.shape[1]:
        print_log('WARNING: There are more cells in anndata than input feature matrix.')
    if aggregate_peak_method == 'group':
        if 'seurat_clusters' not in processed_adata.obs:
            print_log('There is no cluster infomation in anndata. "seurat_clusters in .obs is necessary."')
            sys.exit()
    if bg_iteration * cell_number_per_group < feature_matrix.shape[1] * 10:
        print_log('WARNING: Background iteration number less than 10 times number of cell. May be not enough to estimate background.')
    if re.match(r'chr.*_\d*_\d*', feature_matrix.index[0]) == None:
        print_log("feature_matrix's index should be like this: chr1_222222_333333. Each row is a feature and each column is a cell.")
        sys.exit()
    if cell_cutoff >= cell_number_per_group or peak_confidence >= cell_number_per_group:
        print_log("cell_cutoff or peak_confidence can not greater than cell_number_per_group.")
        sys.exit()
    if search_chip == True:
        if not os.path.exists(chip_index):
            print_log("chip_index do not exist!")
            sys.exit()
    if search_motif == True:
        if not os.path.exists(motif_index):
            print_log("motif_index do not exist!")
            sys.exit()
    if search_chip != True and search_motif != True:
        print_log("One of search_chip or search_motif must be chosen!")
        sys.exit()
    print('~~~~~\nWelcome to use SCRIPT. Here are your settings:')
    print('There is UMAP information in processed_adata. The aggregate peaks method is {aggregate_peak_method}.')
    if aggregate_peak_method == 'group':
        print('This will aggregate {cell_number_per_group} cells from same cluster')
        

def SCRIPT(processed_adata, feature_matrix, 
           aggregate_peak_method='group', cell_number_per_group=50, cell_cutoff=20, peak_confidence=5, 
           fg_bed_path='', fg_map_dict_path='', 
           bg_bed_path='', bg_map_dict_path='', bg_iteration='auto', 
           search_chip=True, chip_index='', chip_index_prefix='', fg_chip_result_path='', bg_chip_result_path='',
           search_motif=True, motif_index='', motif_index_prefix='', fg_motif_result_path='', bg_motif_result_path='', 
           integration=True, confirm=True, clean=True, n_cores=8):
    print_log('Checking paramaters ...')
    if bg_iteration == "auto":
        bg_iteration = int(feature_matrix.shape[1] * 10 / cell_number_per_group) + 1
    check_para(processed_adata, feature_matrix, 
               aggregate_peak_method, cell_number_per_group, cell_cutoff, peak_confidence, 
               fg_bed_path, fg_map_dict_path, 
               bg_bed_path, bg_map_dict_path, bg_iteration, 
               search_chip, chip_index, chip_index_prefix, fg_chip_result_path, bg_chip_result_path,
               search_motif, motif_index, motif_index_prefix, fg_motif_result_path, bg_motif_result_path, 
               integration, confirm, clean, n_cores)
    print_log('Estimating running time ...')
    elapse, future_time = time_estimate(cell_number = feature_matrix.shape[1], bg_iteration_number=bg_iteration, 
                                        peak_methods=aggregate_peak_method, cell_number_per_group=cell_number_per_group, 
                                        chip_process=search_chip, motif_process=search_motif, core=n_cores)
    print_log("It will take about {elapse} to process and finish at {future_time}.".format(elapse = elapse, future_time = future_time))
    if confirm == True:
        print('Type "Y" to continue processing, "N" to abort.')
        while True:
            confirm_info = input()
            if confirm_info == 'y' or confirm_info == 'Y':
                break
            elif confirm_info == 'N' or confirm_info == 'n':
                print_log('Good Bye!')
                sys.exit(0)
            else:
                print('Please type Y / N.')
    
    # generate background peak, if length same as iteration, we consider it has estimated, skip generation.
    # if user generate same length background, but diff depth, may report unaccurate result.
    try:
        tmp_files_list = os.listdir(bg_bed_path)
        if tmp_files_list.__len__() != bg_iteration:
            shutil.rmtree(bg_bed_path)
            print_log('Start generating background beds ...')
            bg_map_dict = generate_background_bed(feature_matrix, bg_bed_path, bg_map_dict_path, cell_number_per_group, bg_iteration, peak_confidence, n_cores)
            print_log('Finished generating background beds!')
        else:
            print_log('WARNING: Using existing results, might wrong.')
            with open(bg_map_dict_path, "rb") as map_dict_file:
                bg_map_dict = pickle.load(map_dict_file)
    except FileNotFoundError:
        print_log('Start generating background beds ...')
        bg_map_dict = generate_background_bed(feature_matrix, bg_bed_path, bg_map_dict_path, cell_number_per_group, bg_iteration, peak_confidence, n_cores)
        print_log('Finished generating background beds!')

    # if user generated foregroud, but diff depth in same folder, may report unaccurate result.
    # this is fast, we just remove it.
    if aggregate_peak_method == "group":
        try:
            shutil.rmtree(fg_bed_path)
        except:
            pass
        print_log('Start generating group beds ...')
        fg_map_dict = generate_cluster_bed(processed_adata, feature_matrix, fg_bed_path, fg_map_dict_path, cell_number_per_group, cell_cutoff, peak_confidence)
        print_log('Finished generating group beds!')
    if aggregate_peak_method == "nearest":
        # generate foreground peaks, if length same as cell number, we consider it has estimated, skip generation.
        # if user generate same length background, but diff depth in same folder, may report unaccurate result.
        try:
            tmp_files_list = os.listdir(fg_bed_path)
            if tmp_files_list.__len__() != feature_matrix.shape[1]:
                shutil.rmtree(fg_bed_path)
                print_log('Start generating nearest neighbor cells beds ...')
                fg_map_dict = generate_neighbor_bed(processed_adata, feature_matrix, fg_bed_path, fg_map_dict_path, 
                                                    cell_number_per_group, peak_confidence, n_cores)
                print_log('Finished generating nearest neighbor cells beds!')
            else:
                print_log('WARNING: Using existing results, might wrong.')
                with open(fg_map_dict_path, "rb") as map_dict_file:
                    fg_map_dict = pickle.load(map_dict_file)
        except FileNotFoundError:
            print_log('Start generating nearest neighbor cells beds ...')
            fg_map_dict = generate_neighbor_bed(processed_adata, feature_matrix, fg_bed_path, fg_map_dict_path, 
                                                cell_number_per_group, peak_confidence, n_cores)
            print_log('Finished generating nearest neighbor cells beds!')
            
    if search_chip == True:
        print_log('Start searching background beds from ChIP-seq index ...')
        search_giggle_batch(bg_bed_path, bg_chip_result_path, chip_index, n_cores)
        print_log('Finished searching background beds from ChIP-seq index!')
        
        print_log('Start searching foreground beds from ChIP-seq index ...')
        search_giggle_batch(fg_bed_path, fg_chip_result_path, chip_index, n_cores)
        print_log('Finished searching factors from ChIP-seq index!')
        
        chip_bg_result = read_giggle_result(bg_chip_result_path, filename_split=".", index_prefix=chip_index_prefix)
        print_log('Finished reading background ChIP-seq index search result!')
        
        chip_result = read_giggle_result(fg_chip_result_path, filename_split=".", index_prefix=chip_index_prefix)
        print_log('Finished reading foreground ChIP-seq index search result!')
        
        chip_result_p = cal_p_table_batch(chip_result, chip_bg_result, n_cores)
        if aggregate_peak_method == "group":
            chip_result_p = extract_by_cell_cluster(chip_result_p.copy(), fg_map_dict)
        chip_result_p = map_factor_on_ChIP(chip_result_p).T
    if search_motif == True:
        print_log('Start searching background beds from motif index ...')
        search_giggle_batch(bg_bed_path, bg_motif_result_path, motif_index, n_cores)
        print_log('Finished searching background beds from motif index!')
        
        print_log('Start searching foreground beds from motif index ...')
        search_giggle_batch(fg_bed_path, fg_motif_result_path, motif_index, n_cores)
        print_log('Finished searching factors from motif index!')
        
        motif_bg_result = read_giggle_result(bg_motif_result_path, filename_split=".", index_prefix=motif_index_prefix)
        print_log('Finished reading background motif index search result!')
        
        motif_result = read_giggle_result(fg_motif_result_path, filename_split=".", index_prefix=motif_index_prefix)
        print_log('Finished reading foreground motif index search result!')
        
        motif_result_p = cal_p_table_batch(motif_result, motif_bg_result, n_cores)
        if aggregate_peak_method == "group":
            motif_result_p = extract_by_cell_cluster(motif_result_p.copy(), fg_map_dict)
        motif_result_p = motif_result_p.T
    if search_chip == True and search_motif == True:
        if integration == True:
            regualtion_adata = merge_giggle_singlecell_experiment(processed_adata, chip_result_p, 'integration', motif_result_p)
        else:
            regualtion_adata = merge_giggle_singlecell_experiment(processed_adata, chip_result_p, 'ChIP-seq')
            regualtion_adata = merge_giggle_singlecell_experiment(regualtion_adata, motif_result_p, 'motif')
    elif search_chip == True and search_motif != True:
        regualtion_adata = merge_giggle_singlecell_experiment(processed_adata, chip_result_p, 'ChIP-seq')
    else:
        regualtion_adata = merge_giggle_singlecell_experiment(processed_adata, motif_result_p, 'motif')
    if clean == True:
        shutil.rmtree(fg_bed_path)
        shutil.rmtree(bg_bed_path)
        if search_chip == True:
            shutil.rmtree(fg_chip_result_path)
            shutil.rmtree(bg_chip_result_path)
        if search_motif == True:
            shutil.rmtree(fg_motif_result_path)
            shutil.rmtree(bg_motif_result_path)
    return regualtion_adata        
        