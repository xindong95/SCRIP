def time_estimate(cell_number, bg_iteration_number, peak_methods, cell_number_per_group, chip_process, motif_process, core):
#   generate bed: 0.5 second per cell,  search index: 10 seconds per cell, cal p: 2 seconds per factor
    if peak_methods == "group":
        bed_number = int(cell_number/cell_number_per_group) + bg_iteration_number
        seconds = 3 * bed_number
        if chip_process == True:
            seconds += 10 * bed_number
        if motif_process == True:
            seconds += 10 * bed_number
    elif peak_methods == "nearest":
        seconds = int(cell_number * core / 50) # single core process, find nearest cells function
        bed_number = cell_number + bg_iteration_number
        seconds += 3 * bed_number
        if chip_process == True:
            seconds += 10 * bed_number
        if motif_process == True:
            seconds += 10 * bed_number
    chip_factor_number = 5069
    motif_factor_number = 1143
    if chip_process == True:
        seconds += 2 * chip_factor_number
    if motif_process == True:
        seconds += 2 * motif_factor_number
    seconds = seconds / core
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
               search_chip, chip_index, fg_chip_result_path, bg_chip_result_path,
               search_motif, motif_index, fg_motif_result_path, bg_motif_result_path, 
               integration, store, result_store_path, confirm, clean, n_cores):
    if processed_adata.shape[0] > feature_matrix.shape[0]:
        print_log('WARNING: There are more cells in anndata than input feature matrix.')
    if aggregate_peak_method == 'group':
        if 'seurat_clusters' not in processed_adata.obs:
            print_log('There is no cluster infomation in anndata. "seurat_clusters in .obs is necessary."')
            sys.exit()
    if bg_iteration * cell_number_per_group < feature_matrix.shape[0] * 5:
        print_log('WARNING: Background iteration number less than 5 times number of cell. May be not enough to estimate background.')
    if re.match(r'chr.*_\d*_\d*', feature_matrix.var_names[0]) == None:
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
    if store == True and (type(result_store_path) != str or not result_store_path.endswith('.h5ad')):
        print_log("The result_store_path must be set as .h5ad file when store is True!")
        sys.exit()
    # prepare information
    information = '\n~~~~~\nWelcome to use SCRIPT. Here are your settings:\n\n'
    information += 'There are UMAP information and cluster information in your processed_adata. '
    information += 'The aggregate peaks method is {aggregate_peak_method}.\n'.format(aggregate_peak_method=aggregate_peak_method)
    if aggregate_peak_method == 'group':
        information += 'This will aggregate {cell_number_per_group} cells from the same cluster, '.format(cell_number_per_group=cell_number_per_group) 
        information += 'the group that cell number is less than {cell_cutoff} will not generate.\n'.format(cell_cutoff=cell_cutoff)
        information += 'For each group, peaks that confidence is less than {peak_confidence} will be deleted.\n'.format(peak_confidence=peak_confidence)
        information += 'This setting is suitable for the background as well. '
    else:
        information += 'This will borrow peaks from {cell_number_per_group} nearest cells for every cell. '.format(cell_number_per_group=cell_number_per_group)
        information += 'For each cell, peaks that confidence is less than {peak_confidence} will be deleted.\n'.format(peak_confidence=peak_confidence)
    information += 'SCRIPT randomly aggregates {cell_number_per_group} cells peaks '.format(cell_number_per_group=cell_number_per_group)
    information += 'and deletes the peaks that confidence less than {peak_confidence}. '.format(peak_confidence=peak_confidence)
    information += 'This will iterate {bg_iteration} times.\n'.format(bg_iteration=bg_iteration)
    information += 'The background peaks will be stored in {bg_bed_path}, '.format(bg_bed_path=bg_bed_path)
    information += 'foreground peaks will be stored in {fg_bed_path}.\n'.format(fg_bed_path=fg_bed_path)
    if search_chip == True and search_motif == True:
        information += 'SCRIPT will search ChIP-seq and motif index both.\n'
        if integration == True:
            information += 'After that, SCRIPT will integrate into one result based on variability from each result.\n'
    elif search_chip == True:
        information += 'SCRIPT will search ChIP-seq index.\n'
    elif search_motif == True:
        information += 'SCRIPT will search motif index.\n'
    if search_chip == True:
        information += 'ChIP-seq index locates at {chip_index}.\n'.format(chip_index=chip_index)
        information += 'ChIP-seq background result will be stored in {bg_chip_result_path}, '.format(bg_chip_result_path=bg_chip_result_path)
        information += 'foreground result will be stored in {fg_chip_result_path}.\n'.format(fg_chip_result_path=fg_chip_result_path)
    if search_motif == True:
        information += 'Motif index locates at {motif_index}.\n'.format(motif_index=motif_index)
        information += 'Motif background result will be stored in {bg_motif_result_path}, '.format(bg_motif_result_path=bg_motif_result_path)
        information += 'foreground result will be stored in {fg_motif_result_path}.\n'.format(fg_motif_result_path=fg_motif_result_path)
    if store == True:
        information += 'Computed result will be stored in {result_store_path}.\n'.format(result_store_path=result_store_path)
    else:
        information += 'Computed result will not be stored.\n'.format(result_store_path=result_store_path)
    information += 'All folders will{not_string} be removed after processing. '.format(not_string='' if clean == True else ' not')
    information += 'All processes will use {n_cores} cores.\n~~~~~\n'.format(n_cores=n_cores)
    print(information)
        
def process(tp, bg_bed_path, bg_result_path, fg_bed_path, fg_result_path, index, n_cores, fg_map_dict, aggregate_peak_method):
    if tp == 'chip':
        log_list = [
            'Start searching background beds from ChIP-seq index ...',
            'Finished searching background beds from ChIP-seq index!',
            'Start searching foreground beds from ChIP-seq index ...',
            'Finished searching factors from ChIP-seq index!',
            'Finished reading background ChIP-seq index search result!',
            'Finished reading foreground ChIP-seq index search result!'
        ]
    else:
        log_list = [
            'Start searching background beds from motif index ...',
            'Finished searching background beds from motif index!',
            'Start searching foreground beds from motif index ...',
            'Finished searching factors from motif index!',
            'Finished reading background motif index search result!',
            'Finished reading foreground motif index search result!'
        ]
    try:
        if os.listdir(bg_bed_path).__len__() != os.listdir(bg_result_path).__len__():
            shutil.rmtree(bg_result_path)
            print_log(log_list[0])
            search_giggle_batch(bg_bed_path, bg_result_path, index, n_cores)
            print_log(log_list[1])
        else:
            print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {bg_result_path}'.format(bg_result_path=bg_result_path))
    except FileNotFoundError:
        print_log(log_list[0])
        search_giggle_batch(bg_bed_path, bg_result_path, index, n_cores)
        print_log(log_list[1])
    try:
        if os.listdir(fg_bed_path).__len__() != os.listdir(fg_result_path).__len__():
            shutil.rmtree(fg_result_path)
            print_log(log_list[2])
            search_giggle_batch(fg_bed_path, fg_result_path, index, n_cores)
            print_log(log_list[3])
        else:
            print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {fg_result_path}'.format(fg_result_path=fg_result_path))
    except FileNotFoundError:
        print_log(log_list[2])
        search_giggle_batch(fg_bed_path, fg_result_path, index, n_cores)
        print_log(log_list[3])
    bg_result = read_giggle_result_batch(bg_result_path, n_cores)
    print_log(log_list[4])
    result = read_giggle_result_batch(fg_result_path, n_cores)
    print_log(log_list[5])
    result_p = cal_p_table_batch(result, bg_result, n_cores)
    if aggregate_peak_method == "group":
        result_p = extract_by_cell_cluster(result_p.copy(), fg_map_dict)
    if tp == 'chip':
        result_p = map_factor_on_ChIP(result_p).T
    else:
        result_p = result_p.T
    return result_p


def SCRIPT(processed_adata, feature_matrix, project='',
           aggregate_peak_method='group', cell_number_per_group=50, cell_cutoff=20, peak_confidence=5, bg_iteration='auto', 
           search_chip=True, chip_index='', search_motif=True, motif_index='',            
           integration=True, store=False, result_store_path='SCRIPT_computed.h5ad',
           confirm=True, clean=True, n_cores=8):
    if project == '':
        tmp_chr_list = [chr(i) for i in range(ord("A"), ord("Z") + 1)] + [chr(i) for i in range(ord("a"), ord("z") + 1)] + [chr(i) for i in range(ord("0"), ord("9") + 1)]
        random.seed(time.time())
        tmp_prefix = str(time.time())[6:13].replace('.','') + '_' + ''.join(random.sample(tmp_chr_list, 4))
        fg_bed_path='.SCRIPT_temp/{tmp_prefix}/fg_bed'.format(tmp_prefix = tmp_prefix)
        fg_map_dict_path='.SCRIPT_temp/{tmp_prefix}/fg_bed.pk'.format(tmp_prefix = tmp_prefix)
        bg_bed_path='.SCRIPT_temp/{tmp_prefix}/bg_bed'.format(tmp_prefix = tmp_prefix)
        bg_map_dict_path='.SCRIPT_temp/{tmp_prefix}/bg_bed.pk'.format(tmp_prefix = tmp_prefix)
        fg_chip_result_path='.SCRIPT_temp/{tmp_prefix}/fg_chip_result'.format(tmp_prefix = tmp_prefix)
        bg_chip_result_path='.SCRIPT_temp/{tmp_prefix}/bg_chip_result'.format(tmp_prefix = tmp_prefix)
        fg_motif_result_path='.SCRIPT_temp/{tmp_prefix}/fg_motif_result'.format(tmp_prefix = tmp_prefix)
        bg_motif_result_path='.SCRIPT_temp/{tmp_prefix}/bg_motif_result'.format(tmp_prefix = tmp_prefix)
    else:
        fg_bed_path='{project}/fg_bed'.format(project = project)
        fg_map_dict_path='{project}/fg_bed.pk'.format(project = project)
        bg_bed_path='{project}/bg_bed'.format(project = project)
        bg_map_dict_path='{project}/bg_bed.pk'.format(project = project)
        fg_chip_result_path='{project}/fg_chip_result'.format(project = project)
        bg_chip_result_path='{project}/bg_chip_result'.format(project = project)
        fg_motif_result_path='{project}/fg_motif_result'.format(project = project)
        bg_motif_result_path='{project}/bg_motif_result'.format(project = project)
    ##################################
    ### pre-check
    ##################################
    print_log('Checking paramaters ...')
    if bg_iteration == "auto":
        bg_iteration = int(processed_adata.shape[0] * 5 / cell_number_per_group) + 1
    check_para(processed_adata, feature_matrix, 
               aggregate_peak_method, cell_number_per_group, cell_cutoff, peak_confidence, 
               fg_bed_path, fg_map_dict_path, 
               bg_bed_path, bg_map_dict_path, bg_iteration, 
               search_chip, chip_index, fg_chip_result_path, bg_chip_result_path,
               search_motif, motif_index, fg_motif_result_path, bg_motif_result_path, 
               integration, store, result_store_path, confirm, clean, n_cores)
    print_log('Estimating running time ...')
    elapse, future_time = time_estimate(cell_number = feature_matrix.shape[0], bg_iteration_number=bg_iteration, 
                                        peak_methods=aggregate_peak_method, cell_number_per_group=cell_number_per_group, 
                                        chip_process=search_chip, motif_process=search_motif, core=n_cores)
    print_log("It will take about {elapse} to process and finish at {future_time}.\n".format(elapse = elapse, future_time = future_time))
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
    ##################################
    ### bed generation
    ##################################
    # generate background peak, if length same as iteration, we consider it has estimated, skip generation.
    # if user generate same length background, but diff depth, may report unaccurate result.
    if os.path.exists(bg_bed_path):
        if os.listdir(bg_bed_path).__len__() != bg_iteration:
            shutil.rmtree(bg_bed_path)
            print_log('Not empty folder, removed existing files!')
            bg_map_dict = generate_background_bed(feature_matrix, bg_bed_path, bg_map_dict_path, cell_number_per_group, bg_iteration, peak_confidence, n_cores)
        else:
            print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {bg_bed_path}'.format(bg_bed_path=bg_bed_path))
            with open(bg_map_dict_path, "rb") as map_dict_file:
                bg_map_dict = pickle.load(map_dict_file)
    else:
        bg_map_dict = generate_background_bed(feature_matrix, bg_bed_path, bg_map_dict_path, cell_number_per_group, bg_iteration, peak_confidence, n_cores)
    # if user generated foregroud, but diff depth in same folder, may report unaccurate result.
    # this is fast, we just remove it.
    if aggregate_peak_method == "group":
        if os.path.exists(fg_bed_path):
            shutil.rmtree(fg_bed_path)
        fg_map_dict = generate_cluster_bed(processed_adata, feature_matrix, fg_bed_path, fg_map_dict_path, cell_number_per_group, cell_cutoff, peak_confidence)
    # generate foreground peaks, if length same as cell number, we consider it has estimated, skip generation.
    # if user generate same length background, but diff depth in same folder, may report unaccurate result.
    if aggregate_peak_method == "nearest":
        if os.path.exists(fg_bed_path):
            if os.listdir(fg_bed_path).__len__() != processed_adata.shape[0]:
                shutil.rmtree(fg_bed_path)
                print_log('Not empty folder, removed existing files!')
                fg_map_dict = generate_neighbor_bed(processed_adata, feature_matrix, fg_bed_path, fg_map_dict_path, 
                                                    cell_number_per_group, peak_confidence, n_cores)
            else:
                print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {fg_bed_path}'.format(fg_bed_path=fg_bed_path))
                with open(fg_map_dict_path, "rb") as map_dict_file:
                    fg_map_dict = pickle.load(map_dict_file)
        else:
            fg_map_dict = generate_neighbor_bed(processed_adata, feature_matrix, fg_bed_path, fg_map_dict_path, 
                                                cell_number_per_group, peak_confidence, n_cores)
    ##################################
    ### Search giggle and compute enrich score
    ##################################
    if search_chip == True:
        chip_result_p = process('chip', bg_bed_path, bg_chip_result_path, fg_bed_path, fg_chip_result_path, 
                                chip_index, n_cores, fg_map_dict, aggregate_peak_method)
    if search_motif == True:
        motif_result_p = process('motif', bg_bed_path, bg_motif_result_path, fg_bed_path, fg_motif_result_path, 
                                motif_index, n_cores, fg_map_dict, aggregate_peak_method)
    ##################################
    ### Summary results
    ##################################
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
    if result_store_path:
        regualtion_adata.write(result_store_path)
    ##################################
    ### Clean files
    ##################################
    if clean == True:
        try:
            shutil.rmtree(fg_bed_path)
            shutil.rmtree(bg_bed_path)
            if search_chip == True:
                shutil.rmtree(fg_chip_result_path)
                shutil.rmtree(bg_chip_result_path)
            if search_motif == True:
                shutil.rmtree(fg_motif_result_path)
                shutil.rmtree(bg_motif_result_path)
        except:
            pass
    return regualtion_adata
