import datetime
from SCRIPT.utils.utils import print_log

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
               integration, result_store_path, confirm, clean, n_cores):
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