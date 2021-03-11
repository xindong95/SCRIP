import os
import sys
import time
import pickle
import shutil
import scanpy as sc
import anndata as ad
from SCRIPT.enrichment.validation import check_para
from SCRIPT.enrichment.utils import EnrichRunInfo, time_estimate
from SCRIPT.utilities.utils import read_config, read_SingleCellExperiment_rds, print_log, excute_info
from SCRIPT.Constants import *


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
    
    # try:
    #     if os.listdir(bg_bed_path).__len__() != os.listdir(bg_result_path).__len__():
    #         shutil.rmtree(bg_result_path)
    #         print_log(log_list[0])
    #         search_giggle_batch(bg_bed_path, bg_result_path, index, n_cores)
    #         print_log(log_list[1])
    #     else:
    #         print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {bg_result_path}'.format(bg_result_path=bg_result_path))
    # except FileNotFoundError:
    #     print_log(log_list[0])
    #     search_giggle_batch(bg_bed_path, bg_result_path, index, n_cores)
    #     print_log(log_list[1])
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


def enrich(processed_adata, cell_feature_adata, project='',
           aggregate_peak_method='group', cell_number_per_group=50, cell_cutoff=20, peak_confidence=5, bg_iteration='auto', 
           chip_index='', motif_index='', reference_method='integration', result_adata_path='SCRIPT_computed.h5ad',
           yes=False, clean=True, n_cores=8, 
           processed_adata_path='NA', cell_feature_adata_path='NA', species='NA'):

    if project == '':
        tmp_chr_list = [chr(i) for i in range(ord("A"), ord("Z") + 1)] + [chr(i) for i in range(ord("a"), ord("z") + 1)] + [chr(i) for i in range(ord("0"), ord("9") + 1)]
        random.seed(time.time())
        tmp_prefix = str(time.time())[6:13].replace('.','') + '_' + ''.join(random.sample(tmp_chr_list, 4))
        project = 'SCRIPT_' + tmp_prefix

    if not os.path.exists(project):
        os.makedirs(project)
    project_abs_path = os.path.abspath(project)
    if bg_iteration != 'auto':
        bg_iteration = int(bg_iteration)
    ##################################
    ### check running information
    ##################################
    params = {
        'project':project_abs_path,
        'aggregate_peak_method':aggregate_peak_method,
        'cell_number_per_group':cell_number_per_group,
        'cell_cutoff':cell_cutoff,
        'peak_confidence':peak_confidence,
        'bg_iteration':bg_iteration,
        'chip_index':chip_index,
        'motif_index':motif_index,
        'reference_method':reference_method,
        'result_adata_path':result_adata_path,
        'clean':clean,
        'processed_adata_path':os.path.abspath(processed_adata_path),
        'cell_feature_adata_path':os.path.abspath(cell_feature_adata_path),
        'species':species
    }
    run_info = EnrichRunInfo(os.path.join(project,'run_info.json'), params)
    if not run_info.check_consist(params):
        print('Not consist! Remove the folder and continue? (Y/N)')
        while True:
            confirm_info = input()
            if confirm_info == 'y' or confirm_info == 'Y':
                shutil.rmtree(project)
                break
            elif confirm_info == 'N' or confirm_info == 'n':
                print_log('Good Bye!')
                sys.exit(0)
            else:
                print('Please type Y / N.')

    fg_bed_path = os.path.join(project, 'fg_files', 'fg_bed')
    fg_map_dict_path = os.path.join(project, 'fg_files', 'fg_bed.pk')
    fg_motif_result_path = os.path.join(project, 'fg_files', 'fg_motif_result')
    fg_chip_result_path = os.path.join(project, 'fg_files', 'fg_chip_result')
    bg_bed_path = os.path.join(project, 'fg_files', 'bg_bed')
    bg_map_dict_path = os.path.join(project, 'fg_files', 'bg_bed.pk')
    bg_chip_result_path = os.path.join(project, 'fg_files', 'bg_chip_result') 
    bg_motif_result_path = os.path.join(project, 'fg_files', 'bg_motif_result')

    ##################################
    ### pre-check
    ##################################
    print_log('Checking parameters ...')
    if bg_iteration == "auto":
        bg_iteration = int(processed_adata.shape[0] * 5 / cell_number_per_group) + 1
        check_para(processed_adata, cell_feature_adata, project,
                    aggregate_peak_method, cell_number_per_group, cell_cutoff, peak_confidence, bg_iteration, 
                    chip_index, motif_index, reference_method, result_adata_path,
                    yes, clean, n_cores)
    print_log('Estimating running time ...')
    if reference_method in ['integration', 'both', 'chip']:
        search_chip = True
    if reference_method in ['integration', 'both', 'motif']:
        search_motif = True
    elapse, future_time = time_estimate(cell_number = cell_feature_adata.shape[0], bg_iteration_number=bg_iteration, 
                                        peak_methods=aggregate_peak_method, cell_number_per_group=cell_number_per_group, 
                                        chip_process=search_chip, motif_process=search_motif, core=n_cores)
    print_log("It will take about {elapse} to process and finish at {future_time}.\n".format(elapse = elapse, future_time = future_time))
    if yes == False:
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
    if run_info['progress']['bg_bed_generation'] == 'No':
        if os.path.exists(bg_bed_path):
            shutil.rmtree(bg_bed_path)
        bg_map_dict = generate_background_bed(cell_feature_adata, bg_bed_path, bg_map_dict_path, 
                                              cell_number_per_group, bg_iteration, peak_confidence, n_cores)
        run_info.finish_stage('bg_bed_generation')
    else:  # run_info['progress']['bg_bed_generation'] == 'Finish'
        with open(bg_map_dict_path, "rb") as map_dict_file:
            bg_map_dict = pickle.load(map_dict_file)
    
    # if os.path.exists(bg_bed_path):
    #     if os.listdir(bg_bed_path).__len__() != bg_iteration:
    #         shutil.rmtree(bg_bed_path)
    #         print_log('Not empty folder, removed existing files!')
    #         bg_map_dict = generate_background_bed(cell_feature_adata, bg_bed_path, bg_map_dict_path, cell_number_per_group, bg_iteration, peak_confidence, n_cores)
    #     else:
    #         print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {bg_bed_path}'.format(bg_bed_path=bg_bed_path))
    #         with open(bg_map_dict_path, "rb") as map_dict_file:
    #             bg_map_dict = pickle.load(map_dict_file)
    # else:
    #     bg_map_dict = generate_background_bed(cell_feature_adata, bg_bed_path, bg_map_dict_path, cell_number_per_group, bg_iteration, peak_confidence, n_cores)

    if run_info['progress']['fg_bed_generation'] == 'No':
        if os.path.exists(fg_bed_path):
            shutil.rmtree(fg_bed_path)
        if aggregate_peak_method == "group":
            fg_map_dict = generate_cluster_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, cell_number_per_group, cell_cutoff, peak_confidence)
            run_info.finish_stage('fg_bed_generation')
        if aggregate_peak_method == "nearest":
            fg_map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, cell_number_per_group, peak_confidence, n_cores)
            run_info.finish_stage('fg_bed_generation')
    else:
        with open(fg_map_dict_path, "rb") as map_dict_file:
            fg_map_dict = pickle.load(map_dict_file)
    # if user generated foregroud, but diff depth in same folder, may report unaccurate result.
    # this is fast, we just remove it.
    # if aggregate_peak_method == "group":
    #     if os.path.exists(fg_bed_path):
    #         shutil.rmtree(fg_bed_path)
    #     fg_map_dict = generate_cluster_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, cell_number_per_group, cell_cutoff, peak_confidence)
    # # generate foreground peaks, if length same as cell number, we consider it has estimated, skip generation.
    # # if user generate same length background, but diff depth in same folder, may report unaccurate result.
    # if aggregate_peak_method == "nearest":
    #     if os.path.exists(fg_bed_path):
    #         if os.listdir(fg_bed_path).__len__() != processed_adata.shape[0]:
    #             shutil.rmtree(fg_bed_path)
    #             print_log('Not empty folder, removed existing files!')
    #             fg_map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, 
    #                                                 cell_number_per_group, peak_confidence, n_cores)
    #         else:
    #             print_log('WARNING: Using existing results, might wrong. If you want to rerun all results, press Ctrl-C to abort and delete {fg_bed_path}'.format(fg_bed_path=fg_bed_path))
    #             with open(fg_map_dict_path, "rb") as map_dict_file:
    #                 fg_map_dict = pickle.load(map_dict_file)
    #     else:
    #         fg_map_dict = generate_neighbor_bed(processed_adata, cell_feature_adata, fg_bed_path, fg_map_dict_path, 
    #                                             cell_number_per_group, peak_confidence, n_cores)
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
            regulation_adata = merge_giggle_singlecell_experiment(processed_adata, chip_result_p, 'integration', motif_result_p)
        else:
            regulation_adata = merge_giggle_singlecell_experiment(processed_adata, chip_result_p, 'ChIP-seq')
            regulation_adata = merge_giggle_singlecell_experiment(regulation_adata, motif_result_p, 'motif')
    elif search_chip == True and search_motif != True:
        regulation_adata = merge_giggle_singlecell_experiment(processed_adata, chip_result_p, 'ChIP-seq')
    else:
        regulation_adata = merge_giggle_singlecell_experiment(processed_adata, motif_result_p, 'motif')
    if result_store_path:
        regulation_adata.write(result_store_path)
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
    return regulation_adata


def run( args ):
    processed_adata_path = args.processed_experiment
    feature_matrix_path = args.feature_matrix
    species = args.species
    project = args.project
    result_adata_path = args.result_adata_path
    aggregate_peak_method = args.aggregate_peak_method
    cell_number_per_group = args.cell_number_per_group
    cell_cutoff = args.cell_cutoff
    peak_confidence = args.peak_confidence
    bg_iteration = args.bg_iter
    reference_method = args.reference
    yes = args.yes
    clean = args.clean
    n_cores = args.n_cores

    CONFIG, _ = read_config()
    if species == 'hs':
        chip_index = CONFIG['index']['human_chip_index']
        motif_index = CONFIG['index']['human_motif_index']
    elif species == 'mm':
        chip_index = CONFIG['index']['mouse_chip_index']
        motif_index = CONFIG['index']['mouse_motif_index']
    else:
        pass

    processed_adata = ad.read_h5ad(processed_adata_path)
    feature_matrix = sc.read_10x_h5(feature_matrix_path, gex_only=False)

    enrich(processed_adata, 
           feature_matrix, 
           project=project,
           aggregate_peak_method=aggregate_peak_method, 
           cell_number_per_group=cell_number_per_group, 
           cell_cutoff=cell_cutoff, 
           peak_confidence=peak_confidence, 
           bg_iteration=bg_iteration, 
           chip_index=chip_index, 
           motif_index=motif_index, 
           reference_method=reference_method, 
           result_adata_path=result_adata_path,
           yes=yes, 
           clean=clean, 
           n_cores=n_cores,
           processed_adata_path=processed_adata_path,
           cell_feature_adata_path=feature_matrix_path,
           species=species
           )

    return