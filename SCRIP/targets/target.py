#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   target.py
@Time    :   2021/05/27 12:46:26
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

# from multiprocessing import Process, Pool
import os
import sys
import scipy
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from SCRIP.utilities.utils import read_config
# from SCRIP.enrichment.bed_generation import find_nearest_cells


def propCalc(genes_list, peaks_list):
    PEAK = 0
    GENE = 1

    w = genes_list + peaks_list
    D = {}
    A = {}
    w.sort()
    for elem in w:
        if elem[2] == PEAK:
            A[elem[-1]] = [elem[0], elem[1]]
            D[elem[-1]] = float("inf")
        else:
            for peak in list(A.keys()):
                p_elem = A[peak]
                if p_elem[0] == elem[0]:
                    D[peak] = elem[1] - p_elem[1]
                    A.pop(peak)
                else:
                    A.pop(peak)

    w.reverse()
    for elem in w:
        if elem[2] == PEAK:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            for peak in list(A.keys()):
                p_elem = A[peak]
                if p_elem[0] == elem[0]:
                    D[peak] = min(p_elem[1] - elem[1], D[peak])
                    A.pop(peak)
                else:
                    A.pop(peak)
    D_values = D.values()
    D_len = float(len(D_values))
    peak_within_1k = [i for i in D.values() if i <= 1000]
    prop = len(peak_within_1k) / D_len

    if prop >= 0.2:
        print("""%.2f%% of peaks are located in the promotor region. 
This is more than the 20%% promoter-type threshold so the the decay distance is set to 1.0kb, 
appropriate for promotor-type analysis. The half decay distance can be specified in the parameters.""" % (100*prop))
    else:
        print("""%.2f%% of peaks are located in the promotor region. 
This is less than the 20%% promoter-type threshold so the the decay distance is set to 10.0kb, 
appropriate for enhancer-type analysis. The half decay distance can be specified in the parameters.""" % (prop*100))
    return prop


def geneInfoSimple(gene_bed):
    genes_info = []
    genes_list = []
    fhd = open(gene_bed, 'rt')
    fhd.readline()  # skip the first line. In our current gene txt file, there is no '#' in the first line. We need to, perhaps, use the 'ExtractGeneInfo' function.
    for line in fhd:
        line = line.strip().split('\t')
        if not line[0].startswith('#'):
            if line[3] == "+":
                genes_info.append((line[2].replace('chr', ''), int(line[4]), 1, "%s@%s@%s" % (line[12], line[2], line[4])))
            else:
                genes_info.append((line[2].replace('chr', ''), int(line[5]), 1, "%s@%s@%s" % (line[12], line[2], line[5])))
                # gene_info [chrom, tss, 1, gene_unique]
    fhd.close()
    genes_info = list(set(genes_info))
    for igene in range(len(genes_info)):
        tmp_gene = list(genes_info[igene])
        genes_list.append(tmp_gene[3])
        tmp_gene[3] = igene
        genes_info[igene] = tmp_gene
    return genes_info, genes_list


def ExtractGeneInfo(gene_bed):
    """Extract gene information from gene bed file."""

    bed = pd.read_csv(gene_bed, sep="\t", header=0, index_col=False)
    bed['transcript'] = [x.strip().split(".")[0] for x in bed['name'].tolist()]
    bed['tss'] = bed.apply(lambda x: x['txStart'] if x['strand'] == '+' else x['txEnd'], axis=1)
    ### adjacent P+GB
    bed["start"] = bed.apply(lambda x: x['txStart']-2000 if x['strand'] == '+' else x['txStart'], axis=1)
    bed["end"] = bed.apply(lambda x: x['txEnd']+2000 if x['strand'] == '-' else x['txEnd'], axis=1)
    bed['promoter'] = bed.apply(lambda x: tuple([x['tss']-2000, x['tss']+2000]), axis=1)
    bed['exons'] = bed.apply(lambda x: tuple([(int(i), int(j)) for i, j in zip(x['exonStarts'].strip(',').split(','), x['exonEnds'].strip(',').split(','))]), axis=1)
    ### exon length
    bed['length'] = bed.apply(lambda x: sum(list(map(lambda i: (i[1]-i[0])/1000.0, x['exons']))), axis=1)
    bed['uid'] = bed.apply(lambda x: "%s@%s@%s" % (x['name2'], x['start'], x['end']), axis=1)
    bed = bed.drop_duplicates(subset='uid', keep="first")

    genes_info = []
    genes_list = []
    for irow, x in bed.iterrows():
        genes_info.append([x['chrom'], x['start'], x['end'], x['tss'], x['promoter'], x['exons'], x['length'], 1, x['uid']])
        genes_list.append(x['uid'])
    ### [chrom_0, start_1, end_2, tss_3, promoter_4, exons_5, length_6, 1_7, uid_8]
    return genes_info, genes_list

# def RP_AddExon(peaks_info, genes_info_full, genes_info_tss, decay):
#     """Multiple processing function to calculate regulation potential."""

#     def Sg(x): return 2**(-x)
#     def checkInclude(x, y): return all([x >= y[0], x <= y[1]])
#     gene_distance = 15 * decay
#     genes_peaks_score_array = scipy.sparse.dok_matrix((len(genes_info_full), len(peaks_info)), dtype=np.float64)

#     w = genes_info_tss + peaks_info
#     A = {}

#     w.sort()
#     for elem in w:
#         if elem[-3] == 1:
#             A[elem[-1]] = elem
#         else:
#             dlist = []
#             for gene_name in list(A.keys()):
#                 g = A[gene_name]
#                 tmp_distance = elem[1] - g[1]
#                 if all([g[0] == elem[0], any(list(map(checkInclude, [elem[1]]*len(g[5]), list(g[5]))))]):
#                     genes_peaks_score_array[gene_name, elem[-1]] = 1.0 / g[-4]
#                 elif all([g[0] == elem[0], tmp_distance <= gene_distance]):
#                     genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
#                 else:
#                     dlist.append(gene_name)
#             for gene_name in dlist:
#                 del A[gene_name]

#     w.reverse()
#     for elem in w:
#         if elem[-3] == 1:
#             A[elem[-1]] = elem
#         else:
#             dlist = []
#             for gene_name in list(A.keys()):
#                 g = A[gene_name]
#                 tmp_distance = g[1] - elem[1]
#                 if all([g[0] == elem[0], any(list(map(checkInclude, [elem[1]]*len(g[5]), list(g[5]))))]):
#                     genes_peaks_score_array[gene_name, elem[-1]] = 1.0 / g[-4]
#                 if all([g[0] == elem[0], tmp_distance <= gene_distance]):
#                     genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
#                 else:
#                     dlist.append(gene_name)
#             for gene_name in dlist:
#                 del A[gene_name]

#     return genes_peaks_score_array


def RP_AddExonRemovePromoter(peaks_info, genes_info_full, genes_info_tss, decay):
    """Multiple processing function to calculate regulation potential."""

    def Sg(x): return 2**(-x)
    def checkInclude(x, y): return all([x >= y[0], x <= y[1]])
    gene_distance = 15 * decay
    genes_peaks_score_array = scipy.sparse.dok_matrix((len(genes_info_full), len(peaks_info)), dtype=np.float64)
    peaks_info_inbody = []
    peaks_info_outbody = []

    w = genes_info_full + peaks_info
    A = {}

    w.sort()
#     print(w[:100])
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                ### NOTE: main change here
                ### if peak center in the gene area
                if all([g[0] == elem[0], elem[1] >= g[1], elem[1] <= g[2]]):
                    ### if peak center in the exons
                    if any(list(map(checkInclude, [elem[1]]*len(g[5]), list(g[5])))):
                        genes_peaks_score_array[gene_name, elem[-1]] = 1.0 / g[-4]
                        peaks_info_inbody.append(elem)
                    ### if peak cencer in the promoter
                    elif checkInclude(elem[1], g[4]):
                        tmp_distance = abs(elem[1]-g[3])
                        genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                        peaks_info_inbody.append(elem)
                    ### intron regions
                    else:
                        continue
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]

    ### remove genes in promoters and exons
    peaks_info_set = [tuple(i) for i in peaks_info]
    peaks_info_inbody_set = [tuple(i) for i in peaks_info_inbody]
    peaks_info_outbody_set = list(set(peaks_info_set)-set(peaks_info_inbody_set))
    peaks_info_outbody = [list(i) for i in peaks_info_outbody_set]

    print("peaks number: ", len(peaks_info_set))
    print("peaks number in gene promoters and exons: ", len(set(peaks_info_inbody_set)))
    print("peaks number out gene promoters and exons:", len(peaks_info_outbody_set))

    w = genes_info_tss + peaks_info_outbody
    A = {}

    w.sort()
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = elem[1] - g[1]
                if all([g[0] == elem[0], tmp_distance <= gene_distance]):
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]

    w.reverse()
    for elem in w:
        if elem[-3] == 1:
            A[elem[-1]] = elem
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = g[1] - elem[1]
                if all([g[0] == elem[0], tmp_distance <= gene_distance]):
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
                else:
                    dlist.append(gene_name)
            for gene_name in dlist:
                del A[gene_name]

    return genes_peaks_score_array


def RP_Simple(peaks_info, genes_info, decay):
    """Multiple processing function to calculate regulation potential."""

    Sg = lambda x: 2**(-x)
    gene_distance = 15 * decay
    genes_peaks_score_array = scipy.sparse.dok_matrix((len(genes_info), len(peaks_info)), dtype=np.float64)

    w = genes_info + peaks_info

    A = {}

    w.sort()
    for elem in w:
        if elem[2] == 1:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = abs(elem[1] - g[1])
                if (g[0] != elem[0]) or (tmp_distance > gene_distance):
                    dlist.append(gene_name)
                else:
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
            for gene_name in dlist:
                del A[gene_name]

    w.reverse()
    for elem in w:
        if elem[2] == 1:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                tmp_distance = abs(g[1] - elem[1])
                if (g[0] != elem[0]) or (tmp_distance > gene_distance):
                    dlist.append(gene_name)
                else:
                    genes_peaks_score_array[gene_name, elem[-1]] = Sg(tmp_distance / decay)
            for gene_name in dlist:
                del A[gene_name]

    return genes_peaks_score_array


def count_to_gene_by_RP(input_adata, decay=10000, model = 'simple', refgene_path=''):
    cells_list = input_adata.obs.index.tolist()
    peaks_list = input_adata.var.index.tolist()

    genes_info, genes_list = geneInfoSimple(refgene_path)
    peaks_info = []
    for ipeak, peak in enumerate(peaks_list):
        peaks_tmp = peak.rsplit("_", maxsplit=2)
        peaks_info.append([peaks_tmp[0][3:], (int(peaks_tmp[1]) + int(peaks_tmp[2])) / 2.0, 0, ipeak]) # chr, center, 0, ipeak

    if decay == 'auto':
        prop = propCalc(genes_info, peaks_info)
        if prop >= 0.2:
            decay = 1000
        else:
            decay = 10000

    if model == 'simple':
        genes_peaks_score_dok = RP_Simple(peaks_info, genes_info, decay)
    elif model == 'enhanced':
        genes_info, genes_list = ExtractGeneInfo(refgene_path)
        genes_info_tss = list()
        genes_info_full = list() ### [chrom, tss, start, end, 1, unique_id]
        for igene in range(len(genes_info)):
            tmp_gene = genes_info[igene]
            genes_info_full.append(tmp_gene + [igene])
            genes_info_tss.append([tmp_gene[0], tmp_gene[3], tmp_gene[1], tmp_gene[2]] + tmp_gene[4:] + [igene])
            ### add index at the end of gene symbol
        # genes = list(set([i.split("@")[0] for i in genes_list]))
        peaks_info = []
        for ipeak, peak in enumerate(peaks_list):
            peaks_tmp = peak.rsplit("_", maxsplit=2)
            peaks_info.append([peaks_tmp[0], (int(peaks_tmp[1])+int(peaks_tmp[2]))/2.0, int(peaks_tmp[1]), int(peaks_tmp[2]), 0, peak, ipeak])
        genes_peaks_score_dok = RP_AddExonRemovePromoter(peaks_info, genes_info_full, genes_info_tss, decay)
    else:
        sys.exit(100)
    
    genes_peaks_score_csr = genes_peaks_score_dok.tocsr()
    genes_cells_score_csr = genes_peaks_score_csr.dot(input_adata.X.T)

    score_cells_dict = {}
    score_cells_sum_dict = {}

    for igene, gene in enumerate(genes_list):
        score_cells_dict[gene] = igene
        score_cells_sum_dict[gene] = genes_cells_score_csr[igene, :].sum()

    score_cells_dict_dedup = {}
    score_cells_dict_max = {}
    genes = list(set([i.split("@")[0] for i in genes_list]))
    for gene in genes:
        score_cells_dict_max[gene] = float("-inf")

    for gene in genes_list:
        symbol = gene.split("@")[0]
        if score_cells_sum_dict[gene] > score_cells_dict_max[symbol]:
            score_cells_dict_dedup[symbol] = score_cells_dict[gene]
            score_cells_dict_max[symbol] = score_cells_sum_dict[gene]
    gene_symbol = sorted(score_cells_dict_dedup.keys())
    matrix_row = []
    for gene in gene_symbol:
        matrix_row.append(score_cells_dict_dedup[gene])

    score_cells_matrix = genes_cells_score_csr[matrix_row, :]

    RP_adata = ad.AnnData(score_cells_matrix.T, obs=pd.DataFrame(index=cells_list), var=pd.DataFrame(index=gene_symbol))
    return RP_adata


def run_target(args):
    feature_matrix_path = args.feature_matrix
    species = args.species
    output = args.output
    decay = args.decay
    model = args.model

    if feature_matrix_path.endswith('.h5'):
        input_mat_adata = sc.read_10x_h5(feature_matrix_path, gex_only=False)
    elif feature_matrix_path.endswith('.h5ad'):
        input_mat_adata = sc.read_h5ad(feature_matrix_path)
    else:
        try:
            input_mat_adata = sc.read_10x_mtx(feature_matrix_path, gex_only=False)
        except KeyError:
            input_mat_adata = sc.read_10x_mtx(feature_matrix_path, var_names='gene_ids', gex_only=False)


    CONFIG, CONFIG_PATH = read_config()

    if species == 'hs':
        refgene_path = os.path.join(CONFIG_PATH, 'GRCh38_refgenes.txt')
    elif species == 'mm':
        refgene_path = os.path.join(CONFIG_PATH, 'GRCm38_refgenes.txt')
    else:
        pass

    if decay != 'auto':
        decay = int(decay)

    rp_adata = count_to_gene_by_RP(input_adata=input_mat_adata,
                                   decay=decay,
                                   model=model,
                                   refgene_path=refgene_path
                                )
    rp_adata.write_h5ad(output)
    return
