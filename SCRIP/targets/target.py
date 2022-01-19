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
import scipy
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from SCRIP.utilities.utils import read_config
# from SCRIP.enrichment.bed_generation import find_nearest_cells

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


def count_to_gene_by_RP(input_adata, decay=10000, refgene_path=''):
    cells_list = input_adata.obs.index.tolist()
    peaks_list = input_adata.var.index.tolist()

    genes_info, genes_list = geneInfoSimple(refgene_path)

    peaks_info = []
    for ipeak, peak in enumerate(peaks_list):
        peaks_tmp = peak.rsplit("_", maxsplit=2)
        peaks_info.append([peaks_tmp[0][3:], (int(peaks_tmp[1]) + int(peaks_tmp[2])) / 2.0, 0, ipeak])

    genes_peaks_score_dok = RP_Simple(peaks_info, genes_info, decay)

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

    if feature_matrix_path.endswith('.h5'):
        input_mat_adata = sc.read_10x_h5(feature_matrix_path, gex_only=False)
    elif feature_matrix_path.endswith('.h5ad'):
        input_mat_adata = sc.read_h5ad(feature_matrix_path)
    else:
        try:
            input_mat_adata = sc.read_10x_mtx(feature_matrix_path, gex_only=False)
        except KeyError:
            input_mat_adata = sc.read_10x_mtx(feature_matrix_path, var_names='gene_ids' ,gex_only=False)


    CONFIG, CONFIG_PATH = read_config()

    if species == 'hs':
        refgene_path = os.path.join(CONFIG_PATH, 'GRCh38_refgenes.txt')
    elif species == 'mm':
        refgene_path = os.path.join(CONFIG_PATH, 'GRCm38_refgenes.txt')
    else:
        pass
    rp_adata = count_to_gene_by_RP(input_adata=input_mat_adata,
                                   decay=decay,
                                   refgene_path=refgene_path
                                )
    rp_adata.write_h5ad(output)
    return
