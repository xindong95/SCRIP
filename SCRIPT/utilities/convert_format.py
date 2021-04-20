#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   convert_format.py
@Time    :   2021/04/16 12:35:19
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import scanpy as sc
import os
import sys
import anndata as ad
import h5py


def read_SingleCellExperiment_rds(input_RDS):
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()
    rscript = 'readRDS("{RDS_file_path}")'.format(RDS_file_path = input_RDS)
    adata = r(rscript)
    adata.var.columns = [str(i) for i in adata.var.columns]
    adata.obs.columns = [str(i) for i in adata.obs.columns]
    return adata