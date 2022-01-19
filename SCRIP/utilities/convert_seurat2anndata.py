#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   convert_seurat2anndata.py
@Time    :   2021/04/16 12:35:25
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import scanpy as sc
# import os
# import sys
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import sklearn
import anndata2ri
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
# #Loading the rpy2 extension enables cell magic to be used
# #This runs R code in jupyter notebook cells
# %load_ext rpy2.ipython

# sc.settings.verbosity = 3
# sc.logging.print_versions()

rscript = '''
obj <- readRDS('/path/to/file.rds')
Seurat::as.SingleCellExperiment(obj$ATAC)
'''

#rscript = "print(1)"
with localconverter(anndata2ri.converter):
    adata = r(rscript)

print(adata)
