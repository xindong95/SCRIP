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
obj <- readRDS('/mnt/Storage2/home/dongxin/Projects/scATAC/10X_ATAC/analysis/atac_pbmc_10k_v1_S1_L001/Result/Analysis/atac_pbmc_10k_v1_S1_L001_scATAC_Object.rds')
Seurat::as.SingleCellExperiment(obj$ATAC)
'''

#rscript = "print(1)"
with localconverter(anndata2ri.converter):
    adata = r(rscript)

print(adata)
