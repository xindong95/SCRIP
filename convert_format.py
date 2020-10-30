import scanpy as sc
import os
import sys
import anndata as ad
import h5py


def singleCellExperiment2h5ad(input_RDS, output_h5ad):
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()
    rscript = 'readRDS("{RDS_file_path}")'.format(RDS_file_path = input_RDS)
    adata = r(rscript)
    adata.write(output_h5ad)