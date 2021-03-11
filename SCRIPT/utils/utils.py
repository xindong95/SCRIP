import os
import sys
# import re
# import pickle
# import random
# import subprocess
import time
# import threading
# import shutil
# from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, wait, ALL_COMPLETED
# from datetime import datetime, timedelta
# from multiprocessing import Process, Pool

import ruamel.yaml
# import numpy as np
# import pandas as pd
# import anndata as ad
# # import h5py
# import Bio
# from Bio import motifs
# import pysam
# import pyranges
# import pybedtools
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# import seaborn as sns
# import sklearn
# from sklearn.preprocessing import quantile_transform
# import scipy as sp
# import scanpy as sc
# from adjustText import adjust_text
# import episcanpy

# from giggle import Giggle



# sc.settings.verbosity = 3
# # sc.logging.print_header()
# sc.logging.print_versions()


# plt.rcParams.update({
#     'font.size' : 15,
#     'figure.figsize': [8.0, 8.0],
#     'font.style' : 'normal',
#     'font.weight':'bold',
#     'figure.titleweight': 'bold',
#     'axes.titleweight': 'bold',
#     'axes.labelweight': 'bold',
#     'axes.spines.right': False,
#     'axes.spines.top': False,
# })

# N = 256
# vals = np.ones((N, 4))
# vals[:, 0] = np.linspace(220/256, 34/256, N)
# vals[:, 1] = np.linspace(220/256, 7/256, N)
# vals[:, 2] = np.linspace(220/256, 141/256, N)
# regulation_cmp = mpl.colors.ListedColormap(vals)


def add_time(func):
    def wrapper(*args, **kw):
        sys.stdout.write('INFO {time} '.format(time=time.strftime("%Y-%m-%d %H:%M:%S")))
        return func(*args, **kw)
    return wrapper

@add_time
def print_log(string, end="\n"):
    print(string, end=end)

def excute_info(start_info=None, end_info=None):
    def call(func):
        def wrapper(*args, **kwargs):
            if start_info:
                print_log(start_info)
            ret = func(*args, **kwargs)
            if end_info:
                print_log(end_info)
            return ret
        return wrapper
    return call

def read_SingleCellExperiment_rds(input_RDS):
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()
    rscript = 'readRDS("{RDS_file_path}")'.format(RDS_file_path = input_RDS)
    adata = r(rscript)
    return adata

def read_config():

    yaml = ruamel.yaml.YAML()
    CONFIG_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','conf','config.yml'))
    with open(CONFIG_PATH,'r') as config_file:
        CONFIG = yaml.load(config_file.read())
    return CONFIG, CONFIG_PATH

