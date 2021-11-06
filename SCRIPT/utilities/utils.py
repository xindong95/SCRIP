#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   utils.py
@Time    :   2021/04/16 12:29:39
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import os
import sys
import pickle
import time

import ruamel.yaml
import pandas as pd
from scipy import io


def add_time(func):
    def wrapper(*args, **kw):
        time_now = time.strftime("%Y-%m-%d %H:%M:%S")
        sys.stdout.write(f'INFO {time_now} ')
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
    rscript = f'readRDS("{input_RDS}")'
    adata = r(rscript)
    return adata

def read_config():
    yaml = ruamel.yaml.YAML()
    CONFIG_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','conf'))
    CONFIG_YML_PATH = os.path.join(CONFIG_PATH,'config.yml')
    with open(CONFIG_YML_PATH, 'r') as config_file:
        CONFIG = yaml.load(config_file.read())
    return CONFIG, CONFIG_PATH

def store_to_pickle(obj, path):
    with open(os.path.join(path), 'wb+') as f:
        pickle.dump(obj, f)


def read_pickle(path):
    with open(os.path.join(path), 'rb') as f:
        ret = pickle.load(f)
    return ret

def safe_makedirs(path):
    if not os.path.exists(path):
        os.makedirs(path)


def write_to_mtx(data, path):
    if not os.path.exists(path):
        os.makedirs(path)
    pd.DataFrame(data.var.index, data.var.index).to_csv(os.path.join(path, "genes.tsv"), sep="\t", header=False)
    pd.DataFrame(data.obs.index).to_csv(os.path.join(path, "barcodes.tsv"), sep = "\t", index=False, header=False)
    data.obs.to_csv(os.path.join(path, "metadata.tsv"), sep = "\t", index=False, header=False)
    io.mmwrite(os.path.join(path, "matrix.mtx"), data.X.T)
