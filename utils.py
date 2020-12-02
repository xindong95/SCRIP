import os
import sys
import re
import pickle
import random
import subprocess
import time
import threading
import shutil
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta
from multiprocessing import Process, Pool

import numpy as np
import pandas as pd
import anndata as ad
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import sklearn
import scipy as sp
import scanpy as sc
import episcanpy.api as epi

# from giggle import Giggle



sc.settings.verbosity = 3
# sc.logging.print_header()
sc.logging.print_versions()


plt.rcParams.update({
    'font.size' : 15,
    'figure.figsize': [8.0, 8.0],
    'font.style' : 'normal',
    'font.weight':'bold',
    'figure.titleweight': 'bold',
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'axes.spines.right': False,
    'axes.spines.top': False,
})

N = 256
vals = np.ones((N, 4))
vals[:, 0] = np.linspace(220/256, 34/256, N)
vals[:, 1] = np.linspace(220/256, 7/256, N)
vals[:, 2] = np.linspace(220/256, 141/256, N)
regulation_cmp = mpl.colors.ListedColormap(vals)


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

