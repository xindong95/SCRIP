import os
import sys
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
import pickle
import random
import subprocess
import time
import threading
from datetime import datetime
from multiprocessing import Process, Pool

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



