#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   impute.py
@Time    :   2021/05/20 09:49:38
@Author  :   Xin Dong 
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import scanpy as sc
import random
import numpy as np
from SCRIPT.enrichment.bed_generation import generate_peak_list


def determine_number_of_cells_per_group(input_mat, start=1, end=70, iteration_time=30, aim_peak_number=10000, peak_confidence=None):
    peak_number = 0
    cell_pool = input_mat.obs_names.to_list()
    cell_number = start - 1 
    while peak_number <= aim_peak_number:
        cell_number += 1
        if peak_confidence == None:
            peak_confidence = np.ceil(0.2*cell_number)
        tmp_peak_number = []
        for i in range(iteration_time):
            cells = random.sample(cell_pool, cell_number)
            tmp_peak_number.append(generate_peak_list(cells, input_mat, peak_confidence).__len__())
        peak_number = sum(tmp_peak_number)/iteration_time
        print(tmp_peak_number)
    print(cell_number)