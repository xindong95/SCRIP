#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   utils.py
@Time    :   2021/04/16 12:34:42
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import os
import sys
from datetime import datetime, timedelta
import json
from SCRIP.Constants import SCRIP_VERSION
from SCRIP.utilities.utils import store_to_pickle, read_pickle


def time_estimate(cell_number, core, chip_factor_number=4499):
#   generate bed: 0.5 second per cell,  search index: 10 seconds per cell, cal p: 2 seconds per factor
    # chip_process = True if reference_method in ['integration', 'both', 'chip'] else False
    # motif_process = True if reference_method in ['integration', 'both', 'motif'] else False
    seconds = int(cell_number * core / 50) # single core process, find nearest cells function
    seconds += 3 * cell_number
    # if chip_process == True:
    seconds += 10 * cell_number
    seconds += 2 * chip_factor_number
    # if motif_process == True:
    #     seconds += 10 * cell_number
    #     seconds += 2 * motif_factor_number

    seconds = seconds / core
    # now_time = datetime.now()
    future_time = (datetime.now() + timedelta(seconds=seconds)).strftime("%Y-%m-%d %H:%M:%S")
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    elapse = "{hour:0>2d}:{minute:0>2d}:{second:0>2d}".format(hour=int(h), minute=int(m), second=int(s))
    return elapse, future_time

class EnrichRunInfo(object):
    def __init__(self, file_path, parameters_dict):
        if os.path.exists(file_path):
            with open(file_path, 'r') as jf:
                info = json.load(jf)
        else:
            info = {}
            info['version'] = SCRIP_VERSION
            info['file_path'] = os.path.abspath(file_path)
            info['project_folder'] = os.path.abspath(os.path.dirname(file_path))
            info['parameters'] = parameters_dict
            info['progress'] = {

                'bed_generation': 'No',
                'bed_search': 'No',

                'dataset_overlap_df_store': 'No',
                'dataset_cell_norm_df_store': 'No',
                'dataset_score_source_df_store': 'No',
                'tf_cell_score_df_store': 'No',

                'result_store': 'No'
            }
        self.version = info['version']
        self.file_path = info['file_path']
        self.parameters = info['parameters']
        self.progress = info['progress']
        self.info = info
        if not os.path.exists(os.path.dirname(file_path)):
            os.makedirs(os.path.dirname(file_path))
        with open(file_path, 'w+') as jf:
            json.dump(self.info, jf, indent=2)

    def finish_stage(self, stage):
        self.progress[stage] = 'Finish'
        self.info['progress'][stage] = 'Finish'
        with open(self.file_path, 'w') as jf:
            json.dump(self.info, jf, indent=2)

    def reset_stage(self, stage):
        self.progress[stage] = 'No'
        self.info['progress'][stage] = 'No'
        with open(self.file_path, 'w') as jf:
            json.dump(self.info, jf, indent=2)

    def check_consist(self, parameters_dict):
        if self.version != SCRIP_VERSION:
            return False
        inconsistent_key = []
        for k in self.parameters.keys():
            try:
                if self.parameters[k] == parameters_dict[k]:
                    continue
                else:
                    inconsistent_key.append(k)
            except KeyError:
                inconsistent_key.append(k)
        if inconsistent_key.__len__() != 0:
            return False
        else:
            return True

    def safe_run_and_store(self, func, func_arg_list, store_path, stage):
        if self.info['progress'][stage] == 'No':
            ret = func(*func_arg_list)
            store_to_pickle(ret, store_path)
            self.finish_stage(stage)
        else:
            ret = read_pickle(store_path)
        return ret

    def safe_run(self, func, func_arg_list, stage):
        if self.info['progress'][stage] == 'No':
            func(*func_arg_list)
            self.finish_stage(stage)
        else:
            pass
        return
