import os
import sys
from datetime import datetime, timedelta
import json
from SCRIPT.Constants import *
from SCRIPT.utilities.utils import print_log

def time_estimate(cell_number, bg_iteration_number, peak_methods, cell_number_per_group, chip_process, motif_process, core):
#   generate bed: 0.5 second per cell,  search index: 10 seconds per cell, cal p: 2 seconds per factor
    if peak_methods == "group":
        bed_number = int(cell_number/cell_number_per_group) + bg_iteration_number
        seconds = 3 * bed_number
        if chip_process == True:
            seconds += 10 * bed_number
        if motif_process == True:
            seconds += 10 * bed_number
    elif peak_methods == "nearest":
        seconds = int(cell_number * core / 50) # single core process, find nearest cells function
        bed_number = cell_number + bg_iteration_number
        seconds += 3 * bed_number
        if chip_process == True:
            seconds += 10 * bed_number
        if motif_process == True:
            seconds += 10 * bed_number
    chip_factor_number = 5069
    motif_factor_number = 1143
    if chip_process == True:
        seconds += 2 * chip_factor_number
    if motif_process == True:
        seconds += 2 * motif_factor_number
    seconds = seconds / core
    now_time = datetime.now()
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
            info['version'] = SCRIPT_VERSION
            info['file_path'] = os.path.abspath(file_path)
            info['project_folder'] = os.path.abspath(os.path.dirname(file_path))
            info['parameters'] = parameters_dict
            info['progress'] = {
                'bg_bed_generation': 'No', 
                'fg_bed_generation': 'No',
                'bg_bed_motif_search': 'No',
                'bg_bed_chip_search': 'No',
                'fg_bed_motif_search': 'No',
                'fg_bed_chip_search': 'No'
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
        self.info['progress'] = 'Finish'
        with open(file_path, 'w') as jf:
            json.dump(self.info, jf, indent=2)

    def reset_stage(self, stage):
        self.progress[stage] = 'No'
        self.info['progress'] = 'No'
        with open(file_path, 'w') as jf:
            json.dump(self.info, jf, indent=2)
    
    def check_consist(self, parameters_dict):
        if self.version != SCRIPT_VERSION:
            return False
        inconsistent_key = []
        for k in self.parameters.keys():
            if self.parameters[k] == parameters_dict[k]:
                continue
            else:
                inconsistent_key.append(k)
        if inconsistent_key.__len__() != 0:
            return False
        else:
            return True



