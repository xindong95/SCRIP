#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   index.py
@Time    :   2022/03/04 14:05:51
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import os
from SCRIP.utilities.utils import print_log, safe_makedirs

def run_index(args):

    input_path = os.path.abspath(args.input)
    output_path = os.path.abspath(args.output)

    safe_makedirs(output_path)
    safe_makedirs(os.path.join(input_path, 'raw_beds'))

    for i in os.listdir(input_path):
        if i.endswith('bed'):
            r_bed_path = os.path.join(input_path, i)
            aim_bed_path = os.path.join(input_path, 'raw_beds', f'{i}.gz')
            cmd = f'sort --buffer-size 2G -k1,1 -k2,2n -k3,3n {r_bed_path} | bedtools merge | bgzip -c > {aim_bed_path}'
            os.system(cmd)

    bed_gz_path = os.path.join(input_path, 'raw_beds')
    peaks_number_path = os.path.join(output_path, 'peaks_number.txt')
    cmd = f'cd {bed_gz_path}\n'
    cmd += 'ulimit -n 4096\n'
    cmd += f"""giggle index -i '*.bed.gz' -o {output_path} -f -s \n"""
    cmd += r"""for i in *.bed.gz; do echo ${i::-7}; number=`zcat $i | wc -l`; echo -e ${i::-7}'\t'$number >>"""
    cmd += f""" {peaks_number_path}; done\n"""
    cmd += f"mv {bed_gz_path} {output_path}\n"
    os.system(cmd)

    print_log("Finished generating the custom index! Please run the SCRIP config function.")

    return
