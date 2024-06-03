#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   setup.py
@Time    :   2021/04/16 12:34:01
@Author  :   Xin Dong
@Contact :   xindong9511@gmail.com
@License :   (C)Copyright 2020-2021, XinDong
'''

import os
import sys
import subprocess


# def install_giggle():
#     prefix = os.environ['PATH'].split(':')[0]
#     cmd = f'cd refpkg/giggle; make; cp bin/giggle {prefix}/; cd ../..'
#     subprocess.run(cmd, shell=True)

def main():
    # install_giggle()
    setup(
        name='SCRIP',
        version='0.1.240527',
        author='Xin Dong',
        author_email='xindong9511@gmail.com',
        description='A package for single cell ATAC-seq analysis',
        packages=['SCRIP', 'SCRIP.enrichment', 'SCRIP.enhancement', 'SCRIP.conf', 'SCRIP.utilities', 'SCRIP.imputation', 'SCRIP.targets', 'SCRIP.index'],
        package_data={
            'SCRIP': [
                'conf/config.yml',
                'conf/GRCh38_refgenes.txt',
                'conf/GRCm38_refgenes.txt',
                ],
        },
        install_requires=[
            'numba>=0.51.2',
            'numpy>=1.18',
            'pandas>=1.3.3',
            'cython>=0.29.22',
            'ruamel.yaml',
            'anndata>=0.7.6',
            'anndata2ri>=1.0.6',
            # 'Bio',
            'pyranges>=0.0.95',
            'pybedtools>=0.8.2',
            'matplotlib>=3.1.2',
            # 'seaborn',
            'scikit-learn',
            'scipy',
            'scanpy>=1.7.1',
        ],
        python_requires='>=3.8.12, !=3.9',
        entry_points={
            'console_scripts': [
                'SCRIP=SCRIP.start:main'
            ]
        },
    )



if __name__ == "__main__":
    try:
        from setuptools import setup, find_packages
        main()
    except ImportError:
        print("Can not load setuptools!")
