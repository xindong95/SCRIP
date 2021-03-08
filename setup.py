#!/usr/bin/env python


def main():
    setup(
        name='SCRIPT',
        version='0.0.210308',
        author='Xin Dong',
        author_email='xindong9511@gmail.com',
        description='A package for single cell ATAC-seq analysis',
        packages=['SCRIPT', 'SCRIPT.enrich', 'SCRIPT.utils'],
        install_requires=[
            'numpy',
            'pandas',
            'cython==0.29.22',
            'anndata',
            'Bio',
            'pyranges==0.0.95',
            'pybedtools==0.8.1',
            'matplotlib',
            'seaborn',
            'sklearn',
            'scipy',
            'scanpy==1.7.1',
        ],
        python_requires='>=3.6.*, !=3.9.*',
        entry_points={
            'console_scripts': [
                'SCRIPT=SCRIPT.start:main'
            ]
        },


    )


if __name__ == "__main__":
    try:
        from setuptools import setup, find_packages
        main()
    except ImportError:
        print("Can not load setuptools!")