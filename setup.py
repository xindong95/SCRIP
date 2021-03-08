#!/usr/bin/env python


try:
    from setuptools import setup, find_packages
except ImportError:
    print("Can not load setuptools!")


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
            'anndata',
            'Bio',
            'pyranges',
            'pybedtools',
            'matplotlib',
            'seaborn',
            'sklearn',
            'scipy',
            'scanpy'
        ],
        entry_points={
            'console_scripts': [
                'SCRIPT=SCRIPT.start:main'
            ]
        },


    )


if __name__ == "__main__":
    main()