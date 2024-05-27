Installation
==============


System requirements
~~~~~~~~~~~~~~~~~~~

* Linux/Unix
* Python >= 3.8
* libpng12-0
* tabix

The reference files can be downloaded from `zenodo <https://zenodo.org/record/5840810>`_.  

We recommend to create an independent conda environment for SCRIP. If users do not have conda, please install Miniconda first:

.. code:: shell

   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh

Install from GitHub
~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   git clone git@github.com:wanglabtongji/SCRIP.git
   cd SCRIP
   python setup.py install

Install from pypi
~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   pip install scrip


Install GIGGLE
~~~~~~~~~~~~~~~~~~~~~

SCRIP fast searching is based on `GIGGLE <https://github.com/ryanlayer/giggle>`_. Please install GIGGLE manually first.

.. code:: shell

   git clone git@github.com:ryanlayer/giggle.git
   cd giggle
   make
   export PATH=$PATH:`pwd` # or cp bin/giggle to your environment
   cd ..

Next, validate the installation:

.. code:: shell

   giggle


It should return:

.. code:: 
   
   giggle, v0.6.3
   usage:   giggle <command> [options]
      index     Create an index
      search    Search an index

Besides manually installation, you can try `SCRIP install_giggle` too.