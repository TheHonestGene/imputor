**********************************
Imputor: Fast imputation library
**********************************

Imputor is a very simple library for imputing missing variants 
in a genotype dataset (i.e. 23andMe)

Users can run Imputor as either a standalone command line tool 
or import it into their python library


Requirements
***********

The imputation uses an LD approach and in order to calculate the LD a reference dataset 
is required (usually the 1000 Genomes HapMap dataset). This HapMap dataset has to be provided as an HDF5 file and should only contain bi-allelic SNPs. 
A version can be downloaded here.


How-To
***********

There are 6 steps to finish in the imputation pipeline.
Some of the steps only have to be done once or once for each genotyping plattform and dataset respectively
Each step maps to a subcommand of the command line script. 
To get information about the various subcommands the user can run: 

.. code:: bash

      imputor -h

Step 1 - Prepare a filtered 1000 genomes HapMap dataset
===============================================================
Currently risk prediction only works reliably with European indviduals. 
So we first need to remove all non-european individuals from the 1000 HapMap dataset 
and also remove the monomorphic SNPs. This can be done with the prepare command. 
This only has to be done once

.. code:: bash

      imputor prepare 1000_genomes.hdf5 1001_genomes_unrelated.hdf5
      

Step 2 - Convert genoptype from text/csv format to HDF5 format
===============================================================
The genotypes from 23andMe and others are usually provided in text or csv format. 
We first need to convert them into an HDF5 format for the downstream analysis:

.. code:: bash

      imputor parse genotype.csv genoptype.hdf5
      
      

Step 3 - Create a nucleotide encoding map
===============================================================
The genotypes have to be converted from the NT form to a binary form. 
To speed up this conversion a nucleotide map is constructed. 
The map is constructed from the HapMap dataset that was created in Step 1. 
This has to be done once for each genotype version and platform respectively. 

.. code:: bash

      imputor nt_map 1000_genomes_unrelated.hdf5 genoptype.hdf5 nt_map.pickled

Step 4 - Create the LD files needed for imputation
===============================================================
For the imputation LD files have to be created. This has to be done once per genotype version 
and platform respectively. The user can specify an optional window_size for the LD calculation

.. code:: bash

      imputor calc_ld 1000_genomes_unrelated.hdf5 nt_map.pickled FOLDER_FOR_LD_FILES


Step 5 - Convert genotype from NT form to binary form
===============================================================
In the last step before imputation the genotype file has to be converted from NT form to binary form. 
This is done using the the nucleotide encoding map that was generated in Step 3::

    $ imputor convert genotype.hdf5 genotype_encoded.hdf5 nt_map.pickled


Step 6 - Imputation
=========================
Using the LD files the missing SNPs can be imputed for an existing genotype::
    
    $ imputor impute genotype_encoded.hdf5 FOLDER_FOR_LD_FILES



Installation
--------------

Of course, the recommended installation method is pip::

    $ pip install imputor

Thank You
-----------

Thanks for checking this library out! We hope you find it useful.

Of course, there's always room for improvement. Feel free to `open an issue <https://github.com/TheHonestGene/imputor/issues>`_ so we can make Imputor better.


