flashpca
========

flashpca performs fast principal component analysis (PCA) of single nucleotide
polymorphism (SNP) data, similar to smartpca from EIGENSOFT
(http://www.hsph.harvard.edu/alkes-price/software/) and shellfish
(https://github.com/dandavison/shellfish). flashpca is based on the randomized
PCA algorithm of Halko et al. 2011 (http://arxiv.org/abs/1007.5510).

Main features:

* Fast: PCA of 15,000 individuals over 43,000 SNPs in &lt;10 min
 (multi-threaded)
* Natively reads PLINK bed/bim/fam files
* Easy to use

Contact
-------

Gad Abraham, gad.abraham@unimelb.edu.au

Citation
--------
G. Abraham and M. Inouye, Fast Principal Component Analysis of Large-Scale
Genome-Wide Data, submitted (2014)
(http://biorxiv.org/content/early/2014/01/30/002238)

License
-------

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Copyright (C) 2014 Gad Abraham. All rights reserved.

Portions of this code are based on SparSNP
(https://github.com/gabraham/SparSNP), Copyright (C) 2011-2012 Gad Abraham
and National ICT Australia (http://www.NICTA.com.au).

Download statically linked version
----------------------------------

Note: we recommend compiling from source for best performance.

[flashpca_x86-64.gz](flashpca_x86-64.gz) for linux 2.6.15 and higher,
   gunzip before use.
   
~~Note: currently OpenMP doesn't work correctly for statically linked version,
so will run at single-thread speed. Compile from source to get the
multiple-thread support.~~ (update: fixed, solution at
http://stackoverflow.com/questions/11165338/segfault-when-using-static-and-fopenmp-in-gfortran)

System requirements
-------------------
* 64-bit linux.
* For large datasets you'll need large amounts of RAM, e.g. for 15,000
   individuals (43K SNPs) you'll need about 14Gb RAM, for 150,000 individuals
   (43K SNPs) you'll need about 145Gb RAM (estimated using
   https://github.com/jhclark/memusg)
   (we are working on a version which will use less RAM).

Requirements for building from source
-------------------------------------

   On Linux:
   * 64-bit OS
   * g++ compiler
   * Eigen (http://eigen.tuxfamily.org), v3.2 or higher
      (if you get a compile error ``error: no match for ‘operator/’ in ‘1 / ((Eigen::MatrixBase...`` you'll need a more recent Eigen)
   * Boost (http://www.boost.org/), specifically boost_system/boost_system-mt,
      boost_iostreams/boost_iostreams-mt,
      boost_filesystem/boost_filesystem-mt,
      boost_program_options/boost_program_options-mt.
   * libgomp for openmp support
   * Recommended: plink2 (https://www.cog-genomics.org/plink2) for SNP
      thinning

   On Mac:
   * Homebrew (http://brew.sh) to install gcc/g++ and boost
   * Eigen, as above
   * Set CXX to whatever g++ version you're using before calling make:
   ```
   CXX=/usr/local/bin/g++-4.7 make all
   ```

Quick start
-----------

To get the latest version:
   ```
   git clone git://github.com/gabraham/flashpca
   ```

To install:
   ```
   cd flashpca
   make
   ```

Note: the compilation process will first look for a local directory named
Eigen. It should contain the file signature_of_eigen3_matrix_library. Next,
it will look for the directory /usr/include/eigen3 (Debian/Ubuntu location
for Eigen), although those available through apt-get tend to be older versions.

First thin the data by LD (highly recommend plink2 for this):
   ```
   plink --bfile data --indep-pairwise 1000 50 0.05
   plink --bfile data --extract plink.prune.in --make-bed --out data_pruned --exclude exclusion_regions.txt
   ```
where [exclusion_regions.txt](exclusion_regions.txt) contains:
   ```
   5 44000000 51500000 r1
   6 25000000 33500000 r2
   8 8000000 12000000 r3
   11 45000000 57000000 r4
   ```
(You may need to change the --indep-pairwise parameters to get a suitable
number of SNPs for you dataset, 10,000-50,000 is usually enough.)

To run on the pruned dataset:
   ```
   ./flashpca --bfile data_pruned
   ```

We highly recommend using multi-threading, to run in multi-threaded mode with 8 threads:
   ```
   ./flashpca --bfile data_pruned --numthreads 8
   ```

For more options:
   ```
   ./flashpca --help
   ```

LD-pruned HapMap3 example data
------------------------------

See the [HapMap3](HapMap3) directory

   
