flashpca
========

flashpca performs fast principal component analysis (PCA) of single nucleotide
polymorphism (SNP) data, similar to smartpca from EIGENSOFT
(http://www.hsph.harvard.edu/alkes-price/software/) and shellfish
(https://github.com/dandavison/shellfish). flashpca is based on the randomized
PCA algorithm of Halko et al. 2011.

Main features:

* Fast: PCA of 15,000 individuals over 43,000 SNPs in &lt;10 min
* Natively reads PLINK bed/bim/fam files
* Easy to use

Contact
-------

Gad Abraham, gad.abraham@unimelb.edu.au

Citation
--------
``Fast Principal Component Analysis of Large-Scale Genome-Wide Data''

License
-------

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Copyright (C) 2014 Gad Abraham and National ICT Australia
All rights reserved.

Portions of this code are based on SparSNP
(https://github.com/gabraham/SparSNP), Copyright (C) 2011-2012 Gad Abraham
and National ICT Australia (http://www.NICTA.com.au).

Requirements
------------

   * Linux OS (might work on Mac OSX), 64-bit
   * g++ compiler
   * Eigen (http://eigen.tuxfamily.org)
   * Boost (http://www.boost.org), specifically boost_system-mt,
      boost_iostreams-mt, boost_filesystem-mt, boost_program_options
   * libgomp for openmp support
   * Recommended: plink2 (https://www.cog-genomics.org/plink2) for SNP
      thinning

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

To run on a PLINK fileset named data (data.bed/data.bim/data.fam):
   ```
   ./flashpca --bfile data
   ```

To run in multi-threaded mode with 8 threads:
   ```
   ./flashpca --bfile data --numthreads 8

For more options
   ```
   ./flashpca --help
   ```


