# flashpca

flashpca performs fast principal component analysis (PCA) of single nucleotide
polymorphism (SNP) data, similar to smartpca from EIGENSOFT
(http://www.hsph.harvard.edu/alkes-price/software/) and shellfish
(https://github.com/dandavison/shellfish). flashpca is based on the randomized
PCA algorithm (Alg. 3) of Halko et al. 2011 (http://arxiv.org/abs/1007.5510).

Main features:

* Fast: PCA of 15,000 individuals over 43,000 SNPs in &lt;10 min
 (multi-threaded)
* Natively reads PLINK bed/bim/fam files
* Easy to use
* Two variants: the original high-memory version, and a slightly slower
   version that uses less RAM (proportional to data size), useful for large datasets with many samples
* R version [available](#R)
* Experimental: [kernel PCA](#kpca), [sparse CCA](#scca)

## Contact

Gad Abraham, gad.abraham@unimelb.edu.au

## Citation
G. Abraham and M. Inouye, Fast Principal Component Analysis of Large-Scale
Genome-Wide Data, PLos ONE 9(4): e93766. [doi:10.1371/journal.pone.0093766](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0093766)

(preprint: http://biorxiv.org/content/early/2014/03/11/002238)

## License
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Copyright (C) 2014 Gad Abraham. All rights reserved.

Portions of this code are based on SparSNP
(https://github.com/gabraham/SparSNP), Copyright (C) 2011-2012 Gad Abraham
and National ICT Australia (http://www.nicta.com.au).

## Download statically linked version (stable versions only)

* We recommend compiling from source for best performance.
* To get the devel version, you'll need to compile yourself

See [Releases](https://github.com/gabraham/flashpca/releases) for statically-linked version for Linux x86-64 &ge; 2.6.15

### System requirements
* 64-bit linux
* Using the original version (`--mem high`, the default), large datasets will require
   large amounts of RAM, e.g. for 15,000 individuals (43K SNPs) you'll need
   about 14GB RAM, for 150,000 individuals (43K SNPs) you'll need about 145GB
   RAM (estimated using https://github.com/jhclark/memusg), roughly 8 &times;
   min(n <sup>2</sup> &times; p + n &times; p, p <sup>2</sup> &times; n + n &times; p\)
   bytes, where p is number of SNPs and n is number of individuals.
   Using `--mem low`, you will only need about 8 &times; n &times; p bytes of
   RAM.

## Building from source

To get the latest version:
   ```
   git clone git://github.com/gabraham/flashpca
   ```

### Requirements

On Linux:

* 64-bit OS
* g++ compiler
* Eigen (http://eigen.tuxfamily.org), v3.2 or higher
   (if you get a compile error ``error: no match for 'operator/' in '1 / ((Eigen::MatrixBase...`` you'll need a more recent Eigen)
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
* Set CXX to whatever g++ version you're using before calling make, e.g.:
```
export CXX=/usr/local/bin/g++-4.7
```

### To install

The [Makefile](Makefile) contains three variables that need to be set according to where you have installed the Eigen
headers and Boost headers and libraries on your system. The default values for these are: 
   ```
   EIGEN_INC=/usr/local/include/eigen
   BOOST_INC=/usr/local/include/boost
   BOOST_LIB=/usr/local/lib
   ```
   
 If your system has these libraries and header files in those locations, you can simply run make:
   ```
   cd flashpca
   make all
   ```
   
 If not, you can override their values on the make command line. For example, if you have the Eigen source in `/opt/eigen-3.2.5` and Boost 1.59.0 installed into `/opt/boost-1.59.0`, you could run: 
   ```
   cd flashpca
   make all EIGEN_INC=/opt/eigen-3.2.5 BOOST_INC=/opt/boost-1.59.0/include BOOST_LIB=/opt/boost-1.59.0/lib
   ```
 
Note: the compilation process will first look for a local directory named
Eigen. It should contain the file signature_of_eigen3_matrix_library. Next,
it will look for the directory /usr/include/eigen3 (Debian/Ubuntu location
for Eigen), although those available through apt-get tend to be older versions.

## Quick start

First thin the data by LD (highly recommend
[plink2](https://www.cog-genomics.org/plink2) for this):
   ```
   plink --bfile data --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions.txt
   plink --bfile data --extract plink.prune.in --make-bed --out data_pruned
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

Eigensoft-scaling of genotypes (default):
   ```
   ./flashpca --stand binom ...
   ```

To use genotype centering (compatible with R prcomp):
   ```
   ./flashpca --stand center ...
   ```

To use the low-memory version:
   ```
   ./flashpca --mem low ...
   ```

To append a custom suffix '_mysuffix.txt' to all output files:
   ```
   ./flashpca --suffix _mysuffix.txt ...
   ```

To see all options
   ```
   ./flashpca --help 
   ```

## Output

flashpca produces the following files:

* `eigenvectors.txt`: the top k eigenvectors of the covariance
   X X<sup>T</sup> / (n - 1), same as matrix U from the SVD of the genotype matrix
   X=UDV<sup>T</sup>.
* `pcs.txt`: the top k principal components (the projection of the data on the
eigenvectors, scaled by the eigenvalues,  same as XV (or UD). This is the file
you will want to plot the PCA plot from.
* `eigenvalues.txt`: the top k eigenvalues of X X<sup>T</sup> / (n - 1). These are the
    square of the singular values D (square of sdev from prcomp).
* `pve.txt`: the proportion of total variance explained by *each of the top k*
   eigenvectors (the total variance is given by the trace of the covariance
   matrix X X<sup>T</sup> / (n - 1), which is the same as the sum of all eigenvalues).
   To get the cumulative variance explained, simply
   do the cumulative sum of the variances (`cumsum` in R).

## Warning

You must perform quality control using PLINK (at least filter using --geno, --mind,
--maf, --hwe) before running flashpca on your data. You will likely get
spurious results otherwise.

## Experimental features

### <a name="kpca"></a>Kernel PCA

* flashpca now experimentally supports low-rank [**kernel
PCA**](http://en.wikipedia.org/wiki/Kernel_principal_component_analysis) using an RBF
(Gaussian) kernel K(x, x') = exp(-||x - x'||<sub>2</sub><sup>2</sup> / sigma<sup>2</sup>) (specify using `--kernel
rbf`).
* The kernel is double-centred.
* The default kernel parameter sigma is the median of the pairwise Euclidean distances of a random subset
   of samples (controlled by `--rbfsample`, default=min(1000, n)), and can also
   be specified using `--sigma`.
* The rest of the options are the same as for standard PCA.
* Currently, the proportion of variation explained is not computed for kPCA.

### <a name="scca"></a>Sparse Canonical Correlation Analysis (SCCA)

* flashpca now experimentally supports sparse CCA
   ([Parkhomenko 2009](http://dx.doi.org/10.2202/1544-6115.1406),
   [Witten 2009](http://dx.doi.org/10.1093/biostatistics/kxp008)),
   between SNPs and multivariate phenotypes.
* The phenotype file is the same as PLINK phenotype file:
   `FID, IID, pheno1, pheno2, pheno3, ...`
   except that there must be no header line. The phenotype file *must be in the same order as
   the FAM file*.
* The L1 penalty for the SNPs is `--lambda1` and for the phenotypes is
 `--lambda2`.
* Two versions: low memory (`--mem low`, roughly 8 x \#Samples x (\#SNPs + \#Phenotypes))
   and high memory (`--mem high`, roughly 8bytes &times; \#SNPs &times; \#Phenotypes).

#### Quick example
   ```
   ./flashpca --scca --bfile data --pheno pheno.txt \
   --lambda1 1e-3 --lambda2 1e-2 --ndim 10 --numthreads 8
   ```

* The file eigenvectorsX.txt are the left eigenvectors of X<sup>T</sup> Y, with size (number of
  SNPs &times; number of dimensions), and eigenvectorsY.txt are the right
  eigenvectors of X<sup>T</sup> Y, with size (number of phenotypes &\times; number of
  dimensions).

#### Example scripts to tune the penalties via split validation

We optimise the penalties by finding the values that maximise the correlation
of the canonical components cor(X U, Y V) in independent test data.

* Wrapper script [scca.sh](scca.sh) ([GNU
   parallel](http://www.gnu.org/software/parallel) is recommended)
* R code for plotting the correlations [scca_pred.R](scca_pred.R)

# <a name="R"></a> Calling flashpca from R

flashpca is now available as an independent R package.

## Requirements for building

R packages Rcpp, RcppEigen, BH, g++ compiler

To install on Mac or Linux, you can use devtools::install_github:
   ```
   library(devtools)
   install_github("gabraham/flashpca/flashpcaR")
   ```

Note: on Mac you will need a GCC/G++ compiler (e.g., from http://brew.sh), and
to set the correct compiler in ~/.R/Makevars to point to that compiler, e.g.,
   ```
   CXX=/usr/local/bin/g++-4.9 -std=c++11
   ```

(issue https://github.com/gabraham/flashpca/issues/5)

Alternatively, after cloning the git archive, install using:
   ```
   R CMD INSTALL flashpcaR
   ```

On Windows, see [Releases](https://github.com/gabraham/flashpca/releases) for a prebuilt Windows binary package.

## PCA

Example usage, assuming `X` is a 100-sample by 1000-SNP matrix in dosage
coding (0, 1, 2) (an actual matrix, not a path to PLINK data)
   ```
   dim(X)
   [1]  100 1000
   library(flashpcaR)
   r <- flashpca(X, do_loadings=TRUE, verbose=TRUE, stand="binom", ndim=10,
   nextra=100)
   ```

PLINK data can be loaded into R either by recoding the data into raw format (`recode A`) or using package [plink2R](https://github.com/gabraham/plink2R).

Output:
   * `values`: eigenvalues
   * `vectors`: eigenvectors
   * `projection`: projection of sample onto eigenvectors (X V)
   * `loadings`: SNP loadings, if using a linear kernel

## Sparse CCA

Sparse CCA of matrices X and Y, with 5 components, penalties lambda1=0.1 and lambda2=0.1:

   ```
   dim(X)
   [1]  100 1000
   dim(Y)
   [1]  100 50
   r <- scca(X, Y, ndim=5, lambda1=0.1, lambda2=0.1)
   ```

## LD-pruned HapMap3 example data

See the [HapMap3](HapMap3) directory

## Changelog (stable versions only)

See [CHANGELOG.txt](CHANGELOG.txt)

