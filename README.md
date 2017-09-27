# FlashPCA2

FlashPCA performs fast principal component analysis (PCA) of single nucleotide
polymorphism (SNP) data, similar to smartpca from EIGENSOFT
(http://www.hsph.harvard.edu/alkes-price/software/) and shellfish
(https://github.com/dandavison/shellfish). FlashPCA is based on the
[https://github.com/yixuan/spectra/](Spectra) library.

Main features:

* Fast: partial PCA (_k_=20 dimensions) of 500,000 individuals with 100,000 SNPs in &lt;6h using 2GB RAM
* Scalable: memory requirements are bounded, scales to at least 1M individuals
* Highly accurate results 
* Natively reads PLINK bed/bim/fam files
* Easy to use; can be called entirely within R (package [flashpcaR](#flashpcaR))

## Help

Google Groups: [https://groups.google.com/forum/#!forum/flashpca-users](https://groups.google.com/forum/#!forum/flashpca-users)

## Contact

Gad Abraham, gad.abraham@unimelb.edu.au

## Citation

version &geq;2: G. Abraham, Y. Qiu, and M. Inouye, ``FlashPCA2: principal
component analysis of biobank-scale genotype datasets'', (2017) Bioinformatics
33(17): 2776-2778.
[doi:10.1093/bioinformatics/btx299](https://doi.org/10.1093/bioinformatics/btx299)
 (bioRxiv preprint [https://doi.org/10.1101/094714](https://doi.org/10.1101/094714))

version &leq;1.2.6: G. Abraham and M. Inouye, ``Fast Principal Component
Analysis of Large-Scale Genome-Wide Data'', (2016) PLOS ONE 9(4): e93766. [doi:10.1371/journal.pone.0093766](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0093766)

## License
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Copyright (C) 2014-2016 Gad Abraham. All rights reserved.

Portions of this code are based on SparSNP
(https://github.com/gabraham/SparSNP), Copyright (C) 2011-2012 Gad Abraham
and National ICT Australia (http://www.nicta.com.au).

## Download statically linked version (stable versions only, not for alpha versions)

* We recommend compiling from source for best performance.
* To get the devel version, you'll need to compile yourself

See [Releases](https://github.com/gabraham/flashpca/releases) for statically-linked version for Linux x86-64 &ge; 2.6.15

### System requirements
* 64-bit Linux or Mac

## Building from source

To get the latest version:
   ```bash
   git clone git://github.com/gabraham/flashpca
   ```

### Requirements

On Linux:

* 64-bit OS
* g++ compiler
* Eigen (http://eigen.tuxfamily.org), v3.2 or higher
   (if you get a compile error ``error: no match for 'operator/' in '1 / ((Eigen::MatrixBase...`` you'll need a more recent Eigen)
* Spectra (https://github.com/yixuan/spectra/)
* Boost (http://www.boost.org/), specifically boost_program_options/boost_program_options-mt.
* libgomp for openmp support
* Recommended: plink2 (https://www.cog-genomics.org/plink2) for SNP
   thinning

On Mac:

* Homebrew (http://brew.sh) to install boost
* Eigen, as above
* Spectra, as above
* clang C++ compiler

### To install

The [Makefile](Makefile) contains three variables that need to be set according to where you have installed the Eigen
headers and Boost headers and libraries on your system. The default values for these are: 
   ```bash
   EIGEN_INC=/usr/local/include/eigen
   BOOST_INC=/usr/local/include/boost
   BOOST_LIB=/usr/local/lib
   SPECTRA_INC=spectra
   ```
   
 If your system has these libraries and header files in those locations, you can simply run make:
   ```bash
   cd flashpca
   make all
   ```
   
 If not, you can override their values on the make command line. For example,
 if you have the Eigen source in `/opt/eigen-3.2.5`, spectra headers in
 `/opt/spectra`, and Boost 1.59.0 installed into `/opt/boost-1.59.0`, you could run: 
   ```bash
   cd flashpca
   make all EIGEN_INC=/opt/eigen-3.2.5 \
      BOOST_INC=/opt/boost-1.59.0/include \
      BOOST_LIB=/opt/boost-1.59.0/lib \
      SPECTRA_INC=/opt/spectra
   ```
 
## Quick start

First thin the data by LD (highly recommend
[plink2](https://www.cog-genomics.org/plink2) for this):
   ```bash
   plink --bfile data --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions_hg19.txt
   plink --bfile data --extract plink.prune.in --make-bed --out data_pruned
   ```
where [exclusion_regions_hg19.txt](exclusion_regions_hg19.txt) contains:
   ```
   5 44000000 51500000 r1
   6 25000000 33500000 r2
   8 8000000 12000000 r3
   11 45000000 57000000 r4
   ```
(You may need to change the --indep-pairwise parameters to get a suitable
number of SNPs for you dataset, 10,000-50,000 is usually enough.)

To run on the pruned dataset:
   ```bash
   ./flashpca --bfile data_pruned
   ```

To append a custom suffix '_mysuffix.txt' to all output files:
   ```bash
   ./flashpca --suffix _mysuffix.txt ...
   ```

To see all options
   ```bash
   ./flashpca --help 
   ```

## Output

By default, flashpca produces the following files:

* `eigenvectors.txt`: the top k eigenvectors of the covariance
   X X<sup>T</sup> / p, same as matrix U from the SVD of the genotype matrix
   X/sqrt(p)=UDV<sup>T</sup> (where p is the number of SNPs).
* `pcs.txt`: the top k principal components (the projection of the data on the
eigenvectors, scaled by the eigenvalues,  same as XV (or UD). This is the file
you will want to plot the PCA plot from.
* `eigenvalues.txt`: the top k eigenvalues of X X<sup>T</sup> / p. These are the
    square of the singular values D (square of sdev from prcomp).
* `pve.txt`: the proportion of total variance explained by *each of the top k*
   eigenvectors (the total variance is given by the trace of the covariance
   matrix X X<sup>T</sup> / p, which is the same as the sum of all eigenvalues).
   To get the cumulative variance explained, simply
   do the cumulative sum of the variances (`cumsum` in R).

## Warning

You must perform quality control using PLINK (at least filter using --geno, --mind,
--maf, --hwe) before running flashpca on your data. You will likely get
spurious results otherwise.

## Projection of new samples onto PCs

flashpca can project new samples onto existing principal components:
   ```bash
   ./flashpca --bfile newdata --project --inmeansd meansd.txt \
      --outproj projections.txt --inload loadings.txt -v
   ```

To project data, you must ensure:

* The old and new PLINK files contain _exactly_ the same SNPs and alleles (you
can use `plink --reference-allele ...` to ensure consistent allele ordering).
* You have previously run flashpca and saved the SNP loadings
(`--outload loadings.txt`) and their means and standard deviations
(`--outmeansd meansd.txt`).
* You are using the same standardisation (`--standx`) for the old and new
data, as well as the same divisor (`--div`; by default `p`). 


## Checking accuracy of results

flashpca can check how accurate a decomposition is, where accuracy is defined
as || X X<sup>T</sup> / p - U D<sup>2</sup> ||<sub>F</sub><sup>2</sup> / (n
&times; k).

This is done using

   ```bash
   ./flashpca --bfile data --check \
   --outvec eigenvectors.txt --outval eigenvalues.txt
   ```

The final mean squared error should be low (e.g., <1e-8).

## Outlier removal in PCA

Unlike EIGENSOFT/smartpca, flashpca does not remove outliers automatically
(`numoutlieriter` in EIGENSOFT). We recommend inspecting the PCA plot
manually, and if you wish to remove outliers and repeat PCA on the remaining
samples, use `plink --remove` to create a new bed/bim/fam fileset and run
flashpca on the new data.

## <a name="scca"></a>Sparse Canonical Correlation Analysis (SCCA)

* flashpca now supports sparse CCA
   ([Parkhomenko 2009](http://dx.doi.org/10.2202/1544-6115.1406),
   [Witten 2009](http://dx.doi.org/10.1093/biostatistics/kxp008)),
   between SNPs and multivariate phenotypes.
* The phenotype file is the same as PLINK phenotype file:
   `FID, IID, pheno1, pheno2, pheno3, ...`
   except that there must be no header line. The phenotype file *must be in the same order as
   the FAM file*.
* The L1 penalty for the SNPs is `--lambda1` and for the phenotypes is
 `--lambda2`.

#### Quick example
   ```bash
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

# <a name="flashpcaR"></a>flashpcaR: flashpca in R

FlashPCA can be called (almost) entirely within R.

## Installation

* See [Releases](https://github.com/gabraham/flashpca/releases) for pre-packaged
flashpcaR packages for Windows, Mac, and Linux.

* Compiling from source

   ```R
   devtools::install_github("gabraham/flashpca/flashpcaR")
   ```

   This may fail on Windows due to symbolic link issues; in that case use 
   ```bash
   git clone https://github.com/gabraham/flashpca
   cd flashpca
   R CMD build flashpcaR
   R CMD INSTALL flashpcaR_2.0.tar.gz
   ```

Depends on packages:
   
   * [Rcpp](https://cran.r-project.org/package=Rcpp) (>= 0.11.1)
   * [RcppEigen](https://cran.r-project.org/package=RcppEigen) (>= 0.3.2.5.1)
   * [BH](https://cran.r-project.org/package=BH)
   * [RSpectra](https://cran.r-project.org/package=RSpectra)
   * [abind](https://cran.r-project.org/package=abind)

Suggests:
   
   * [foreach](https://cran.r-project.org/package=foreach)
   * [testthat](https://cran.r-project.org/package=testthat) (for unit tests)
   * [knitr](https://cran.r-project.org/package=knitr)
   * [rmarkdown](https://cran.r-project.org/package=rmarkdown)
   * [qqman](https://cran.r-project.org/package=qqman)

## PCA
   
### On a numeric matrix

   ```R
   data(hm3.chr1)
   X <- scale2(hm3.chr1$bed)
   dim(X)
   f <- flashpca(X, ndim=10, scale="none")
   ```

### On PLINK data

You can supply a path to a PLINK dataset (with extensions .bed/.bim/.fam, all
lowercase):
   ```R
   fn <- gsub("\\.bed", "",
      system.file("extdata", "data_chr1.bed", package="flashpcaR"))
   fn
   f <- flashpca(fn, ndim=10)
   ```

## UCCA (aka univariate canonical correlation analysis; aka ANOVA of each SNP on multiple phenotypes)

### On a numeric matrix

Use HapMap3 genotypes, standardise them, simulate some phenotypes, and test each
SNP for association with all phenotypes:

   ```R
   data(hm3.chr1)
   X <- scale2(hm3.chr1$bed)
   k <- 10
   B <- matrix(rnorm(ncol(X) * k), ncol=k)
   Y <- X %*% B + rnorm(nrow(X) * k)
   f1 <- ucca(X, Y, standx="none", standy="sd")
   head(f1$result)
   ```

### On PLINK data

   ```R
   fn <- gsub("\\.bed", "",
      system.file("extdata", "data_chr1.bed", package="flashpcaR"))
   fn
   f2 <- ucca(fn, Y, standx="binom2", standy="sd")
   head(f2$result)
   ```

## Sparse Canonical Correlation Analysis (SCCA)

### On a numeric matrix

Use HapMap3 genotypes, standardise them, simulate some phenotypes, and run
sparse canonical correlation analysis over all SNPs and all phenotypes:

   ```R
   data(hm3.chr1)
   X <- scale2(hm3.chr1$bed)
   k <- 10
   B <- matrix(rnorm(ncol(X) * k), ncol=k)
   Y <- X %*% B + rnorm(nrow(X) * k)
   f1 <- scca(X, Y, standx="none", standy="sd", lambda1=1e-2, lambda2=1e-3)
   diag(cor(f1$Px, f1$Py))

   # 3-fold cross-validation
   cv1 <- cv.scca(X, Y, standx="sd", standy="sd",
      lambda1=seq(1e-3, 1e-1, length=10), lambda2=seq(1e-6, 1e-3, length=5),
      ndim=3, nfolds=3)

   # Plot the canonical correlations over the penalties, for the 1st dimension
   plot(cv1, dim=1)
   ```

### On PLINK data

   ```R
   fn <- gsub("\\.bed", "",
      system.file("extdata", "data_chr1.bed", package="flashpcaR"))
   fn
   f2 <- scca(fn, Y, standx="binom2", standy="sd", lambda1=1e-2, lambda2=1e-3)
   diag(cor(f2$Px, f2$Py))

   # Cross-validation isn't yet supported for PLINK data

   ```

# LD-pruned HapMap3 and 1000Genomes example data

See the [HapMap3](HapMap3) directory

# Changelog (stable versions only)

See [CHANGELOG.txt](CHANGELOG.txt)

