# flashpcaR

FlashPCA performs fast principal component analysis (PCA) of single nucleotide
polymorphism (SNP) data, similar to smartpca from EIGENSOFT
(<http://www.hsph.harvard.edu/alkes-price/software/>) and shellfish
(<https://github.com/dandavison/shellfish>). FlashPCA is based on the
[https://github.com/yixuan/spectra/](Spectra) library.

Main features:

* Fast: partial PCA (_k_ = 20 dimensions) of 500,000 individuals with 100,000 SNPs in &lt;6h using 2GB RAM
* Scalable: memory requirements are bounded, scales to at least 1M individuals
* Highly accurate results
* Natively reads PLINK bed/bim/fam files
* Easy to use; can be called entirely within R (package [flashpcaR](#flashpcaR))

## Installation

Install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("umr1283/flashpcaR")
```

## Example

### PCA

#### On a numeric matrix

```R
data(hm3.chr1)
X <- scale2(hm3.chr1$bed)
dim(X)
f <- flashpca(X, ndim = 10, scale = "none")
```

#### On PLINK data

You can supply a path to a PLINK dataset (with extensions .bed/.bim/.fam, all lowercase):

```r
(fn <- gsub("\\.bed", "", system.file("extdata", "data_chr1.bed", package = "flashpcaR")))
f <- flashpca(fn, ndim = 10)
```

### UCCA (aka univariate canonical correlation analysis; aka ANOVA of each SNP on multiple phenotypes)

#### On a numeric matrix

Use HapMap3 genotypes, standardise them, simulate some phenotypes, and test each
SNP for association with all phenotypes:

```r
data(hm3.chr1)
X <- scale2(hm3.chr1$bed)
k <- 10
B <- matrix(rnorm(ncol(X) * k), ncol = k)
Y <- X %*% B + rnorm(nrow(X) * k)
f1 <- ucca(X, Y, standx = "none", standy = "sd")
head(f1$result)
```

#### On PLINK data

```r
(fn <- gsub("\\.bed", "", system.file("extdata", "data_chr1.bed", package = "flashpcaR")))
f2 <- ucca(fn, Y, standx = "binom2", standy = "sd")
head(f2$result)
```

### Sparse Canonical Correlation Analysis (SCCA)

#### On a numeric matrix

Use HapMap3 genotypes, standardise them, simulate some phenotypes, and run
sparse canonical correlation analysis over all SNPs and all phenotypes:

```r
data(hm3.chr1)
X <- scale2(hm3.chr1$bed)
k <- 10
B <- matrix(rnorm(ncol(X) * k), ncol = k)
Y <- X %*% B + rnorm(nrow(X) * k)
f1 <- scca(X, Y, standx = "none", standy = "sd", lambda1 = 1e-2, lambda2 = 1e-3)
diag(cor(f1$Px, f1$Py))

# 3-fold cross-validation
cv1 <- cv.scca(
   X, Y,
   standx = "sd",
   standy = "sd",
   lambda1 = seq(1e-3, 1e-1, length = 10),
   lambda2 = seq(1e-6, 1e-3, length = 5),
   ndim = 3,
   nfolds = 3
)

# Plot the canonical correlations over the penalties, for the 1st dimension
plot(cv1, dim = 1)
```

#### On PLINK data

```r
fn <- gsub("\\.bed", "", system.file("extdata", "data_chr1.bed", package = "flashpcaR"))
fn
f2 <- scca(fn, Y, standx = "binom2", standy = "sd", lambda1 = 1e-2, lambda2 = 1e-3)
diag(cor(f2$Px, f2$Py))
# Cross-validation isn't yet supported for PLINK data
```

## Help

Google Groups: [https://groups.google.com/forum/#!forum/flashpca-users](https://groups.google.com/forum/#!forum/flashpca-users)

## Contact

Gad Abraham, gad.abraham@baker.edu.au

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

Copyright (C) 2014-2020 Gad Abraham. All rights reserved.

Portions of this code are based on SparSNP
(<https://github.com/gabraham/SparSNP>), Copyright (C) 2011-2012 Gad Abraham
and National ICT Australia (<http://www.nicta.com.au>).
