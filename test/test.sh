#!/bin/bash

################################################################################
# PCA on HapMap3 data

PLINK=~/Software/plink-ng/2015-10-17/plink
HM3=../HapMap3/data
TOL=1e-4
SEED=4387
NDIM=50

## High-memory
../flashpca \
   --bfile $HM3 \
   --tol $TOL \
   --mem high \
   --seed $SEED \
   --ndim $NDIM \
   --suffix _highmem.txt \
   --v

# Low memory
../flashpca \
   --bfile $HM3 \
   --tol $TOL \
   --mem low \
   --seed $SEED \
   --ndim $NDIM \
   --suffix _lowhmem.txt \
   --v

$PLINK \
   --bfile $HM3 \
   --pca $NDIM

Rscript test_pca.R && echo "PCA OK!"

#################################################################################
# Kernel PCA on HapMap3 data
../flashpca \
   --bfile $HM3 \
   --tol $TOL \
   --kernel rbf \
   --seed $SEED \
   --ndim $NDIM \
   --suffix _rbf.txt \
   --v

################################################################################
# Sparse CCA on HapMap3 data plus simulated phenotypes

Rscript simulate_pheno.R

PHENO=pheno.txt
L1=1e-2
L2=1e-3

../flashpca \
   --bfile $HM3 \
   --pheno $PHENO \
   --tol $TOL \
   --scca \
   --seed $SEED \
   --lambda1 $L1 \
   --lambda2 $L2 \
   --ndim 10 \
   --suffix _scca.txt \
   --v

Rscript test_scca.R


