#!/bin/bash

################################################################################
# This script will 
# 1) split your PLINK data into training and testing (with proportions
#    $PROP and 1 - $PROP)
# 2) Run SCCA on the training data
# 3) Project the test SNPs data onto the eigenvectors, producing "principal
#    components" in the test data
#
# You must specify ROOT and PHENO, where ROOT is the PLINK data rootname
# (without the bed/bim,fam), and PHENO is a text phenotype file in the format
# 
# FID, IID, pheno1, pheno2, ...
# 
# The phenotype file *must be in the same order* as the PLINK FAM file.
# 
################################################################################

ROOT=
PHENO=
PROP=0.8
pen1="1.1e-3 1.3e-3 1.5e-3 1.7e-3 1.9e-3 2.1e-3 \
   2.3e-3 2.5e-3 2.7e-3 2.9e-3 3.1e-3 \
   3.3e-3 3.5e-3 3.7e-3 3.9e-3 \
   4.1e-3 4.3e-3 4.5e-3 4.7e-3 4.9e-3 \
   4.4e-3 4.6e-3 4.8e-3 5.1e-3 5.2e-3"
pen2="1e-5 1e-4 1e-3"
NUMPROC=10

################################################################################

BN=$(basename $ROOT)
PHENOBN=$(basename $PHENO)

# manual cross-validation split
PHENOTRAIN=${PHENOBN/.txt/_train}.txt

awk -vp=$PROP 'BEGIN{srand(1)} {if(rand() < p) print}' $PHENO > $PHENOTRAIN
awk '{print $1, $2}' $PHENOTRAIN > samples_train.txt

plink --bfile $ROOT \
   --keep samples_train.txt \
   --make-bed \
   --out ${BN}_train

plink --bfile $ROOT \
   --remove samples_train.txt \
   --make-bed \
   --out ${BN}_test

pen1="1.1e-3 1.3e-3 1.5e-3 1.7e-3 1.9e-3 2.1e-3 \
   2.3e-3 2.5e-3 2.7e-3 2.9e-3 3.1e-3 \
   3.3e-3 3.5e-3 3.7e-3 3.9e-3 \
   4.1e-3 4.3e-3 4.5e-3 4.7e-3 4.9e-3 \
   4.4e-3 4.6e-3 4.8e-3 5.1e-3 5.2e-3"
pen2="1e-5 1e-4 1e-3"


echo $pen1 > lambda1.txt
echo $pen2 > lambda2.txt

parallel -P$NUMPROC "flashpca \
      --scca \
      --mem low \
      --ndim 10 \
      --lambda1 {1} --lambda2 {2} \
      --bfile ${BN}_train \
      --pheno $PHENOTRAIN \
      --stand sd \
      --outpcx pcsX_{1}_{2}.txt \
      --outpcy pcsY_{1}_{2}.txt \
      --outval eigenvalues_{1}_{2}.txt \
      --outvecx eigenvectorsX_{1}_{2}.txt \
      --outvecy eigenvectorsY_{1}_{2}.txt \
      --numthreads 10 2>& 1| tee log_{1}_{2}" \
      ::: $pen1 ::: $pen2

parallel -P$NUMPROC "predict \
      --stand sd \
      --bfile ${BN}_test \
      --w eigenvectorsX_{1}_{2}.txt \
      --numthreads 8 \
      --out predX_{1}_{2}.txt" \
      ::: $pen1 ::: $pen2
   



