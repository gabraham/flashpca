#!/bin/bash

# Requires the following files from the PCCA pipeline
# (https://github.com/pachterlab/PCCA):
#
# geuvadis.bed
# geuvadis.bim
# geuvadis.fam
# expression_plink.tsv
# fin_ids.txt
# gbr_ids.txt
# tsi_ids.txt
# yri_ids.txt
#

PLINK=~/Software/plink1.9/2019-11-30/plink

for pop in fin gbr tsi yri
do
   eval $PLINK --bfile geuvadis \
      --keep ${pop}_ids.txt \
      --maf 0.1 \
      --mind 0.1 \
      --geno 0.1 \
      --hwe 1e-6 \
      --make-bed \
      --out geuvadis_${pop}
done

# Get the SNPs that remain in all populations
awk '{print $2}' geuvadis_*.bim  | sort | uniq -c | awk '$1 == 4 {print $2}' \
   > snps.txt

eval $PLINK --bfile geuvadis \
   --extract snps.txt \
   --indep-pairwise 1000 50 0.1

eval $PLINK --bfile geuvadis \
   --extract plink.prune.in \
   --make-bed \
   --out geuvadis_thinned


