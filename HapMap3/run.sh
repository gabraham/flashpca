#!/bin/bash

export LD_LIBRARY_PATH=/software/plink-ng/zlib-1.2.8
export PATH=~/Software/shellfish:~/Software/plink-ng:$PATH

DAT=data

~/Code/flashpca/flashpca --bfile ${DAT} \
   --ndim 100 \
   --numthreads 8 2>&1 | tee log

# smartpca won't accept missing phenotypes so make some up
awk '{print $1, $2, $3, $4, $5, "2"}' data.fam > tmp

mv tmp data.fam

cat > ${DAT}.par<<EOF
genotypename: $DAT.bed
snpname: $DAT.bim
indivname: $DAT.fam
evecoutname: ${DAT}.pca.evec
evaloutname: ${DAT}.eval
altnormstyle: NO
numoutevec: 100
numoutlieriter: 0
qtmode: 0
EOF

time /software/EIG4.2/src/smartpca -p ${DAT}.par > eigen.log

# original shellfish doesn't run with current gtool
# use this slightly modified version
time /software/shellfish/shellfish_gad.py \
   --file ${DAT} --pca --numpcs 100 -v > shellfish.log

Rscript plot.R

