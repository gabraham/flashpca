
# Datasets


## HapMap3 data

Files: `HM3_thinned_autosomal_overlap.{bed,bim,fam}`

This dataset consists of 957 founders extracted from the HapMap3 phase III data
(http://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format),
hapmap3_r2_b36_fwd.consensus.qc.poly.{bed,bim,fam},
using PLINK options:
   ```
   --bfile hapmap3_r2_b36_fwd.consensus.qc.poly
   --maf 0.01
   --geno 0.01
   --mind 0.01
   --hwe 5e-6
   --filter-founders
   --autosome
   ```

The data was then LD-thinned with PLINK using
   ```
   --indep-pairwise 1000 10 0.02
   ```

The final dataset has 14,079 autosomal SNPs.

## 1000 Genomes data, phase 1 release 3

Files: `1kg.ref.phase1_release_v3.20101123_thinned_autosomal_overlap.{bed,bim,fam}`

This data contains genotypes for 1,092 individuals from the 1000 Genomes
project, for the same 14,079 SNPs from the HapMap3 dataset above.


