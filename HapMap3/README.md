HapMap3 data
============

This dataset consists of 957 founders extracted from the HapMap3 phase iii data
http://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format
hapmap3_r2_b36_fwd.consensus.qc.poly.{bed,bim,fam}
using PLINK options:
   ```
   --bfile hapmap3_r2_b36_fwd.consensus.qc.poly
   --maf 0.01
   --geno 0.01
   --mind 0.01
   --hwe 5e-6
   --filter-founders
   ```

The data was then LD-thinned with PLINK using
   ```
   --indep-pairwise 1000 10 0.02
   ```

