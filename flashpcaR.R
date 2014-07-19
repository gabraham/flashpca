
library(plink2R)
library(flashpcaR)

set.seed(1)

dat <- read_plink("HapMap3/data", impute="random")
X <- dat$bed

r1 <- flashpca(X, do_loadings=TRUE, verbose=TRUE, stand="binom", ndim=10)

