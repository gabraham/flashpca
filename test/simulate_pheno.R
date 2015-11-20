
library(plink2R)

set.seed(12345)

dat <- read_plink("../HapMap3/data", impute="random")

n <- nrow(dat$bed)
p <- ncol(dat$bed)
k <- 100
sp <-  sample(0:1, p * k, prob=c(0.999, 0.001), replace=TRUE)
B <- matrix(rnorm(p * k) * sp, p, k)
Y <- dat$bed %*% B + matrix(rnorm(n * k), n, k)

d <- data.frame(dat$fam[, 1:2], Y)

write.table(format(d, digits=3), file="pheno.txt",
   col.names=FALSE, row.names=FALSE, quote=FALSE)

