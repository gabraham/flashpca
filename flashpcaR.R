
library(plink2R)
library(flashpcaR)

set.seed(1)

dat <- read_plink("HapMap3/data", impute="random")
X <- dat$bed

################################################################################
# Principal component analysis
system.time({
   r1 <- flashpca(X, do_loadings=TRUE, verbose=FALSE, stand="binom", ndim=10)
})


################################################################################
# Sparse canonical correlation analysis with simulated Y 
p <- ncol(X)
k <- 50
B <- matrix(rnorm(p * k), p, k)
Y <- X %*% B

system.time({
   r2 <- scca(X, Y, ndim=5, lambda1=0.01, lambda2=0.01)
})

diag(cor(r2$Px, r2$Py))

m <- read.table("HapMap3/relationships_w_pops_121708.txt", header=TRUE, sep="",
   stringsAsFactors=FALSE)
rownames(m) <- m[,2]
pop <- factor(m[dat$fam[,2], 7])


par(pty="s", mar=c(1, 1, 1, 1), mfrow=c(5, 5))
for(i in 1:5) {
   for(j in 1:5) {
      if(i == j) {
	 plot(NULL, xlim=c(0, 1), ylim=c(0, 1))
      } else {
	 plot(r2$Px[,i], r2$Py[,j], col=pop)
      }
   }
}

