rm(list=ls())

library(RSpectra)
library(plink2R)
library(abind)

dat <- read_plink("data", impute="none")

scale2 <- function(X)
{
   p <- colSums(X, na.rm=TRUE) / (2 * colSums(!is.na(X)))
   S <- sweep(
      sweep(X, MARGIN=2, STATS=2 * p, FUN="-"),
      MARGIN=2, STATS=sqrt(2 * p * (1 - p)), FUN="/"
   )
   S[is.na(S)] <- 0
   attr(S, "scaled:center") <- 2 * p
   attr(S, "scaled:scale") <- sqrt(2 * p * (1 - p))
   S
}

X <- scale2(dat$bed)

X <- X / sqrt(ncol(X))

k <- 10
tol <- 1e-6
s1 <- svd(X, nu=k, nv=k)
s2 <- svds(X, k=k, tol=tol)

# Spectra / batch
system(paste0("../flashpca --bfile data --ndim ", k, " --suffix .batch.txt ",
   "--tol ", tol, " --outload loadings.batch.txt",
   " --outmeansd meansd.batch.txt --precision 10"))

# Spectra / online
system(paste0("../flashpca --bfile data --online --ndim ", k,
   " --suffix .online.txt --tol ", tol, " --outload loadings.online.txt",
   " --outmeansd meansd.online.txt --precision 10"))

# Randomised / batch
system(paste0("../flashpca --rand --bfile data --ndim ", k, " --suffix .rand.batch.txt ",
   "--tol ", tol, " --outload loadings.rand.batch.txt",
   " --outmeansd meansd.rand.batch.txt --precision 10"))

# Randomised / online
system(paste0("../flashpca --rand --bfile data --online --ndim ", k,
   " --suffix .rand.online.txt --tol ", tol, " --outload loadings.rand.online.txt",
   " --outmeansd meansd.rand.online.txt --precision 10"))

# Projection
system(paste0("../flashpca --bfile data --project --inmeansd meansd.online.txt",
   " --outproj projections.online.txt --inload loadings.online.txt -v",
   " --precision 10"))

evec3 <- as.matrix(read.table("eigenvectors.batch.txt", header=FALSE))
evec4 <- as.matrix(read.table("eigenvectors.online.txt", header=FALSE))
evec5 <- as.matrix(read.table("eigenvectors.rand.batch.txt", header=FALSE))
evec6 <- as.matrix(read.table("eigenvectors.rand.online.txt", header=FALSE))

eval3 <- scan("eigenvalues.batch.txt")
eval4 <- scan("eigenvalues.online.txt")
eval5 <- scan("eigenvalues.rand.batch.txt")
eval6 <- scan("eigenvalues.rand.online.txt")

load3 <- as.matrix(read.table("loadings.batch.txt", header=TRUE, row.names=1))
load4 <- as.matrix(read.table("loadings.online.txt", header=TRUE, row.names=1))
load5 <- as.matrix(read.table("loadings.rand.batch.txt", header=TRUE, row.names=1))
load6 <- as.matrix(read.table("loadings.rand.online.txt", header=TRUE, row.names=1))

pc3 <- as.matrix(read.table("pcs.batch.txt", header=TRUE))
pc4 <- as.matrix(read.table("pcs.online.txt", header=TRUE))
pc5 <- as.matrix(read.table("pcs.rand.batch.txt", header=TRUE))
pc6 <- as.matrix(read.table("pcs.rand.online.txt", header=TRUE))

pve3 <- scan("pve.batch.txt")
pve4 <- scan("pve.online.txt")
pve5 <- scan("pve.rand.batch.txt")
pve6 <- scan("pve.rand.online.txt")

msd3 <- read.table("meansd.batch.txt", header=TRUE, row.names=1)
msd4 <- read.table("meansd.online.txt", header=TRUE, row.names=1)
msd5 <- read.table("meansd.rand.batch.txt", header=TRUE, row.names=1)
msd6 <- read.table("meansd.rand.online.txt", header=TRUE, row.names=1)

proj4 <- read.table("projections.online.txt", header=TRUE)

err.tol <- 1e-9

################################################################################
# Scaling
scl.mean <- cbind(attr(X, "scaled:center"),
   msd3[,1], msd4[,1], msd5[,1], msd6[,1])
scl.mean.rmse <- sapply(1:ncol(scl.mean), function(i) {
   sapply(1:ncol(scl.mean), function(j) {
      sqrt(mean((scl.mean[,i] - scl.mean[,j])^2))
   })
})
stopifnot(all(scl.mean.rmse < err.tol))

scl.sd <- cbind(attr(X, "scaled:scale"),
   msd3[,2], msd4[,2], msd5[,2], msd6[,2])
scl.sd.rmse <- sapply(1:ncol(scl.sd), function(i) {
   sapply(1:ncol(scl.sd), function(j) {
      sqrt(sd((scl.sd[,i] - scl.sd[,j])^2))
   })
})
stopifnot(all(scl.sd.rmse < err.tol))

################################################################################
# Eigenvalues
eval <- cbind(svd=s1$d[1:k]^2, RSpectra=s2$d^2,
   FlashPCA_batch=eval3, FlashPCA_online=eval4,
   FlashPCA_rand_batch=eval5, FlashPCA_rand_online=eval6)
eval.rmse <- sapply(1:ncol(eval), function(i) {
   sapply(1:ncol(eval), function(j) {
      sqrt(mean((eval[,i] - eval[,j])^2))
   }) 
})
stopifnot(all(eval.rmse < err.tol))

################################################################################
# Eigenvectors
evec <- abind(list(s1$u, s2$u, evec3, evec4, evec5, evec6), along=3)
evec.rmse <- sapply(1:dim(evec)[3], function(i) {
   sapply(1:dim(evec)[3], function(j) {
      r <- sapply(1:k, function(m) {
	 min(
      	    mean((evec[,m,i] - evec[,m,j]))^2,
      	    mean((evec[,m,i] + evec[,m,j]))^2
      	 )
      })
      sqrt(sum(r))
   })
})
stopifnot(all(evec.rmse < err.tol))

################################################################################
# Principal components
pc <- abind(list(X %*% s1$v, X %*% s2$v, pc3, pc4, pc5, pc6), along=3)
pc.rmse <- sapply(1:dim(pc)[3], function(i) {
   sapply(1:dim(pc)[3], function(j) {
      r <- sapply(1:k, function(m) {
	 min(
	    mean((pc[,m,i] - pc[,m,j]))^2,
	    mean((pc[,m,i] + pc[,m,j]))^2
	 )
      })
      sqrt(sum(r))
   })
})
stopifnot(all(pc.rmse < err.tol))

################################################################################
# Proportion of variance explained
tr <- sum(X^2)
pve <- cbind(svd=s1$d[1:k]^2 / tr, RSpectra=s2$d^2 / tr,
   FlashPCA_batch=pve3, FlashPCA_online=pve4,
   FlashPCA_rand=pve5, FlashPCA_rand_online=pve6)
pve.rmse <- sapply(1:ncol(pve), function(i) {
   sapply(1:ncol(pve), function(j) {
      sqrt(mean((pve[,i] - pve[,j])^2))
   }) 
})
stopifnot(all(pve.rmse < err.tol))

################################################################################
# SNP loadings
load <-  abind(list(s1$v, s2$v, load3, load4, load5, load6), along=3)
load.rmse <- sapply(1:dim(load)[3], function(i) {
   sapply(1:dim(load)[3], function(j) {
      r <- sapply(1:k, function(m) {
	 min(
	    mean((load[,m,i] - load[,m,j]))^2,
	    mean((load[,m,i] + load[,m,j]))^2
	 )
      })
      sqrt(sum(r))
   })
})
stopifnot(all(load.rmse < err.tol))

################################################################################
# Projection (=PCs if same samples used for deriving eigenvectors)
proj <- abind(list(X %*% s1$v, X %*% s2$v, pc4, proj4), along=3)
proj.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      r <- sapply(1:k, function(m) {
	 min(
	    mean((proj[,m,i] - proj[,m,j]))^2,
	    mean((proj[,m,i] + proj[,m,j]))^2
	 )
      })
      sqrt(sum(r))
   })
})
stopifnot(all(proj.rmse < err.tol))

