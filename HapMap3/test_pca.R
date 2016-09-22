rm(list=ls())

library(RSpectra)
library(plink2R)
library(abind)

dat1 <- read_plink(
   "HM3_thinned_autosomal_overlap", impute="none")
dat2 <- read_plink(
   "1kg.ref.phase1_release_v3.20101123_thinned_autosomal_overlap",
   impute="none")

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

X <- scale2(dat1$bed)

X <- X / sqrt(ncol(X))

k <- 10
tol <- 1e-6
s1 <- svd(X, nu=k, nv=k)
s2 <- svds(X, k=k, tol=tol)

maf <- colMeans(dat1$bed, na.rm=TRUE) / 2
write.table(data.frame(SNP=colnames(X), MAF=round(maf, 20)),
   file="maf.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Spectra
system(paste0("../flashpca ",
   " --bfile HM3_thinned_autosomal_overlap --ndim ", k,
   " --tol ", tol, " --outload loadings.txt",
   " --outmeansd meansd.txt --precision 20"))

# Projection onto same data
system(paste0(
   "../flashpca",
   " --bfile HM3_thinned_autosomal_overlap",
   " --project --inmeansd meansd.txt",
   " --outproj projections.txt --inload loadings.txt -v",
   " --precision 20"))

# Projection onto other data
system(paste0(
   "../flashpca ",
   " --bfile 1kg.ref.phase1_release_v3.20101123_thinned_autosomal_overlap",
   " --project --inmeansd meansd.txt",
   " --outproj projections.1kg.txt --inload loadings.txt -v",
   " --precision 20"))

# Projection using MAF instead of means+sd
system(paste0(
   "../flashpca",
   " --bfile HM3_thinned_autosomal_overlap",
   " --project --inmaf maf.txt",
   " --outproj projections.maf.txt --inload loadings.txt -v",
   " --precision 20"))

# PCA checking mode
d.chk <- read.table(pipe(paste0(
   "../flashpca ",
   " --bfile HM3_thinned_autosomal_overlap",
   " --check",
   " --outval eigenvalues.txt --outvec eigenvectors.txt -v",
   " --precision 20 --notime | awk -F',| ' '/eval/{print $2, $7}'")),
   header=FALSE, sep="", stringsAsFactors=FALSE)

evec.d <- read.table("eigenvectors.txt", header=TRUE)
evec <- as.matrix(evec.d[, -(1:2)])
eval <- scan("eigenvalues.txt")
loadings.d <- read.table("loadings.txt", header=TRUE)
loadings <- as.matrix(loadings.d[, -(1:2)])
pcs.d <- read.table("pcs.txt", header=TRUE)
pcs <- as.matrix(pcs.d[, -(1:2)])
pve <- scan("pve.txt")
msd.d <- read.table("meansd.txt", header=TRUE)
msd <- as.matrix(msd.d[, -(1:2)])
proj.d <- read.table("projections.txt", header=TRUE)
proj <- as.matrix(proj.d[, -(1:2)])
proj.1kg.d <- read.table("projections.1kg.txt", header=TRUE)
proj.1kg <- as.matrix(proj.1kg.d[, -(1:2)])

stopifnot(all(evec.d[,1] == dat1$fam[,1]))
stopifnot(all(evec.d[,2] == dat1$fam[,2]))

stopifnot(all(pcs.d[,1] == dat1$fam[,1]))
stopifnot(all(pcs.d[,2] == dat1$fam[,2]))

stopifnot(all(proj.d[,1] == dat1$fam[,1]))
stopifnot(all(proj.d[,2] == dat1$fam[,2]))

stopifnot(all(msd.d[,1] == dat1$bim[,2]))
stopifnot(all(msd.d[,2] == dat1$bim[,5]))

stopifnot(all(loadings.d[,1] == dat1$bim[,2]))
stopifnot(all(loadings.d[,2] == dat1$bim[,5]))

# X is already divided by sqrt(p)
XXU <- X %*% crossprod(X, evec)
UD <- evec %*% diag(eval)
sse <- colSums((XXU - UD)^2)
rmse <- sqrt(sum(sse) / (nrow(XXU) * ncol(XXU)))
colnames(d.chk) <- c("eigenvalue", "sse_observed")
d.chk$sse_expected <- sse

# Project 1KG data onto HM3 PCA
proj.s1.1kg <- scale(dat2$bed,
   center=attr(X, "scaled:center"),
   scale=attr(X, "scaled:scale")) %*% s1$v / sqrt(ncol(dat2$bed))

err.tol <- 1e-6

################################################################################
# Scaling
scl.mean <- cbind(attr(X, "scaled:center"), msd[,1])
scl.mean.rmse <- sapply(1:ncol(scl.mean), function(i) {
   sapply(1:ncol(scl.mean), function(j) {
      sqrt(mean((scl.mean[,i] - scl.mean[,j])^2))
   })
})
stopifnot(all(scl.mean.rmse < err.tol))

scl.sd <- cbind(attr(X, "scaled:scale"), msd[,2])
scl.sd.rmse <- sapply(1:ncol(scl.sd), function(i) {
   sapply(1:ncol(scl.sd), function(j) {
      sqrt(sd((scl.sd[,i] - scl.sd[,j])^2))
   })
})
stopifnot(all(scl.sd.rmse < err.tol))

################################################################################
# Eigenvalues
eval <- cbind(svd=s1$d[1:k]^2, RSpectra=s2$d^2, FlashPCA=eval)
eval.rmse <- sapply(1:ncol(eval), function(i) {
   sapply(1:ncol(eval), function(j) {
      sqrt(mean((eval[,i] - eval[,j])^2))
   }) 
})
stopifnot(all(eval.rmse < err.tol))

################################################################################
# Eigenvectors
evec <- abind(list(s1$u, s2$u, evec), along=3)
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
pc <- abind(list(X %*% s1$v, X %*% s2$v, pcs), along=3)
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
pve <- cbind(svd=s1$d[1:k]^2 / tr, RSpectra=s2$d^2 / tr, FlashPCA=pve)
pve.rmse <- sapply(1:ncol(pve), function(i) {
   sapply(1:ncol(pve), function(j) {
      sqrt(mean((pve[,i] - pve[,j])^2))
   }) 
})
stopifnot(all(pve.rmse < err.tol))

################################################################################
# SNP loadings
load <- abind(list(s1$v, s2$v, loadings), along=3)
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
proj.r <- abind(list(X %*% s1$v, X %*% s2$v, pcs, proj), along=3)
proj.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      r <- sapply(1:k, function(m) {
	 min(
	    mean((proj.r[,m,i] - proj.r[,m,j]))^2,
	    mean((proj.r[,m,i] + proj.r[,m,j]))^2
	 )
      })
      sqrt(sum(r))
   })
})
stopifnot(all(proj.rmse < err.tol))

################################################################################
# Projection of 1KG data on HM3 data
proj.1kg.r <- abind(list(proj.s1.1kg, proj.1kg), along=3)
proj.1kg.rmse <- sapply(1:dim(proj.1kg.r)[3], function(i) {
   sapply(1:dim(proj.1kg.r)[3], function(j) {
      r <- sapply(1:k, function(m) {
	 min(
	    mean((proj.1kg.r[,m,i] - proj.1kg.r[,m,j]))^2,
	    mean((proj.1kg.r[,m,i] + proj.1kg.r[,m,j]))^2
	 )
      })
      sqrt(sum(r))
   })
})
stopifnot(all(proj.1kg.rmse < err.tol))

################################################################################
# Check the PCA-checking, i.e., the sum squared errors for
# (X X' U / div - U D)^2 in each dimension.
check.mse <- with(d.chk, (sse_observed - sse_expected)^2)
stopifnot(all(check.mse < err.tol))


