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
   S
}

X <- scale2(dat$bed)

X <- X / sqrt(ncol(X))

k <- 50
tol <- 1e-6
s1 <- svd(X, nu=k, nv=k)
s2 <- svds(X, k=k, tol=tol)
system(paste0("../flashpca --bfile data --ndim ", k, " --suffix .batch.txt ",
   "--tol ", tol, " --outload loadings.batch.txt"))
system(paste0("../flashpca --bfile data --online --ndim ", k,
   " --suffix .online.txt --tol ", tol, " --outload loadings.online.txt"))
system(paste0("../flashpca --bfile data --suffix .batch.txt ",
   "--tol ", tol, " --outload loadings.batch.txt"))
system(paste0("../flashpca --bfile data --online --ndim ", k,
   " --suffix .online.txt --tol ", tol, " --outload loadings.online.txt"))

evec3 <- as.matrix(read.table("eigenvectors.batch.txt", header=FALSE))
evec4 <- as.matrix(read.table("eigenvectors.online.txt", header=FALSE))
eval3 <- scan("eigenvalues.batch.txt")
eval4 <- scan("eigenvalues.online.txt")
load3 <- as.matrix(read.table("loadings.batch.txt", header=FALSE))
load4 <- as.matrix(read.table("loadings.online.txt", header=FALSE))
pc3 <- as.matrix(read.table("pcs.batch.txt", header=TRUE))
pc4 <- as.matrix(read.table("pcs.online.txt", header=TRUE))
pve3 <- scan("pve.batch.txt")
pve4 <- scan("pve.online.txt")

################################################################################
# Eigenvalues
eval <- cbind(svd=s1$d[1:k]^2, RSpectra=s2$d^2,
   FlashPCA_batch=eval3, FlashPCA_online=eval4)
eval.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      sqrt(mean((eval[,i] - eval[,j])^2))
   }) 
})

################################################################################
# Eigenvectors
evec <- abind(list(s1$u, s2$u, evec3, evec4), along=3)
evec.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      sqrt(mean((evec[,,i] - evec[,,j]))^2)
   })
})

################################################################################
# Principal components
pc <- abind(list(X %*% s1$v, X %*% s2$v, pc3, pc4), along=3)
pc.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      sqrt(mean((pc[,,i] - pc[,,j]))^2)
   })
})

################################################################################
# Proportion of variance explained
tr <- sum(X^2)
pve <- cbind(svd=s1$d[1:k]^2 / tr, RSpectra=s2$d^2 / tr,
   FlashPCA_batch=pve3, FlashPCA_online=pve4)
pve.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      sqrt(mean((pve[,i] - pve[,j])^2))
   }) 
})

################################################################################
# SNP loadings
load <-  abind(list(s1$v, s2$v, load3, load4), along=3)
load.rmse <- sapply(1:4, function(i) {
   sapply(1:4, function(j) {
      sqrt(mean((load[,,i] - load[,,j]))^2)
   })
})

################################################################################
# Projection (=PCs if same samples used for deriving eigenvectors)


