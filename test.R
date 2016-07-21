
set.seed(32312)

library(flashpcaR)
library(RSpectra)
library(plink2R)

#n <- 1000
#p <- 2000
#X <- matrix(rnorm(n * p), n, p)
#X <- scale(X, center=TRUE, scale=FALSE)
dat <- read_plink("HapMap3/data", impute="random")

X <- scale(dat$bed, center=TRUE, scale=FALSE)

system.time({
   #S <- tcrossprod(X) / (nrow(X) - 1)
   S <- tcrossprod(X)
   e <- eigen(S)
})

system.time({
   s <- svd(X)
})

nthr <- 1
tol <- 1e-3
k <- 20
system.time({
   f1 <- flashpca(X, ndim=k, stand="center", transpose=FALSE, tol=tol,
      num_threads=nthr, maxiter=100)
})
f2 <- flashpca(t(X), ndim=k, stand="center", transpose=TRUE, tol=tol,
   num_threads=nthr, maxiter=100)
f3 <- flashpca(X, ndim=k, stand="none", transpose=FALSE, tol=tol,
   num_threads=nthr, maxiter=100)
f4 <- flashpca(t(X), ndim=k, stand="none", transpose=TRUE, tol=tol,
   num_threads=nthr, maxiter=100)
f5 <- svds(X, k=10)

fit <- function(evec, eval, div=1)
{
   x1 <- X %*% crossprod(X, evec) / div
   x2 <- evec %*% diag(eval)
   mean((x1 - x2)^2)
}

fit(e$vectors[,1:10], e$values[1:10])
fit(f1$vectors[,1:10], f1$values[1:10], div=nrow(X) - 1)

#stop()
#f5 <- gklBidiag(X, runif(ncol(X)), reorth=0, maxit=100)
#
#(r <- cbind(
#   e=e$val[1:k],
#   f1=f1$val,
#   f2=f2$val,
#   f3=f3$val,
#   f4=f4$val,
#   s=s$d[1:k]^2 / (nrow(X) - 1)
#))
#
#cor(r)
#
