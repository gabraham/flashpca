
set.seed(32312)

library(flashpcaR)

n <- 1000
p <- 5000
X <- matrix(rnorm(n * p), n, p)
X <- scale(X, center=TRUE, scale=FALSE)

system.time({
   S <- tcrossprod(X) / (nrow(X) - 1)
   e <- eigen(S)
})

system.time({
   s <- svd(X)
})

nthr <- 1
tol <- 1e-9
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

(r <- cbind(
   e=e$val[1:k],
   f1=f1$val,
   f2=f2$val,
   f3=f3$val,
   f4=f4$val,
   s=s$d[1:k]^2 / (nrow(X) - 1)
))

