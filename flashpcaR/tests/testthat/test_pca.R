context("Checking PCA")

n <- 500
p <- 1000
ndim <- 50
nextra <- 100
tol <- 1e-3

compare_scales <- function(S, ...)
{
   l <- list(...)
   for(i in 1:length(l)) {
      expect_equal(attr(S, "scaled:center"), l[[i]]$center, tolerance=tol)
      expect_equal(attr(S, "scaled:scale"), l[[i]]$scale, tolerenace=tol)
   }
}

compare_projections <- function(...)
{
   l <- list(...)

   # Check PCs are correlated across the methods
   s <- sapply(2:length(l), function(i) {
      abs(diag(cor(l[[1]], l[[i]])))
   })
   expect_equal(as.numeric(s), rep(1.0, length(s)), tolerance=tol)

   # We can't directly compare the eigenvalues / principal components from two
   # methods since they're only defined up to sign, so compare the dot
   # products of each PCs, and check that they are all very similar (low
   # variance across the methods)
   r <- sapply(l, function(x) {
      diag(crossprod(x))
   })
   v <- apply(r, 1, var)
   expect_equal(as.numeric(v), numeric(length(v)), tolerance=tol)
}

test_that("Testing PCA with stand='binom'", {
   X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
   q <- colMeans(X) / 2
   S <- scale(X, center=2 * q, scale=sqrt(q * (1 - q)))

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="binom")

   # Don't check prcomp scales because it returns the wrong scale
   compare_scales(S, f2)

   compare_projections(
      f1$x[, 1:ndim], f2$projection
   )
})

test_that("Testing PCA with stand='binom2'", {
   X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
   q <- colMeans(X) / 2
   S <- scale(X, center=2 * q, scale=sqrt(2 * q * (1 - q)))

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="binom2")

   compare_scales(S, f1, f2)

   compare_projections(
      f1$x[, 1:ndim], f2$projection
   )
})


test_that("Testing PCA with stand='sd'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=TRUE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="sd")

   compare_scales(S, f2)

   compare_projections(
      f1$x[, 1:ndim], f2$projection
   )
})

test_that("Testing PCA with stand='none'", {
   X <- matrix(rnorm(n * p), n, p)

   f1 <- prcomp(X, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="none")

   compare_projections(
      f1$x[, 1:ndim], f2$projection
   )
})

test_that("Testing PCA with stand='center'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=FALSE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="center")

   attr(S, "scaled:scale") <- rep(1, p)
   compare_scales(S, f2)

   compare_projections(
      f1$x[, 1:ndim], f2$projection
   )
})

