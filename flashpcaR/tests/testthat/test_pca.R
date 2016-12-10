context("Testing PCA")

n <- 500
p <- 1000
ndim <- 50
nextra <- 100
tol <- 1e-4

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

compare_scales <- function(S, ...)
{
   l <- list(...)
   for(i in 1:length(l)) {
      expect_equal(attr(S, "scaled:center"),
	 l[[i]]$center, tolerance=tol, check.names=FALSE)
      expect_equal(attr(S, "scaled:scale"),
	 l[[i]]$scale, tolerenace=tol, check.names=FALSE)
   }
}

compare_eigenvecs <- function(...)
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
   S <- scale2(hm3.chr1$bed, type="1")

   #f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f1 <- eigen(tcrossprod(S) / ncol(S), symmetric=TRUE)
   f1$projection <- with(f1, vectors %*% diag(sqrt(values)))
   f2 <- flashpca(S, ndim=ndim, stand="none")
   f3 <- flashpca(bedf, ndim=ndim, stand="binom")

   # Don't check prcomp scales because it returns the wrong scale (S is
   # already standardised). Also don't compare f2 since we didn't standardise
   # in that call.
   compare_scales(S, f3)

   compare_eigenvecs(
      f1$vectors[, 1:ndim], f2$vectors, f3$vectors
   )

   compare_eigenvecs(
      f1$projection[, 1:ndim], f2$projection, f3$projection
   )
})

test_that("Testing PCA with stand='binom2'", {
   S <- scale2(hm3.chr1$bed, type="2")

   f1 <- eigen(tcrossprod(S) / ncol(S), symmetric=TRUE)
   f1$projection <- with(f1,
      vectors[, 1:ndim] %*% diag(sqrt(values[1:ndim])))
   f2 <- flashpca(S, ndim=ndim, stand="none")
   f3 <- flashpca(bedf, ndim=ndim, stand="binom2")

   compare_scales(S, f3)

   compare_eigenvecs(
      f1$projection[, 1:ndim], f2$projection, f3$projection
   )
})


test_that("Testing PCA with stand='sd'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=TRUE)

   f1 <- eigen(tcrossprod(S) / ncol(S), symmetric=TRUE)
   f1$projection <- with(f1,
      vectors[, 1:ndim] %*% diag(sqrt(values[1:ndim])))
   f2 <- flashpca(X, ndim=ndim, stand="sd")

   compare_scales(S, f2)

   compare_eigenvecs(
      f1$projection[, 1:ndim], f2$projection
   )
})

test_that("Testing PCA with stand='none'", {
   X <- matrix(rnorm(n * p), n, p)

   f1 <- eigen(tcrossprod(X) / ncol(X), symmetric=TRUE)
   f1$projection <- with(f1, vectors %*% diag(sqrt(values)))
   f2 <- flashpca(X, ndim=ndim, stand="none")

   compare_eigenvecs(
      f1$projection[, 1:ndim], f2$projection
   )
})

test_that("Testing PCA with stand='center'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=FALSE)

   f1 <- eigen(tcrossprod(S) / ncol(S), symmetric=TRUE)
   f1$projection <- with(f1,
      vectors[, 1:ndim] %*% diag(sqrt(values[1:ndim])))
   f2 <- flashpca(X, ndim=ndim, stand="center")

   attr(S, "scaled:scale") <- rep(1, p)
   compare_scales(S, f2)

   compare_eigenvecs(
      f1$projection[, 1:ndim], f2$projection
   )
})

