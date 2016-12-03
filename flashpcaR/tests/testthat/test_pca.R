context("Testing PCA")

n <- 500
p <- 1000
ndim <- 50
nextra <- 100
tol <- 1e-3

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))
datf <- system.file("extdata", "data_chr1.rds", package="flashpcaR")

scale2 <- function(x, type=c("2", "1"))
{
   type <- match.arg(type)
   mult <- ifelse(type == "1", 1, 2)

   sum2 <- nrow(x) - colSums(apply(x, 2, is.na))
   p <- colSums(x, na.rm=TRUE) / (2 * sum2)
   xsd <- sqrt(mult * p * (1 - p))
   names(p) <- names(xsd) <- colnames(x)

   s <- sweep(
      sweep(x, MARGIN=2, STATS=2 * p, FUN="-"),
	 MARGIN=2, STATS=xsd, FUN="/"
   )
   s[is.na(s)] <- 0
   attr(s, "scaled:center") <- 2 * p
   attr(s, "scaled:scale") <- xsd
   s
}

compare_scales <- function(S, ...)
{
   l <- list(...)
   for(i in 1:length(l)) {
      expect_equal(attr(S, "scaled:center"), l[[i]]$center, tolerance=tol,
	 check.names=FALSE)
      expect_equal(attr(S, "scaled:scale"), l[[i]]$scale, tolerenace=tol,
	 check.names=FALSE)
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
   #X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
   #S <- scale2(X, type="1")
   dat <- readRDS(datf)
   S <- scale2(dat$bed, type="1")

   #f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f1 <- eigen(tcrossprod(S) / ncol(S), symmetric=TRUE)
   f1$projection <- f1$vectors %*% S
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
   X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
   S <- scale2(X, type="2")

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="binom2")

   compare_scales(S, f1, f2)

   compare_eigenvecs(
      f1$x[, 1:ndim], f2$projection
   )
})


test_that("Testing PCA with stand='sd'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=TRUE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="sd")

   compare_scales(S, f2)

   compare_eigenvecs(
      f1$x[, 1:ndim], f2$projection
   )
})

test_that("Testing PCA with stand='none'", {
   X <- matrix(rnorm(n * p), n, p)

   f1 <- prcomp(X, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="none")

   compare_eigenvecs(
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

   compare_eigenvecs(
      f1$x[, 1:ndim], f2$projection
   )
})

