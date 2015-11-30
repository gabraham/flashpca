context("Checking PCA")

n <- 200
p <- 1000
ndim <- 50
nextra <- 50

test_that("Testing PCA with stand='binom'", {
   X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
   q <- colMeans(X) / 2
   S <- scale(X, center=TRUE, scale=sqrt(q * (1 - q)))

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, mem="low", nextra=nextra)
   f3 <- flashpca(X, ndim=ndim, mem="high", nextra=nextra)

   expect_equal(attr(S, "scaled:center"), f2$center)
   expect_equal(attr(S, "scaled:scale"), f2$scale)
   expect_equal(attr(S, "scaled:center"), f3$center)
   expect_equal(attr(S, "scaled:scale"), f3$scale)

   r0 <- abs(diag(cor(f2$projection, f3$projection)))
   expect_equal(r0, rep(1.0, ndim), tolerance=1e-3)

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})


test_that("Testing PCA with stand='sd'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=TRUE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="sd", mem="low", nextra=nextra)
   f3 <- flashpca(X, ndim=ndim, stand="sd", mem="high", nextra=nextra)

   expect_equal(attr(S, "scaled:center"), f2$center)
   expect_equal(attr(S, "scaled:scale"), f2$scale)
   expect_equal(attr(S, "scaled:center"), f3$center)
   expect_equal(attr(S, "scaled:scale"), f3$scale)

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})

test_that("Testing PCA with stand='none'", {
   X <- matrix(rnorm(n * p), n, p)

   f1 <- prcomp(X, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="none", mem="low", nextra=nextra)
   f3 <- flashpca(X, ndim=ndim, stand="none", mem="high", nextra=nextra)

   expect_identical(numeric(0), f2$center)
   expect_identical(numeric(0), f2$scale)
   expect_identical(numeric(0), f3$center)
   expect_identical(numeric(0), f3$scale)

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})

test_that("Testing PCA with stand='center'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=FALSE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, mem="low", stand="center", nextra=nextra)
   f3 <- flashpca(X, ndim=ndim, mem="high", stand="center", nextra=nextra)

   expect_equal(attr(S, "scaled:center"), f2$center)
   expect_equal(rep(1, p), f2$scale)
   expect_equal(attr(S, "scaled:center"), f3$center)
   expect_equal(rep(1, p), f3$scale)

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})

