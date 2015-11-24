context("Checking PCA")

n <- 200
p <- 1000
ndim <- 50

test_that("Testing PCA with stand='binom'", {
   X <- matrix(sample(0:2, n * p, replace=TRUE), n, p)
   p <- colMeans(X) / 2
   S <- scale(X, center=TRUE, scale=sqrt(2 * p * (1 - p)))

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, mem="low")
   f3 <- flashpca(X, ndim=ndim, mem="high")

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})


test_that("Testing PCA with stand='sd'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=TRUE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="sd", mem="low")
   f3 <- flashpca(X, ndim=ndim, stand="sd", mem="high")

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})

test_that("Testing PCA with stand='none'", {
   X <- matrix(rnorm(n * p), n, p)

   f1 <- prcomp(X, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, stand="none", mem="low")
   f3 <- flashpca(X, ndim=ndim, stand="none", mem="high")

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})

test_that("Testing PCA with stand='center'", {
   X <- matrix(rnorm(n * p), n, p)
   S <- scale(X, center=TRUE, scale=FALSE)

   f1 <- prcomp(S, center=FALSE, scale.=FALSE)
   f2 <- flashpca(X, ndim=ndim, mem="low", stand="center")
   f3 <- flashpca(X, ndim=ndim, mem="high", stand="center")

   r1 <- abs(diag(cor(f1$x[, 1:ndim], f2$projection)))
   expect_equal(r1, rep(1.0, ndim), tolerance=1e-3)

   r2 <- abs(diag(cor(f1$x[, 1:ndim], f3$projection)))
   expect_equal(r2, rep(1.0, ndim), tolerance=1e-3)
})
