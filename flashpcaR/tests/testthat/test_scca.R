context("Checking SCCA")

n <- 500
p <- 100
k <- 50
m <- 5
M <- matrix(rnorm(n * m), n, m)
Bx <- matrix(rnorm(m * p), m, p)
By <- matrix(rnorm(m * k), m, k)
X0 <- scale(M %*% Bx + rnorm(n * p))
Y0 <- scale(M %*% By + rnorm(n * k))

test_that("Testing SCCA", {

   X <- apply(X0, 2, function(x) {
      round(2 * (x - min(x)) / (max(x) - min(x)))
   })
   Y <- apply(Y0, 2, function(x) {
      round(2 * (x - min(x)) / (max(x) - min(x)))
   })
   
   s1 <- scca(X, Y, lambda1=1e-2, lambda2=1e-2, ndim=5, stand="binom")
})


test_that("Testing SCCA with stand='sd'", {
   s <- scca(X0, Y0, lambda1=1e-2, lambda2=1e-2, ndim=5, stand="sd")
})

test_that("Testing SCCA with stand='none'", {
   s <- scca(X0, Y0, lambda1=1e-2, lambda2=1e-2, ndim=5, stand="none")
})

test_that("Testing SCCA with stand='center'", {
   s <- scca(X0, Y0, lambda1=1e-2, lambda2=1e-2, ndim=5, stand="center")
})

