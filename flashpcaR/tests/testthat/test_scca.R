context("Checking SCCA")

## It's kind of hard to test that SCCA is working, but we can at least test
## that SCCA of X with X gives 
## - canonical correlations, i.e., diag(cor(Px, Py)), are equal to 1.
## - and that the canonical covariances ``d'' are the same as the eigenvalues
## of X^T X.

n <- 500
p <- 100
k <- 50
m <- 10
M <- matrix(rnorm(n * m), n, m)
Bx <- matrix(rnorm(m * p), m, p)
By <- matrix(rnorm(m * k), m, k)
X0 <- scale(M %*% Bx + rnorm(n * p))
Y0 <- scale(M %*% By + rnorm(n * k))

l1 <- 1e-6
l2 <- 1e-6
test.tol <- 1e-4
ndim <- min(n, p, k, 5)

test_that("Testing self-self SCCA (X with X)", {

   # Round the continuous values into {0, 1, 2}
   Xb <- apply(X0, 2, function(x) {
      round(2 * (x - min(x)) / (max(x) - min(x)))
   })

   # "binomial" standardization
   q <- colMeans(Xb) / 2
   Xbs <- scale(Xb, center=TRUE, scale=sqrt(q * (1 - q)))
   eval <- eigen(crossprod(Xbs))$val[1:ndim]
   
   s <- scca(Xb, Xb, lambda1=l1, lambda2=l2, ndim=ndim, stand="binom")
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA with stand='sd'", {
   s <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim, stand="sd")
   eval <- eigen(crossprod(X0))$val[1:ndim]
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA with stand='none'", {
   s <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim, stand="none")
   eval <- eigen(crossprod(X0))$val[1:ndim]
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA with stand='center'", {
   s <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim, stand="center")
   eval <- eigen(crossprod(X0))$val[1:ndim]
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

