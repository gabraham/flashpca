context("Checking SCCA")

## It's kind of hard to test that SCCA is working, but we can at least test
## that SCCA of X with X gives 
## - canonical correlations, i.e., diag(cor(Px, Py)), are equal to 1.
## - and that the canonical covariances ``d'' are the same as the eigenvalues
## of X^T X.
##
## We use very small penalisation so as to give very close results to
## classic eigen-decomposition.


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
   
   # Essentially power method for eigen-decomposition of XX'
   s <- scca(Xb, Xb, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom", standy="binom")
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA with stand='sd'", {
   # X0 is already unit-variance so scaling within scca is moot but still
   # useful as a sanity check
   s <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
      standx="sd", standy="sd")
   eval <- eigen(crossprod(X0))$val[1:ndim]
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA with stand='none'", {
   # X0 is already unit-variance so not scaling within scca is ok
   s <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
      standx="none", standy="none")
   eval <- eigen(crossprod(X0))$val[1:ndim]
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA with stand='center'", {
   # X0 is already zero-mean so scaling within scca is moot but still
   # useful as a sanity check
   s <- scca(X0, X0, lambda1=l1, lambda2=l2, ndim=ndim,
      standx="center", standy="center")
   eval <- eigen(crossprod(X0))$val[1:ndim]
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

test_that("Testing self-self SCCA (X with X), initialising V0", {

   # Round the continuous values into {0, 1, 2}
   Xb <- apply(X0, 2, function(x) {
      round(2 * (x - min(x)) / (max(x) - min(x)))
   })

   # "binomial" standardization
   q <- colMeans(Xb) / 2
   Xbs <- scale(Xb, center=TRUE, scale=sqrt(q * (1 - q)))
   sv <- svd(Xbs)
   eval <- sv$d[1:ndim]^2
   
   # Essentially power method for eigen-decomposition of XX'
   s <- scca(Xb, Xb, lambda1=l1, lambda2=l2, ndim=ndim,	 
      standx="binom", standy="binom", V=sv$v[, 1:ndim])
   r <- diag(cor(s$Px, s$Py))
   expect_equal(rep(1, length(r)), r, tol=test.tol)
   expect_equal(eval, s$d, tol=test.tol)
})

