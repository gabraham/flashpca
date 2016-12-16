context("Testing PCA checking")

n <- 500
p <- 1000
ndim <- 50
nextra <- 100
tol <- 1e-3

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

test_that("Testing PCA with stand='binom'", {
   S <- scale2(hm3.chr1$bed, type="1")

   f2 <- flashpca(S, ndim=ndim, stand="none")
   f3 <- flashpca(bedf, ndim=ndim, stand="binom")

   XXU <- S %*% crossprod(S, f2$vectors) / ncol(S)
   UD2 <- f2$vectors %*% diag(f2$values)
   err <- colSums((XXU - UD2)^2)
   mse <- sum(err) / (nrow(S) * ndim)

   c1 <- check(S, stand="none", evec=f2$vectors, eval=f2$values)
   c2 <- check(bedf, stand="binom", evec=f2$vectors, eval=f2$values)

   expect_equal(err, c1$err)
   expect_equal(err, c2$err)
   expect_equal(numeric(ndim), c1$err)
   expect_equal(mse, c1$mse)
   expect_equal(mse, c2$mse)
})

test_that("Testing input checking", {

   X <- scale2(hm3.chr1$bed, type="1")

   ndim <- 5
   evec <- matrix(rnorm((nrow(X) + 3) * ndim), nrow(X) + 3, ndim)
   eval<- rnorm(ndim)^2

   # Test incompatible number of rows in evec and X
   expect_error(
      check(X, stand="none", evec=evec, eval=eval)
   )
   expect_error(
      check(bedf, stand="none", evec=evec, eval=eval)
   )

   # Test incompatible number of dimensions between eval and evec
   eval <- eval[1:3]
   expect_error(
      check(X, stand="none", evec=evec, eval=eval)
   )
   expect_error(
      check(bedf, stand="none", evec=evec, eval=eval)
   )
})

