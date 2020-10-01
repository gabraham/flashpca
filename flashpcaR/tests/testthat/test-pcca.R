
context("Testing PCCA")

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

k <- 50
m <- ncol(hm3.chr1$bed)
X <- scale2(hm3.chr1$bed)
n <- nrow(X)
By <- matrix(rnorm(m * k), m, k)
Y <- scale(X %*% By + rnorm(n * k))
colnames(Y) <- paste0("marker", 1:k)

test_that("Testing pcca", {
   kx <- 3
   ky <- 5
   ndim <- min(kx, ky)
   nfolds <- 3
   s1 <- pcca(X, Y, kx=kx, ky=ky)

   expect_equal(s1$ndim, ndim)
   expect_true(all(s1$r >= 0))
   expect_equivalent(diag(cor(s1$Px, s1$Py)), s1$r)

   s2 <- pcca(X, X, kx=kx, ky=ky)
   expect_equivalent(s2$r, rep(1, ndim))

   expect_error(pcca(X, NULL))
   expect_error(pcca(NULL, Y))
   expect_error(pcca(X[-1,], Y))
})

test_that("Testing cv.pcca", {
   kx <- 3
   ky <- 5
   ndim <- min(kx, ky)
   s1 <- cv.pcca(X, Y, kx=kx, ky=ky, nfolds=nfolds)
   s2 <- cv.pcca(X, Y, kx=kx, ky=ky, folds=s1$folds)

   expect_equal(s1$ndim, ndim)
   expect_equal(s1$nfolds, nfolds)
   expect_equivalent(diag(cor(s1$Px, s1$Py)), s1$r)
   expect_equal(s1$folds, s2$folds)

   s3 <- cv.pcca(X, X, kx=kx, ky=ky, folds=s1$folds)
   expect_equivalent(diag(cor(s3$Px, s3$Py)), rep(1, ndim))

   expect_error(cv.pcca(X, NULL))
   expect_error(cv.pcca(NULL, Y))
   expect_error(cv.pcca(X[-1,], Y))
})


