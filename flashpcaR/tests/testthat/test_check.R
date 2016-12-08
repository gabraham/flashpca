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

   c1 <- check(S, stand="none", evec=f2$vectors, eval=f2$values)

})

