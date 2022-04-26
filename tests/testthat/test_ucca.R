context("Testing UCCA")

n <- 500
p <- 100
k <- 15
m <- 10
M <- matrix(rnorm(n * m), n, m)
Bx <- matrix(rnorm(m * p), m, p)
By <- matrix(rnorm(m * k), m, k)
X0 <- scale(M %*% Bx + rnorm(n * p))
Y0 <- scale(M %*% By + rnorm(n * k))

data(hm3.chr1)

bedf <- gsub("\\.bed", "",
   system.file("extdata", "data_chr1.bed", package="flashpcaR"))

test.tol <- 1e-4

test <- function(s, X, Y)
{
   p <- ncol(X)

   r2.obs <- s$result[,"R"]^2
   r2.exp <- numeric(p)

   F.obs <- s$result[, "Fstat"]
   F.exp <- numeric(p)
   
   p.obs <- s$result[, "P"]
   p.exp <- numeric(p)

   # expect_equal also tests names
   names(r2.exp) <- names(F.exp) <- names(p.exp) <- colnames(X)

   for(i in 1:p) {
      l <- lm(X[,i] ~ Y)
      l2 <- summary(l)
      r2.exp[i] <- l2$r.squared
      l3 <- anova(l)
      F.exp[i] <- l3[1, "F value"]
      p.exp[i] <- l3[1, "Pr(>F)"]
   }

   expect_equal(r2.exp, r2.obs, tol=test.tol)
   expect_equal(F.exp, F.obs, tol=test.tol)
   expect_equal(p.exp, p.obs, tol=test.tol)
   expect_equal(log(p.exp), log(p.obs), tol=test.tol)
}

test_that("Testing UCCA (binomial) with matrices and PLINK", {
   X <- scale2(hm3.chr1$bed, type="1")
	
   n <- nrow(X)
   p <- ncol(X)
   By <- matrix(rnorm(p * k), p, k)
   Y <- scale(X %*% By + rnorm(n * k))

   s1 <- ucca(X, Y, standx="none", standy="none")
   s2 <- ucca(bedf, Y, standx="binom", standy="none")
   test(s1, X, Y)
   test(s2, X, Y)
})

test_that("Testing UCCA (binomial2) with matrices and PLINK", {
   X <- scale2(hm3.chr1$bed, type="2")
	
   n <- nrow(X)
   p <- ncol(X)
   By <- matrix(rnorm(p * k), p, k)
   Y <- scale(X %*% By + rnorm(n * k))

   s1 <- ucca(X, Y, standx="none", standy="none")
   s2 <- ucca(bedf, Y, standx="binom2", standy="none")
   test(s1, X, Y)
   test(s2, X, Y)
})

test_that("Testing input checking", {

   X <- scale2(hm3.chr1$bed, type="1")

   # Test incompatible number of rows
   Z <- matrix(rnorm((nrow(X) + 3) * 10), nrow(X) + 3, 10)

   expect_error(
      ucca(X, Z, standx="none", standy="none", mem="high")
   )
   expect_error(
      ucca(X, Z, standx="none", standy="none", mem="low")
   )
   expect_error(
      ucca(bedf, Z, standx="binom2", standy="none")
   )

   # Test having too many phenotypes
   Z <- matrix(rnorm(nrow(X) * (nrow(X) + 3)), nrow(X), nrow(X) + 3)
   expect_error(
      ucca(X, Z, standx="none", standy="none", mem="high")
   )
   expect_error(
      ucca(X, Z, standx="none", standy="none", mem="low")
   )
   expect_error(
      ucca(bedf, Z, standx="binom2", standy="none")
   )

})


