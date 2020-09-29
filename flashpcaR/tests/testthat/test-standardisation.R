context("Testing standardisation and mean-imputation")


test_that("Testing standardisation", {
   n <- 50
   m <- 10

   X <- matrix(rbinom(n * m, size=2, prob=0.3), n, m)
   storage.mode(X) <- "numeric"

   ################################################################################
   # No missing values

   # No missing values, no standardisation
   s1 <- standardise_impute(X, 0)
   expect_equal(s1, X)

   # No missing values, just standardisation to unit sd
   s2 <- standardise_impute(X, 1)
   expect_equal(numeric(m), colMeans(s2))
   expect_equal(numeric(m) + 1, apply(s2, 2, sd))

   # No missing values, just standardisation to unit sd
   s3 <- standardise_impute(X, 1)
   expect_equal(numeric(m), colMeans(s3))

   # No missing values, just standardisation to unit sd (binom)
   s4 <- standardise_impute(X, 2)
   expect_equal(numeric(m), colMeans(s4))
   s4.exp <- scale2(X, type="1")
   expect_equal(as.numeric(s4), as.numeric(s4.exp))

   # No missing values, just standardisation to unit sd (binom2)
   s5 <- standardise_impute(X, 3)
   expect_equal(numeric(m), colMeans(s5))
   s5.exp <- scale2(X, type="2")
   expect_equal(as.numeric(s5), as.numeric(s5.exp))

   # No missing values, just center
   s6 <- standardise_impute(X, 4)
   expect_equal(numeric(m), colMeans(s6))
   s6.exp <- scale(X, center=TRUE, scale=FALSE)
   expect_equal(as.numeric(s6), as.numeric(s6.exp))

   
   ################################################################################
   # Missing values, impute to mean of each column
   
   # Missing values, no standardisation, just impute
   # to mean of each column
   X2 <- X
   diag(X2) <- NA
   Xmean <- colMeans(X2, na.rm=TRUE)
   s7 <- standardise_impute(X2, 0)
   expect_equal(Xmean, diag(s7))
   expect_equal(Xmean, colMeans(s7)) # shouldn't be changed
   expect_that(sum(is.na(s7)), equals(0))

   # Missing values, standardise to SD (binom), impute
   # to mean of each column
   s8 <- standardise_impute(X2, 1)
   expect_equal(numeric(m), colMeans(s8))
   s8.exp <- scale(X2)
   s8.exp[is.na(s8.exp)] <- 0
   expect_equal(as.numeric(s8), as.numeric(s8.exp))

   # Missing values, standardise to SD (binom), impute
   # to mean of each column
   s8 <- standardise_impute(X2, 2)
   expect_equal(numeric(m), colMeans(s8))
   s8.exp <- scale2(X2, type="1")
   expect_equal(as.numeric(s8), as.numeric(s8.exp))

   # Missing values, standardise to SD (binom2), impute
   # to mean of each column
   s9 <- standardise_impute(X2, 3)
   expect_equal(numeric(m), colMeans(s9))
   s9.exp <- scale2(X2, type="2")
   expect_equal(as.numeric(s9), as.numeric(s9.exp))

   ## No missing values, just center
   s10 <- standardise_impute(X2, 4)
   expect_equal(numeric(m), colMeans(s10))
   s10.exp <- scale(X2, center=TRUE, scale=FALSE)
   s10.exp[is.na(s10.exp)] <- 0
   expect_equal(as.numeric(s10), as.numeric(s10.exp))
})

