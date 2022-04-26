#' Performs univariate canonical correlation analysis, i.e., ANOVA of all
#' phenotypes on each SNP.
#'
#' @param X An n by p numeric matrix, or a character string pointing to a
#' PLINK dataset
#'
#' @param Y An n by k numeric matrix of phenotypes.
#' 
#' @param standx Character. One of "binom" (zero mean, unit variance
#' where variance is p * (1 - p), for SNP data), "binom2" (zero mean, unit
#' variance where variance is 2 * p * (1 - p),
#' "sd" (zero-mean, unit Gaussian variance), "center" (zero mean), or "none". Note
#' that if you use "binom" for non-SNP data you may get garbage. Note that
#' currently the same standardisation is applied to both X and Y. If you
#' require different standardisations, it's best to standardise yourself
#' and then choose standx="none".
#'
#' @param standy Character. Stanardisation of Y.
#'
#' @param verbose Logical.
#'
#' @param check_geno Logical. Whether to explicitly check if the matrices
#' X and Y contain values other than {0, 1, 2}, when standx or standy is one 
#' "binom" or "binom2". This can
#' be set to FALSE if you are sure your matrices only contain these values
#' (only matters when using "binom"/"binom2").
#'
#' @param check_fam Logical. Whether to check that the number of row in 
#' the PLINK fam file (if X is a character string) matches the number of
#' rows in the eigenvectors.
#' 
#' @details This is an efficient implementation of Ferreira and Purcell's
#' plink.multivariate test for association between multiple phenotypes and one SNP at a time.
#'
#' This test is equivalent to the F-test for a linear regression of each SNP on the
#' phenotypes \code{lm(SNP ~ pheno1 + pheno2 + ...)},
#' hence the number of phenotypes that can be tested is limited by
#' the sample size (the model is not penalised).
#'
#' By default, \code{ucca} performs mean-imputation of missing values.
#'
#' @return \code{ucca} returns a list containing the following components:
#'
#' \describe{  
#'    \item{result:}{A numeric matrix of the association results, in the original ordering.}
#'    \item{npheno:}{The number of phenotypes used.}
#'    \item{nsnps:}{The number of SNPs used.}
#' }
#'
#' @references M. A. R. Ferreira and S. M. Purcell (2009) _A multivariate test
#' of association_, Bioinformatics 25(1) 132-133
#'
#' @examples
#'
#' data(hm3.chr1)
#' n <- nrow(hm3.chr1$bed)
#' p <- ncol(hm3.chr1$bed)
#' k <- 10
#' X <- scale2(hm3.chr1$bed)
#' By <- matrix(rnorm(p * k), p, k)
#' Y <- scale(X %*% By + rnorm(n * k))
#' 
#' # Call on a numeric matrix X;
#' # X has already been standardised
#' s1 <- ucca(X, Y, standx="none", standy="sd")
#' 
#' # Call on a PLINK dataset
#' bedf <- gsub("\\.bed", "",
#'     system.file("extdata", "data_chr1.bed", package="flashpcaR"))
#' s2 <- ucca(bedf, Y, standx="binom2", standy="sd")
#'
#' head(s1$result)
#' head(s2$result)
#'
#' # Same result (to within numerical tolerance)
#' mean((s1$result[, "Fstat"] - s2$result[, "Fstat"])^2)
#'
#' @importFrom utils read.table
#'
#' @export
ucca <- function(X, Y,
   standx=c("binom2", "binom", "sd", "center", "none"),
   standy=c("binom2", "binom", "sd", "center", "none"),
   check_geno=TRUE, check_fam=TRUE, verbose=FALSE)
{
   standx <- match.arg(standx)
   standy <- match.arg(standy)

   if(!is.numeric(Y)) {
      stop("Y must be a numeric matrix")
   } else if(any(is.na(Y))) {
       warning("Y contains missing values, will be mean-imputed")
   }

   Y <- cbind(Y)
   # If Y is an integer matrix, Rcpp will throw exception
   storage.mode(Y) <- "numeric"

   if(is.numeric(X)) {
      X <- cbind(X)
      storage.mode(X) <- "numeric"
      if(any(is.na(X))) {
	 warning("X contains missing values, will be mean-imputed")
      }
   } else if(is.character(X)) {
      if(!standx %in% c("binom", "binom2")) {
	 stop("When using PLINK data, you must use standx='binom' or 'binom2'")
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
   }

   if(is.character(X) && check_fam) {
      fam <- read.table(paste0(X, ".fam"), header=FALSE, sep="",
         stringsAsFactors=FALSE)
      if(ncol(Y) > nrow(fam)) {
         stop(paste(
            "The phenotype matrix Y cannot have more columns than",
            "the sample size"))
      } else if(nrow(Y) != nrow(fam)) {
         stop(paste0("The number of rows in ", X, ".fam and Y don't match"))
      }
      rm(fam)
   } else {
      if(ncol(Y) > nrow(X)) {
	 stop(paste(
	    "The phenotype matrix Y cannot have more columns than",
	    "the sample size"))
      } else if(nrow(Y) != nrow(X)) {
	 stop("The number of rows in X and Y don't match")
      }
   }

   if(is.numeric(X) && standx %in% c("binom", "binom2") && check_geno) {
      wx <- X %in% c(0:2, NA)
      if(sum(wx) != length(X)) {
         stop(
            paste("Your X matrix contains values other than {0, 1, 2},
               standx='binom'/'binom2' can't be used here"))
      }
   }

   if(standy %in% c("binom", "binom2") && check_geno) {
      wy <- Y %in% c(0:2, NA)
      if(sum(wy) != length(Y)) {
         stop(
            paste("Your Y matrix contains values other than {0, 1, 2},
               standy='binom'/'binom2' can't be used here"))
      }
   }

   std <- c(
      "none"=0L,
      "sd"=1L,
      "binom"=2L,
      "binom2"=3L,
      "center"=4L
   )

   standx_i <- std[standx]
   standy_i <- std[standy]

   res <- try(
      if(is.character(X)) {
	 ucca_plink_internal(X, Y, standx_i, standy_i, verbose)
      } else {
	 ucca_internal(X, Y, standx_i, standy_i, verbose)
      }
   )
   class(res) <- "ucca"
   res$npheno <- ncol(Y)
   if(is(res, "try-error")) {
      NULL
   } else {
      if(!is.character(X)) {
	 rownames(res$result) <- colnames(X)
      }
      res
   }
}

#' Prints a UCCA object
#'
#' @param x A flashpca object to be printed
#' @param ... Ignored
#'
#' @export 
print.ucca <- function(x, ...)
{
   cat("ucca object; #phenotypes=",
      x$npheno, ", #SNPs=", nrow(x$result), "\n")
   invisible(x)
}

