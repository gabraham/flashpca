#' Check the accuracy of Principal Component Analysis
#'
#' @param X A numeric matrix to project onto the PCs, or a
#' character string pointing to a PLINK dataset. 
#' 
#' @param evec A numeric matrix of eigenvectors (samples on rows,
#'   ndim dimensions on columns).
#'
#' @param eval A numeric vector of eigenvalues.
#'
#' @param stand A character string indicating how to standardise X before PCA,
#' one of "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style),
#' "sd" (zero-mean unit-variance), "center" (zero mean), or "none".
#'
#' @param divisor A character string indicating whether to divide the
#' eigenvalues by number of columns of X ("p"), the number of 
#' rows of X minus 1 ("n1") or none ("none").
#' 
#' @param block_size Integer. Block size for PCA on PLINK files.
#' 
#' @param verbose Logical. Verbose output.
#' 
#' @param check_geno Logical. Whether to explicitly check if the matrix X
#' contains values other than {0, 1, 2}, when stand="binom". This can be
#' be set to FALSE if you are sure your matrix only contains these values
#' (only matters when using stand="binom").
#'
#' @param check_fam Logical. Whether to check that the number of row in 
#' the PLINK fam file (if X is a character string) matches the number of
#' rows in the eigenvectors.
#' 
#' @details
#' \code{check} computes the accuracy of the eigen-decomposition of XX'/m,
#' defined as
#'    \deqn{[1/(nK) \sum_{k=1}^K ||(1/m) XX' U_k - U_k d^2_k||_F^2]^{1/2}}
#'
#' Note: this definition is sensitive to the magnitude of the eigenvalues and 
#' therefore requires specifiying the same standardisation that was used in
#' the original eigen-decomposition, and the same divisor of XX' (usually the
#' number of SNPs, p).
#'
#' @return \code{check} returns a list containing the following components:
#' \describe{  
#'    \item{err}{A numeric vector. The squared error along each dimension k.}
#'    \item{mse}{A numeric value. The mean squared error; the sum of the
#'	 squared errors divided by N * ndim.}
#'    \item{rmse}{A numeric value. The root mean squared error, the square
#'    root of mse.}
#' }
#'
#' @export
check <- function(X, evec, eval,
   stand=c("binom2", "binom", "sd", "center", "none"),
   divisor="p", block_size=1000, verbose=FALSE,
   check_geno=TRUE, check_fam=TRUE)
{
   stand <- match.arg(stand)
   divisor <- match.arg(divisor)

   eval <- as.numeric(eval)
   
   if(is.numeric(X)) {
      if(any(is.na(X))) {
	 warning("X contains missing values, will be mean imputed")
      }

      if(stand %in% c("binom", "binom2") && check_geno) {
	 wx <- X %in% c(0:2, NA)
      	 if(sum(wx) != length(X)) {
      	    stop(
      	       "Your data contains values other than {0, 1, 2}, ",
      	       "stand='binom'/'binom2' can't be used here")
      	 }
      }

      if(nrow(evec) != nrow(X)) {
	 stop("The number of rows in X and evec don't match")
      }

   } else if(is.character(X)) {
      if(!stand %in% c("binom", "binom2")) {
	 stop("When using PLINK data, you must use stand='binom' or 'binom2'")
      }
      if(check_fam) {
	 fam <- read.table(paste0(X, ".fam"), header=FALSE, sep="",
	    stringsAsFactors=FALSE)
	 if(nrow(fam) != nrow(evec)) {
	    stop(paste0("The number of rows in ", X, ".fam",
	       " and evec don't match"))
	 }
	 rm(fam)
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
   }

   if(ncol(evec) != length(eval)) {
       stop(paste0("The number of columns of evec doesn't",
	    "match the number of eigenvalues eval"))
   }

   divisors <- c(
      "p"=2,
      "n1"=1,
      "none"=0
   )
   div <- divisors[divisor]

   std <- c(
      "none"=0L,
      "sd"=1L,
      "binom"=2L,
      "binom2"=3L,
      "center"=4L
   )
   stand_i <- std[stand]

   # If the matrix is integer, Rcpp will throw an exception
   if(is.numeric(X)) {
      storage.mode(X) <- "numeric"
   }

   res <- try(
      if(is.character(X)) {
	 check_plink_internal(X, stand_i, evec, eval, block_size, div,
	    verbose)
      } else {
	 check_internal(X, stand_i, evec, eval, div, verbose)
      }
   )
   class(res) <- "flashpca.check"
   if(is(res, "try-error")) {
      NULL
   } else {
      res
   }
}

