#' Principal Component Analysis using FlashPCA
#'
#' @param X A numeric matrix to perform PCA on, or a
#' character string pointing to a PLINK dataset. 
#'
#' @param ndim Integer. How many dimensions to return in results.
#'
#' @param stand A character string indicating how to standardise X before PCA,
#' one of "binom" (old Eigenstrat-style), "binom2" (new Eigenstrat-style),
#' "sd" (zero-mean unit-variance), "center" (zero mean), or "none".
#' 
#' @param divisor A character string indicating whether to divide the
#' eigenvalues by number of columns of X ("p"), the number of 
#' rows of X minus 1 ("n1") or none ("none").
#'
#' @param maxiter Integer. How many iterations to use in PCA.
#'
#' @param tol Numeric. Tolerance for determining PCA convergence.
#'
#' @param seed Integer. Seed for random number generator.
#'
#' @param block_size Integer. Block size for PCA on PLINK files.
#'
#' @param verbose logical. More verbose output.
#'
#' @param do_loadings Logical. Whether to compute loadings (matrix V from SVD).
#'
#' @param check_geno Logical. Whether to explicitly check if the matrix X
#' contains values other than {0, 1, 2}, when stand="binom". This can be
#' be set to FALSE if you are sure your matrix only contains these values
#' (only matters when using stand="binom").
#' 
#' @param return_scale Logical. Whether to return the means and standard
#' deviations used in standardizing the matrix X.
#'
#' @details
#'    The default decomposition is of X X' / m, where m is the number of SNPs
#'    (the denominator can be changed using the 'divisor' argument).
#'
#'    The data is standardised by default. \code{stand = "binom"} uses the old Eigensoft
#'    (Price 2006) formulation of
#'	       \deqn{p_j = sum_{i=1}^n X_{i,j} / (2 * n)}
#'	       \deqn{mean_j = 2 * p}
#'	       \deqn{sd_j = sqrt(p * (1 - p)),}
#'    where j is the index for the SNP and i is the index for the sample.
#'    Alternatively, `stand = "binom2"' uses the newer formula, which is
#'    similar to the above except for
#'	       \deqn{sd_j = sqrt(2 * p * (1 - p)).}
#' 
#' @return \code{flashpca} returns a list containing the following components:
#' \describe{  
#'    \item{values:}{a numeric vector. The eigenvalues of X X' / m.}
#'    \item{vectors:}{a numeric matrix. The eigenvectors of X X' / m.}
#'    \item{projection:}{a numeric matrix. Equivalent to X V.}
#'    \item{loadings:}{a numeric matrix. The matrix of variable loadings, i.e., V
#'    from SVD.}
#'    \item{scale:}{a list of two elements, ``center'' and ''scale'', which was
#'	 used to standardise the input matrix X.}
#' }
#' 
#' @examples
#' 
#' set.seed(123)
#' 
#' #######################
#' ## HapMap3 chr1 example
#' data(hm3.chr1)
#' ndim <- 10
#' X <- scale2(hm3.chr1$bed)
#' f1 <- flashpca(X, ndim=ndim, stand="none")
#'
#' # prcomp's is too slow for this example
#' r <- eigen(tcrossprod(X) / ncol(X), symmetric=TRUE)
#' 
#' ## Compare eigenvalues
#' eval <- cbind(r$values[1:ndim], f1$values)
#' cor(eval)
#' mean((eval[,1] - eval[,2])^2)
#' 
#' ## Compare eigenvectors
#' diag(cor(r$vectors[, 1:ndim], f1$vectors))
#'
#' ####################################
#' # HapMap3 chr1 example, PLINK format
#'
#' bedf <- gsub("\\.bed", "",
#'    system.file("extdata", "data_chr1.bed", package="flashpcaR"))
#'
#' f2 <- flashpca(bedf, ndim=ndim) 
#' 
#' eval <- cbind(r$values[1:ndim], f2$values)
#' cor(eval)
#' mean((eval[,1] - eval[,2])^2)
#' 
#' ## Compare eigenvectors
#' diag(cor(r$vectors[, 1:ndim], f2$vectors))
#'
#' @export
flashpca <- function(X, ndim=10,
   stand=c("binom2", "binom", "sd", "center", "none"),
   divisor=c("p", "n", "none"),
   maxiter=1e2, tol=1e-4, seed=1, block_size=1000, verbose=FALSE,
   do_loadings=FALSE, check_geno=TRUE, return_scale=TRUE)
{
   stand <- match.arg(stand)
   divisor <- match.arg(divisor)
   
   if(is.numeric(X)) {
      if(any(is.na(X))) {
	 warning("X contains missing values, will be mean imputed")
      }

      if(ncol(X) < 2) {
	 stop("X must have at least two columns")
      }

      if(nrow(X) < 2) {
	 stop("X must have at least two rows")
      }

      if(stand %in% c("binom", "binom2") && check_geno) {
	 wx <- X %in% c(0:2, NA)
      	 if(sum(wx) != length(X)) {
      	    stop(
      	       "Your data contains values other than {0, 1, 2}, ",
      	       "stand='binom'/'binom2' can't be used here")
      	 }
      }

      block_size <- min(block_size, ncol(X))
   } else if(is.character(X)) {
      if(!stand %in% c("binom", "binom2")) {
	 stop("When using PLINK data, you must use stand='binom' or 'binom2'")
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
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

   if(is.numeric(X)) {
      maxdim <- min(dim(X))
      ndim <- min(maxdim, ndim)
   }

   # If the matrix is integer, Rcpp will throw an exception
   if(is.numeric(X)) {
      storage.mode(X) <- "numeric"
   }

   res <- try(
      if(is.character(X)) {
	 flashpca_plink_internal(X, stand_i, ndim, div,
	    maxiter, block_size, tol, seed,
	    verbose, do_loadings, return_scale)
      } else {
	 flashpca_internal(X, stand_i, ndim, div, maxiter,
	    tol, seed, verbose, do_loadings, return_scale)
      }
   )
   class(res) <- "flashpca"
   if(is(res, "try-error")) {
      NULL
   } else {
      res
   }
}

#' Prints a flashpca object
#'
#' @param x A flashpca object to be printed
#' @param ... Ignored.
#'
#' @export 
print.flashpca <- function(x, ...)
{
   cat("flashpca object; ndim=", length(x$values), "\n")
   invisible(x)
}

