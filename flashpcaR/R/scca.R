#' Performs sparse canonical correlation analysis.
#'
#' @param X An n by p numeric matrix
#'
#' @param Y An n by k numeric matrix
#' 
#' @param lambda1 Numeric. Non-negative L1 penalty on canonical vectors of X.
#'
#' @param lambda2 Numeric. Non-negative L1 penalty on canonical vectors of Y.
#"
#' @param standx Character. One of "binom" (zero mean, unit variance
#' where variance is p * (1 - p), for SNP data), "binom2" (zero mean, unit
#' variance where variance is 2 * p * (1 - p),
#' "sd" (zero-mean, unit Gaussian variance), "center" (zero mean), or "none". Note
#' that if you use "binom" for non-SNP data you may get garbage. Note that
#' currently the same standardization is applied to both X and Y. If you
#' require different standardizations, it's best to standardise yourself
#' and then choose standx="none". When X is a PLINK dataset name, standx must be one of "binom2" 
#' or "binom".
#'
#' @param standy Character. Stanardisation of Y.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param maxiter Integer. Positive, maximum number of iterations to perform.
#'
#' @param tol Numeric. Tolerance for convergence.
#'
#' @param seed Integer. Random seed for initialisation of CCA iterations.
#'
#' @param verbose Logical.
#'
#' @param num_threads Integer. Number of OpenMP threads to
#' use (only supported under Linux and on Mac if compiled with GCC).
#'
#' @param mem Character. One of "low" or "high", where "low"
#' doesn't explicitly compute X^T Y, and "high" does.
#' This is useful for large X and Y.
#'
#' @param check_geno Logical. Whether to explicitly check if the matrices
#' X and Y contain values other than {0, 1, 2}, when standx/standy="binom"/"binom2". This can
#' be set to FALSE if you are sure your matrices only contain these values
#' (only matters when using "binom" or "binom2").
#' 
#' @param V Numeric. A vector to initialise "v" in SCCA iterations. By
#' default, it will be a vector of normally distributed variates.
#'
#' @param blocksize Integer. Size of blocks for reading PLINK data.
#'
#' @return \code{scca} returns a list containing the following components:
#'
#' \describe{  
#'    \item{U}{Top ndim left canonical vectors of X^T Y.}
#'    \item{V}{Top ndim right canonical vectors of X^T Y.}
#'    \item{d}{Top ndim canonical covariances, i.e., diag(cov(X U, Y V)). 
#'	 Note that we don't divide by n-1.}
#'    \item{Px}{X times U.}
#'    \item{Py}{Y times V.}
#' }
#'
#' @examples
#'
#' #######################
#' ## HapMap3 chr1 example
#' data(hm3.chr1)
#' X <- scale2(hm3.chr1$bed)
#' n <- nrow(X)
#' m <- ncol(X)
#' k <- 10
#' B <- matrix(rnorm(m * k), m, k)
#' Y <- X %*% B + rnorm(n * k)
#'
#' s <- scca(X, Y, lambda1=1e-2, lambda2=1e-2, ndim=5, standy="sd")
#'
#' ## The canonical correlations
#' diag(cor(s$Px, s$Py))
#'
#' @export
scca <- function(X, Y, lambda1=0, lambda2=0,
   standx=c("binom2", "binom", "sd", "center", "none"),
   standy=c("binom2", "binom", "sd", "center", "none"),
   ndim=10, maxiter=1e3, tol=1e-4, seed=1L, verbose=FALSE, num_threads=1,
   mem=c("low", "high"), check_geno=TRUE, V=NULL, block_size=500)
{
   standx <- match.arg(standx)
   standy <- match.arg(standy)
   mem <- match.arg(mem)

   if(!is.numeric(Y)) {
      stop("Y must be a numeric matrix")
   } else if(any(is.na(Y))) {
       stop("Y cannot contain any missing values")
   }

   Y <- cbind(Y)
   # If Y is an integer matrix, Rcpp will throw exception
   storage.mode(Y) <- "numeric"

   if(is.numeric(X)) {
      if(any(is.na(X))) {
	 stop("X cannot contain any missing values")
      }

      if(ncol(X) < 2) {
	 stop("X must have at least two columns")
      }

      if(nrow(X) < 2) {
	 stop("X must have at least two rows")
      }

      if(standx %in% c("binom", "binom2") && check_geno) {
	 wx <- X %in% 0:2
      	 if(sum(wx) != length(X)) {
      	    stop(
      	       "Your data contains values other than {0, 1, 2}, ",
      	       "standx='binom'/'binom2' can't be used here")
      	 }
      }
   } else if(is.character(X)) {
      if(!standx %in% c("binom", "binom2")) {
	 stop("When using PLINK data, you must use standx='binom' or 'binom2'")
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
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

   #maxdim <- min(dim(X), dim(Y))
   #ndim <- min(maxdim, ndim)

   if(mem == "high") {
      mem_i <- 2L
   } else {
      mem_i <- 1L
   }

   if(!is.null(V)) {
      V <- cbind(V)
      if(nrow(V) < ncol(Y) || ncol(V) != ndim) {
         stop("dimensions of V must be (ncol(Y) x (ndim))")
      }
      useV <- TRUE
   } else 
   {
      V <- matrix(0.)
      useV <- FALSE
   }


   res <- try(
      if(is.character(X)) {
	 scca_plink_internal(X, Y, lambda1, lambda2, ndim,
	    standx_i, standy_i, mem_i, seed, maxiter, tol,
	    verbose, num_threads, block_size, useV, V)
      } else {
	 scca_internal(X, Y, lambda1, lambda2, ndim,
	    standx_i, standy_i, mem_i, seed, maxiter, tol,
	    verbose, num_threads, useV, V)
      }
   )
   class(res) <- "ucca"
   if(is(res, "try-error")) {
      NULL
   } else {
      #if(!is.character(X)) {
	 #rownames(res$result) <- colnames(X)
      #}
      res
   }
}

#' Prints an SCCA object
#'
#' @param x A flashpca object to be printed
#' @export 
print.scca <- function(x, ...)
{
   cat("scca object; ndim=", length(x$d), "\n")
   invisible(x)
}

cv.scca <- function(X, Y, nfolds=10, parallel=FALSE, ...)
{
   n <- nrow(Y)
   if(is.character(X)) {
      stop("Cross-validation currently only supported when X is a matrix")
   }

   folds <- sample(1:nfolds, n, replace=TRUE)

   if(parallel) {
      require(foreach)
      res <- foreach(fold=1:nfolds) %dopar% {
         w <- folds != fold
         s <- scca(X[w,], Y[w,], ...)
         px <- X[!w,] %*% s$U
         py <- Y[!w,] %*% s$V
         diag(cor(px, py))
      }
   } else {
      res <- sapply(1:nfolds, function(fold) {
         w <- folds != fold
         s <- scca(X[w,], Y[w,], ...)
         px <- X[!w,] %*% s$U
         py <- Y[!w,] %*% s$V
         diag(cor(px, py))
      })
   }
   res
}

