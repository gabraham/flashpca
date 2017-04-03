#' Performs sparse canonical correlation analysis.
#'
#' @param X An n by p numeric matrix, or a character string pointing to a
#' PLINK dataset.
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
#' @param check_fam Logical. Whether to check that the number of rows in 
#' the PLINK fam file (if X is a character string) matches the number of
#' rows in the eigenvectors.
#' 
#' @param V Numeric. A vector to initialise "v" in SCCA iterations. By
#' default, it will be a vector of normally distributed variates.
#'
#' @param block_size Integer. Size of blocks for reading PLINK data.
#'
#' @param simplify Logical. Whether to return a single \code{scca} object or a
#' list when only one model is fitted.
#'
#'
#' @return \code{scca} returns a list containing the following components:
#'
#' \describe{  
#'    \item{U:}{Top ndim left canonical vectors of X^T Y.}
#'    \item{V:}{Top ndim right canonical vectors of X^T Y.}
#'    \item{d:}{Top ndim canonical covariances, i.e., diag(cov(X U, Y V)). 
#'	 Note that we don't divide by n-1.}
#'    \item{Px:}{X * U.}
#'    \item{Py:}{Y * V.}
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
#' s <- scca(X, Y, lambda1=1e-2, lambda2=1e-2, ndim=5,
#'   standx="none", standy="sd")
#'
#' ## The canonical correlations
#' diag(cor(s$Px, s$Py))
#'
#' @importFrom utils read.table
#'
#' @export
scca <- function(X, Y, lambda1=0, lambda2=0,
   standx=c("binom2", "binom", "sd", "center", "none"),
   standy=c("binom2", "binom", "sd", "center", "none"),
   ndim=10, maxiter=1e3, tol=1e-4, seed=1L, verbose=FALSE, num_threads=1,
   mem=c("low", "high"), check_geno=TRUE, check_fam=TRUE,
   V=NULL, block_size=500, simplify=TRUE)
{
   standx <- match.arg(standx)
   standy <- match.arg(standy)
   mem <- match.arg(mem)

   if(!is.numeric(Y)) {
      stop("Y must be a numeric matrix")
   } else if(any(is.na(Y))) {
       warning("Y cantains missing values, will be mean imputed")
   }

   Y <- cbind(Y)
   # If Y is an integer matrix, Rcpp will throw exception
   storage.mode(Y) <- "numeric"

   if(is.numeric(X)) {
      if(any(is.na(X))) {
	 warning("X cantains missing values, will be mean imputed")
      }

      if(ncol(X) < 2) {
	 stop("X must have at least two columns")
      }

      if(nrow(X) < 2) {
	 stop("X must have at least two rows")
      }

      if(standx %in% c("binom", "binom2") && check_geno) {
	 wx <- X %in% c(0:2, NA)
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
      if(check_fam) {
	 fam <- read.table(paste0(X, ".fam"), header=FALSE, sep="",
	    stringsAsFactors=FALSE)
	 if(nrow(Y) != nrow(fam)) {
	    stop(paste0("The number of rows in ", X,
		  ".fam and Y don't match"))
	 }
	 rm(fam)
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
   }

   if(is.character(X)) {
      fam <- read.table(paste0(X, ".fam"), header=FALSE, sep="",
         stringsAsFactors=FALSE)
      if(nrow(Y) != nrow(fam)) {
	 stop("The number of rows in X and Y don't match")
      }
      rm(fam)
   } else {
      if(nrow(Y) != nrow(X)) {
	 stop("The number of rows in X and Y don't match")
      }
   }

   if(any(lambda1 < 0)) {
      stop("lambda1 must be non-negative")
   }
   if(any(lambda2 < 0)) {
      stop("lambda2 must be non-negative")
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
	 lapply(lambda1, function(l1) {
	    lapply(lambda2, function(l2) {
	       x <- scca_plink_internal(X, Y, l1, l2, ndim,
		  standx_i, standy_i, mem_i, seed, maxiter, tol,
		  verbose, num_threads, block_size, useV, V)
	    })
	 })
      } else {
	 lapply(lambda1, function(l1) {
	    lapply(lambda2, function(l2) {
	       scca_internal(X, Y, l1, l2, ndim,
		  standx_i, standy_i, mem_i, seed, maxiter, tol,
		  verbose, num_threads, useV, V)
	    })
	 })
      }
   )
   if(is(res, "try-error")) {
      NULL
   } else {
      #if(!is.character(X)) {
	 #rownames(res$result) <- colnames(X)
      #}
      if(simplify && length(lambda1) == 1 && length(lambda2) == 1) {
	 x <- res[[1]][[1]]
	 class(x) <- "scca"
	 x
      } else {
	 class(res) <- "scca-list"
	 res
      }
   }
}

#' Prints an SCCA object
#'
#' @param x A flashpca object to be printed
#' @param ... Ignored
#' @export 
print.scca <- function(x, ...)
{
   cat("scca object; ndim=", length(x$d), "\n")
   invisible(x)
}

#' Cross-validated grid search over SCCA penalties
#'
#' @param X A numeric matrix. The use of PLINK datasets is currently not
#' supported.
#'
#' @param Y A numeric matrix
#'
#' @param lambda1 A numeric vector of non-negative penalties
#'
#' @param lambda2 A numeric vector of non-negative penalties
#'
#' @param ndim Integer. The number of dimensions to infer.
#'
#' @param nfolds Integer. The number of cross-validation folds.
#'
#' @param parallel Logical. Whether to parallelise the cross-validation using
#' the foreach package.
#'
#' @param ... Other arguments that will be passed to \code{scca}.
#' 
#' @details 
#' Note that the default penalties may not be appropriate for every dataset
#' and for some values the algorithm may not converge (especially for small
#' penalties).
#'
#' @return \code{cv.scca} returns an array containing the average
#' cross-validated canonical Pearson correlations between X and Y, as follows:
#'    Dimension 1: the ndim different canonical dimensions;
#'    Dimension 2: along the lambda1 penalties;
#'    Dimension 3: along the lambda2 penalties.
#'
#' @importFrom stats cor
#'
#' @export
cv.scca <- function(X, Y,
   lambda1=seq(1e-6, 1e-3, length=5), lambda2=seq(1e-6, 1e-3, length=5),
   ndim=3, nfolds=10, parallel=FALSE, ...)
{
   n <- nrow(Y)
   if(is.character(X)) {
      stop(
	 "Cross-validation currently only supported when X is a numeric matrix")
   }

   folds <- sample(1:nfolds, n, replace=TRUE)
   xpred <- array(0, c(n, ndim, length(lambda1), length(lambda2)))
   ypred <- array(0, c(n, ndim, length(lambda1), length(lambda2)))
   nzx <- array(0, c(ndim, length(lambda1), length(lambda2)))
   nzy <- array(0, c(ndim, length(lambda1), length(lambda2)))

   # https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Suggested-packages
   if(parallel && requireNamespace("foreach", quietly=TRUE)) {
      res <- foreach::`%dopar%`(foreach::foreach(fold=1:nfolds), {
         w <- folds != fold
         scca(X[w,], Y[w,], ndim=ndim,
	    lambda1=lambda1, lambda2=lambda2,
	    simplify=FALSE, ...)
      })
   } else {
      res <- lapply(1:nfolds, function(fold) {
         w <- folds != fold
	 scca(X[w,], Y[w,], ndim=ndim,
	    lambda1=lambda1, lambda2=lambda2,
	    simplify=FALSE, ...)
      })
   }

   # Collect all the folds into one result (same way as glmnet), instead
   # of averaging over k folds
   for(fold in 1:nfolds) {
      for(i in seq(along=lambda1)) {
	 for(j in seq(along=lambda2)) {
	    x <- res[[fold]][[i]][[j]]
	    w <- folds == fold
	    xpred[w, ,i, j] <- X[w,] %*% x$U
	    ypred[w, ,i, j] <- Y[w,] %*% x$V
	 }
      }
   }

   for(i in seq(along=lambda1)) {
      for(j in seq(along=lambda2)) {
	 mx <- sapply(res, function(x) colSums(x[[i]][[j]]$U != 0, na.rm=TRUE))
	 my <- sapply(res, function(x) colSums(x[[i]][[j]]$V != 0, na.rm=TRUE))
	 nzx[,i,j] <- rowMeans(rbind(mx))
	 nzy[,i,j] <- rowMeans(rbind(my))
      }
   }

   r <- lapply(1:ndim, function(k) {
      sapply(seq(along=lambda1), function(i) {
	 sapply(seq(along=lambda2), function(j) {
	    cor(xpred[,k,i,j], ypred[,k,i,j])
	 })
      })
   })
   r <- abind(r, along=3)
   r <- aperm(r, c(3, 2, 1))

   res2 <- list(
      ndim=ndim,
      lambda1=lambda1,
      lambda2=lambda2,
      corr=r,
      nzero.x=nzx,
      nzero.y=nzy
   )
   class(res2) <- "cv.scca"
   res2
}

#' Plot the results of SCCA cross-validation
#'
#' @param x An object of class "cv.scca"
#'
#' @param dim Integer. Which dimension to plot (all will be plotted by
#' default).
#' 
#' @param ... Other arguments that will be passed to \code{plot}.
#'
#' @details
#'    Plots the cross-validated Pearson correlation, as a function of the
#'    number of non-zero entries in the canonical vectors U (on the x-axis),
#'    with a separate curve separately for each lambda2 penalty (which affects
#'    the non-zero entries in the canonical vector V. For example, if
#'    \code{cv.scca} was run with 5x lambda1 values and 7x lambda2 values,
#'    there will be 7 separate curves each with 5 points on the x-axis.
#'
#' @examples
#'
#' #######################
#' ## HapMap3 chr1 example
#' data(hm3.chr1)
#' X <- scale2(hm3.chr1$bed)
#' n <- nrow(X)
#' m <- ncol(X)
#' k <- 5
#' B <- matrix(rnorm(m * k), m, k)
#' Y <- X %*% B + rnorm(n * k)
#'
#' r <- cv.scca(X, Y, standx="sd", standy="sd", nfolds=3, ndim=2,
#'    lambda1=seq(1e-3, 1e-1, length=10),
#'    lambda2=seq(1e-4, 0.5, length=3))
#'
#' par(mfrow=c(1, 2))
#' plot(r, dim=1)
#' plot(r, dim=2)
#' 
#' @importFrom abind abind
#' @importFrom graphics matplot
#' @importFrom graphics legend
#'
#' @export
plot.cv.scca <- function(x, dim=NULL, ...)
{
   if(class(x) != "cv.scca") {
      stop("x is not of class 'cv.scca'")
   }

   if(is.null(dim)) {
      for(i in 1:x$ndim) {
         matplot(x$nzero.x[1,,], x$corr[i,,],
            type="b", xlab="# of variables in X with non-zero weight",
               ylab="Pearson correlation", log="x",
               main=paste0("Canonical dimension ", i))
      }
   } else if(is.numeric(dim) && dim > 0 && dim <= x$ndim) {
      matplot(x$nzero.x[dim,,], x$corr[dim,,],
	 type="b", xlab="# of variables in X with non-zero weight",
	 ylab="Pearson correlation", log="x",
	 main=paste0("Canonical dimension ", dim))
   } else {
      stop(paste("dim must be a positive integer, <=", x$ndim))
   }
   legend("bottomright", legend=round(x$lambda2, 6),
      title="lambda2", lwd=2, lty=1:length(x$lambda2),
      col=1:length(x$lambda2))
   invisible(x)
}

