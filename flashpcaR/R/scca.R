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
#'
#' @param standx Character. One of "binom" (zero mean, unit variance
#' where variance is p * (1 - p), for SNP data), "binom2" (zero mean, unit
#' variance where variance is 2 * p * (1 - p),
#' "sd" (empricial standard deviation), "center" (zero mean), or "none". Note
#' that if you use "binom" for non-SNP data you may get garbage.
#' When X is a PLINK dataset name, standx must be one of "binom3" 
#' or "binom".
#'
#' @param standy Character. Stanardisation of Y.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param divisor A character string indicating whether to divide the
#' eigenvalues by number of rows of X minus 1 ("n1") or none ("none").
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
#' @param check_geno Logical. Whether to explicitly check if the matrices
#' X and Y contain values other than {0, 1, 2}, when standx/standy="binom"/"binom2". This can
#' be set to FALSE if you are sure your matrices only contain these values
#' (only matters when using "binom" or "binom2").
#'
#' @param check_fam Logical. Whether to check that the number of rows in 
#' the PLINK fam file (if X is a character string) matches the number of
#' rows in the eigenvectors.
#'
#' @param check_bim Logical. Whether to use the PLINK bim file
#' (if X is a character string) information to name the vectors.
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
#' s1 <- scca(X, Y, lambda1=1e-2, lambda2=1e-2, ndim=5,
#'   standx="none", standy="sd")
#'
#' ## The canonical correlations
#' diag(cor(s1$Px, s1$Py))
#'
#' ## Example with PLINK-format SNP genotype data
#' f <- system.file("extdata", "data_chr1.bed", package="flashpcaR")
#' s2 <- scca(gsub("\\.bed", "", f), Y, lambda1=1e-2, lambda2=1e-2,
#'    ndim=3, standx="binom2", standy="sd")
#'
#' @importFrom utils read.table
#' @importFrom stats rnorm
#' 
#' @seealso \code{\link{cv.scca}}
#'
#' @export
scca <- function(X, Y, lambda1=0, lambda2=0,
   standx=c("binom2", "binom", "sd", "center", "none"),
   standy=c("binom2", "binom", "sd", "center", "none"),
   ndim=10, divisor=c("n1", "none"),
   maxiter=1e3, tol=1e-4, seed=1L, verbose=FALSE, num_threads=1,
   check_geno=TRUE, check_fam=TRUE, check_bim=TRUE,
   V=NULL, block_size=500, simplify=TRUE)
{
   standx <- match.arg(standx)
   standy <- match.arg(standy)
   divisor <- match.arg(divisor)

   divisors <- c(
      "n1"=1L,
      "none"=0L
   )
   div <- divisors[divisor]

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

      if(nrow(Y) != nrow(X)) {
         stop("The number of rows in X and Y don't match")
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
      }
      if(check_bim) {
	 bim <- read.table(paste0(X, ".bim"), header=FALSE, sep="",
	    stringsAsFactors=FALSE)
	 #if(nrow(Y) != nrow(fam)) {
	 #   stop(paste0("The number of rows in ", X,
	 #	  ".fam and Y don't match"))
	 #}
	 #rm(fam)
      }
   } else {
      stop("X must be a numeric matrix or a string naming a PLINK fileset")
   }

   #if(is.character(X)) {
   #   fam <- read.table(paste0(X, ".fam"), header=FALSE, sep="",
   #      stringsAsFactors=FALSE)
   #   if(nrow(Y) != nrow(fam)) {
   #      stop("The number of rows in X and Y don't match")
   #   }
   #   rm(fam)
   #} else {
   #   if(nrow(Y) != nrow(X)) {
   #      stop("The number of rows in X and Y don't match")
   #   }
   #}

   if(is.null(lambda1) || any(lambda1 < 0)) {
      stop("lambda1 must be non-negative")
   }
   if(is.null(lambda2) || any(lambda2 < 0)) {
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

   if(ndim < 1) {
      stop("ndim can't be less than 1")
   }
   max_dim <- (min(ncol(X), nrow(X), ncol(Y), nrow(Y)))
   if(ndim > max_dim) {
      msg <- paste("You asked for ", ndim," dimensions, but only ",
	 as.integer(max_dim), " allowed", sep="")
      stop(msg)
   }

   if(!is.null(V)) {
      V <- cbind(V)
      if(nrow(V) != ncol(Y) || ncol(V) != ndim) {
         stop("dimensions of V must be (ncol(Y) x (ndim))")
      }
      useV <- TRUE
   } else {
      # Initialise the requested SCCA with a very mildly penalised SCCA, whih
      # is equivalent to SVD of crossprod(X, Y)
      if(verbose) {
	 cat("initialising V\n")
      }
      V <- matrix(0.) # initialised to Gaussian internally
      useV <- FALSE
      l1 <- 1e-9
      l2 <- 1e-9
      if(is.character(X)) {
	 x0 <- try(scca_plink_internal(X, Y, l1, l2, ndim,
	    standx_i, standy_i, div, seed, maxiter, tol,
	     verbose, num_threads, block_size, useV, V))
      } else {
	 x0 <- try(scca_internal(X, Y, l1, l2, ndim,
	    standx_i, standy_i, div, seed, maxiter, tol,
	    verbose, num_threads, useV, V))
      }
      if(is(x0, "try-error")) {
	 return(NULL)
      }
      V <- x0$V
   }

   useV <- TRUE

   # Convenience function to set correct dimnames
   scca_plink_internal2 <- function(X, Y, l1, l2, ndim, standx_i, standy_i,
      div, seed, maxiter, tol, verbose, num_threads, useV, V) {
      s <- scca_plink_internal(X, Y, l1, l2, ndim,
	 standx_i, standy_i, div, seed, maxiter, tol,
	 verbose, num_threads, block_size, useV, V)
      if(check_bim) {
	 rownames(s$U) <- bim$V2
      }
      if(check_fam) {
	 rownames(s$Px) <- paste0(fam$V1, ":", fam$V2)
      }
      if(!is.null(colnames(Y))) {
	 rownames(s$V) <- colnames(Y)
      }
      if(!is.null(rownames(Y))) {
	 rownames(s$Py) <- rownames(Y)
      }
      s
   }

   # Convenience function to set correct dimnames
   scca_internal2 <- function(...) {
      s <- scca_internal(...)
      if(!is.null(colnames(X))) {
	 rownames(s$U) <- colnames(X)
      }
      if(!is.null(colnames(Y))) {
	 rownames(s$V) <- colnames(Y)
      }
      if(!is.null(rownames(X))) {
	 rownames(s$Px) <- rownames(X)
      }
      if(!is.null(rownames(Y))) {
	 rownames(s$Py) <- rownames(Y)
      }
      s
   }

   func <- ifelse(is.character(X), scca_plink_internal2, scca_internal2)

   res <- try(
      lapply(lambda1, function(l1) {
         lapply(lambda2, function(l2) {
            s <- func(X, Y, l1, l2, ndim,
               standx_i, standy_i, div, seed, maxiter, tol,
               verbose, num_threads, useV, V)
	    class(s) <- "scca"
	    s
         })
      })
   )
   if(is(res, "try-error")) {
      NULL
   } else {
      #if(!is.character(X)) {
	 #rownames(res$result) <- colnames(X)
      #}
      if(simplify && length(lambda1) == 1 && length(lambda2) == 1) {
	 #x <- res[[1]][[1]]
	 #class(x) <- "scca"
	 #x
	 res[[1]][[1]]
      } else {
	 class(res) <- "scca-list"
	 res
      }
   }
}

#' Prints an scca object
#'
#' @param x An scca object to be printed
#' @param ... Ignored
#' @export 
print.scca <- function(x, ...)
{
   cat("scca object; ndim=", x$ndim, "\n")
   if(!x$converged) {
      cat("warning: model has not converged\n")
   }
   invisible(x)
}

#' Prints a cv.scca object
#'
#' @param x A cv.scca object to be printed
#' @param ... Ignored
#' @export 
print.cv.scca <- function(x, ...)
{
   cat(paste0(
      "cv.scca object; ndim=", x$ndim, "; ", x$nfolds,
      "-fold cross-validation\n"))
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
#' @param folds Integer vector. The fold identifiers.
#'
#' @param opt.dim Integer. Which dimension to select the optimal penalties
#' for.
#'
#' @param parallel Logical. Whether to parallelise the cross-validation using
#' the foreach package.
#' 
#' @param init Logical. Whether to initialise SCCA with the SVD of X'Y.
#'
#' @param verbose Logical. Whether to print more information.
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
#'    lambda2=seq(1e-4, 0.5, length=8))
#'
#' par(mfrow=c(1, 2))
#' plot(r, dim=1)
#' plot(r, dim=2)
#' 
#' @importFrom abind abind
#' @importFrom graphics matplot
#' @importFrom graphics legend
#' @importFrom stats rnorm
#'
#' @seealso \code{\link{scca}}
#'
#' @export
cv.scca <- function(X, Y,
   lambda1=seq(1e-6, 1e-3, length=5), lambda2=seq(1e-6, 1e-3, length=5),
   ndim=1, nfolds=10, folds=NULL, opt.dim=1, parallel=FALSE, init=TRUE,
   verbose=FALSE, return.models=FALSE, return.folds=FALSE, ...)
{
   n <- nrow(Y)
   if(is.character(X)) {
      stop(
	 "Cross-validation currently only supported when X is a numeric matrix")
   }

   if(nfolds > nrow(Y)) {
      stop("nfolds is too large for the number of samples")
   }

   if(opt.dim <=0 || opt.dim > ndim) {
      stop("opt.dim must be between 1 and ndim")
   }

   if(!is.logical(init)) stop("init muct be TRUE or FALSE")

   if(!is.null(folds)) {
      folds <- as.integer(folds)
      if(length(folds) != n) {
	 stop(
	    "'folds' must be of same number of rows as X and Y")
      }
      if(any(diff(sort(folds)) > 1)) {
	 stop("'folds' must be a set of contiguous integers from 1 to nfolds")
      }
      warning("'folds' will override 'nfolds' parameter")
      nfolds <- max(folds)
   }
   else {
      folds <- sample(1:nfolds, n, replace=TRUE)
   }

   xpred <- array(0, c(n, ndim, length(lambda1), length(lambda2)))
   ypred <- array(0, c(n, ndim, length(lambda1), length(lambda2)))
   nzx <- array(0, c(ndim, length(lambda1), length(lambda2)))
   nzy <- array(0, c(ndim, length(lambda1), length(lambda2)))
   converged <- array(FALSE, c(nfolds, length(lambda1), length(lambda2)))

   # https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Suggested-packages
   if(parallel && requireNamespace("foreach", quietly=TRUE)) {
      # Avoid dopar warnings
      if(!foreach::getDoParRegistered()) {
	 foreach::registerDoSEQ()
      }
      res <- foreach::`%dopar%`(foreach::foreach(fold=1:nfolds), {
         w <- folds != fold
	 V0 <- NULL
	 if(init) {
	    V0 <- matrix(rnorm(ncol(Y) * ndim), ncol(Y), ndim)
	    s0 <- scca(X[w,], Y[w,], ndim=ndim, verbose=verbose,
	       lambda1=1e-9, lambda2=1e-9, simplify=FALSE, V=V0, ...)
	    V0 <- s0$V
	 }
         scca(X[w,], Y[w,], ndim=ndim, verbose=verbose,
	    lambda1=lambda1, lambda2=lambda2, simplify=TRUE, V=V0, ...)
      })
   } else {
      res <- lapply(1:nfolds, function(fold) {
         w <- folds != fold
	 V0 <- NULL
	 if(verbose) {
	    cat("-> fold: ", fold, "\n")
	 }
	 if(init) {
	    if(verbose) {
	       cat("start init\n")
	    }
	    V0 <- matrix(rnorm(ncol(Y) * ndim), ncol(Y), ndim)
	    s0 <- scca(X[w,], Y[w,], ndim=ndim, verbose=verbose,
	       lambda1=1e-12, lambda2=1e-12, simplify=TRUE, V=V0, ...)
	    V0 <- s0$V 
	    if(verbose) {
	       cat("end init\n")
	    }
	 }
	 scca(X[w,], Y[w,], ndim=ndim, verbose=verbose,
	    lambda1=lambda1, lambda2=lambda2, simplify=FALSE, V=V0, ...)
      })
   }

   # Collect all the folds into one result (same way as glmnet), instead
   # of averaging over k folds.
   #
   # If the model has not converged for one dimension, then
   # we consider the whole model to have not converged, since all dimensions
   # are evaluated sequentially (depend on previous ones)
   for(fold in 1:nfolds) {
      for(i in seq(along=lambda1)) {
	 for(j in seq(along=lambda2)) {
	    x <- res[[fold]][[i]][[j]]
	    w <- folds == fold
	    converged[fold, i, j] <- x$converged
	    if(x$converged) {
	       xpred[w, , i, j] <- X[w,] %*% x$U
	       ypred[w, , i, j] <- Y[w,] %*% x$V
	    } else {
	       xpred[w, , i, j] <- NA
	       ypred[w, , i, j] <- NA
	    }
	 }
      }
   }

   # It's possible that for the same penalties, a model from one set of
   # training folds will converge but another will not.
   for(i in seq(along=lambda1)) {
      for(j in seq(along=lambda2)) {
	 mx <- sapply(res, function(x) colSums(x[[i]][[j]]$U != 0, na.rm=TRUE))
	 my <- sapply(res, function(x) colSums(x[[i]][[j]]$V != 0, na.rm=TRUE))
	 nzx[,i,j] <- rowMeans(rbind(mx))
	 nzy[,i,j] <- rowMeans(rbind(my))
      }
   }

   r <- lapply(1:ndim, function(k) {
      r <- sapply(seq(along=lambda1), function(i) {
	 sapply(seq(along=lambda2), function(j) {
	    m <- suppressWarnings(cor(xpred[,k,i,j], ypred[,k,i,j]))
	    cbind(m)
	 })
      })
      cbind(r) # otherwise abind breaks for 1-length penalties
   })

   r <- abind(r, along=3)
   r <- aperm(r, c(3, 2, 1))

   r.mx <- max(r[opt.dim,,], na.rm=TRUE)
   w <- which(r[opt.dim,,] == r.mx, arr.ind=TRUE)

   res2 <- list(
      ndim=ndim,
      lambda1=lambda1,
      lambda2=lambda2,
      opt.dim=opt.dim,
      best.lambda1=lambda1[w[1]],
      best.lambda2=lambda2[w[2]],
      best.corr=r.mx,
      corr=r,
      nzero.x=nzx,
      nzero.y=nzy,
      nfolds=nfolds,
      converged=converged,
      models=if(return.models) {res} else {NULL},
      folds=if(return.folds) {folds} else {NULL}
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
#'    lambda2=seq(1e-4, 0.5, length=8))
#'
#' par(mfrow=c(1, 2))
#' plot(r, dim=1)
#' plot(r, dim=2)
#' 
#' @importFrom abind abind
#' @importFrom graphics matplot
#' @importFrom graphics legend
#' @importFrom stats rnorm
#'
#' @seealso \code{\link{cv.scca}}
#'
#' @export
plot.cv.scca <- function(x, dim=NULL, type=c("nzero", "penalty"), ...)
{
   if(class(x) != "cv.scca") {
      stop("x is not of class 'cv.scca'")
   }

   type <- match.arg(type)

   if(is.null(dim)) {
      if(type == "nzero") {
	 for(i in 1:x$ndim) {
      	    matplot(x$nzero.x[1,,], x$corr[i,,],
      	       type="b", xlab="# of variables in X with non-zero weight",
      	          ylab="Pearson correlation", log="x",
      	          main=paste0("Canonical dimension ", i))
      	 }
      } else {
	 for(i in 1:x$ndim) {
      	    matplot(x$lambda1, x$corr[i,,],
      	       type="b", xlab="lambda1",
      	          ylab="Pearson correlation", log="x",
      	          main=paste0("Canonical dimension ", i))
	 }
      }
   } else if(is.numeric(dim) && dim > 0 && dim <= x$ndim) {
      if(type == "nzero") {
	 matplot(x$nzero.x[dim,,], x$corr[dim,,],
	    type="b", xlab="# of variables in X with non-zero weight",
	    ylab="Pearson correlation", log="x",
	    main=paste0("Canonical dimension ", dim))
      } else {
	 matplot(x$lambda1, x$corr[dim,,],
	    type="b", xlab="lambda1",
	    ylab="Pearson correlation", log="x",
	    main=paste0("Canonical dimension ", dim))
      }
   } else {
      stop(paste("dim must be a positive integer, <=", x$ndim))
   }
   legend("bottomright", legend=round(x$lambda2, 6),
      title="lambda2", lwd=2, lty=1:length(x$lambda2),
      col=1:length(x$lambda2))
   invisible(x)
}

#' @export plot2d
plot2d <- function(x, ...)
{
   UseMethod("plot2d")
}

#' Plot the results of SCCA cross-validation (Pearson correlation) as a 2D
#' surface. The optimal position is shown as a cross ('+').
#'
#' @param x An object of class "cv.scca"
#'
#' @param dim Integer. Which dimension to plot (all will be plotted by
#' default).
#' 
#' @param plot Logical. Whether to plot (print) the ggplot2 object or just
#' return it.
#' 
#' @details
#'    Plots the cross-validated Pearson correlation, as a 2D surface with the
#'    lambda1 and lambda2 penalties on the axes.
#'
#' @return the ggplot2 object
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
#'    lambda2=seq(1e-4, 0.5, length=8))
#'
#' g <- plot2d(r, dim=1, plot=FALSE)
#' print(g)
#' 
#' @importFrom abind abind
#' @importFrom ggplot2 ggplot aes geom_raster geom_contour scale_fill_viridis_c scale_x_continuous scale_y_continuous geom_point theme_bw
#' @importFrom stats rnorm
#'
#' @method plot2d cv.scca
#'
#' @seealso \code{\link{cv.scca}}
#'
#' @export
plot2d.cv.scca <- function(x, dim=1, plot=FALSE)
{
   dim <- as.integer(dim)
   if(!is.integer(dim) || dim < 1 || dim > x$ndim) {
      stop(paste("dim must be positive integer, <=", x$ndim)) 
   }
   rd <- data.frame(r=as.numeric(x$corr[dim,,]),
      nz.x=as.numeric(x$nzero.x[dim,,]),
      nz.y=as.numeric(x$nzero.y[dim,,]),
      lambda1=rep(x$lambda1, length(x$lambda2)),
      lambda2=rep(x$lambda2, each=length(x$lambda1)))
   g1 <- ggplot(rd, aes(x=lambda1, y=lambda2, fill=r))
   g1 <- g1 + geom_raster()
   g1 <- g1 + geom_contour(aes(z=r), bins=50, colour="white") + theme_bw()
   g1 <- g1 + scale_fill_viridis_c("Correlation")
   g1 <- g1 + scale_x_continuous(expression(lambda[1]))
   g1 <- g1 + scale_y_continuous(expression(lambda[2]))
   g1 <- g1 + geom_point(aes(x=x.best.lambda1, y=x.best.lambda2),
      data=data.frame(x$best.lambda1, x$best.lambda2),
      shape="+", size=8, fill="black")
   if(plot) {
      print(g1)
   }
   invisible(g1)
}

#' Performs a low-rank whitening of a matrix X.
#' 
#' @param X numeric matrix.
#'
#' @param ndim Integer. How many dimensions to use in the whitening.
#'
#' @param keep.svd Logical. Whether to keep the singular value decomposition
#' results.
#' 
#' @param stand character. Which standardisation to use on columns of X.
#'
#' @details A k-rank singular value decomposition of X is done
#' \deqn{X = U_{1:k} D_{1:k} V_{1:k}^T} 
#' The 
#'
#' @examples
#'
#' data(hm3.chr1)
#' X <- scale2(hm3.chr1$bed)
#' Xw <- whiten(X, ndim=10)
#'
#' names(attr(Xw, "svd"))
#'
#' @seealso \code{\link{validate.rank}}
#'
#' @export
#'
whiten <- function(X, ndim=20, keep.svd=TRUE, stand="none")
{
   f <- flashpca(X, ndim=ndim, do_loadings=TRUE, stand=stand, divisor="n1")
   u <- f$vectors
   v <- f$loadings
   d <- sqrt(f$values)
   
   R <- tcrossprod(u, v)
   colnames(R) <- colnames(X)
   rownames(R) <- rownames(X)
   if(keep.svd) {
      attr(R, "svd") <- list(u=u, d=d, v=v)
   }
   R
}

#' Estimates the effective number of dimensions of a matrix X, based on
#' reconstruction error from singular value decomposition.
#' 
#' @param X numeric matrix.
#'
#' @param maxdim Integer. How many dimensions to try.
#'
#' @param test.prop numeric. What proportion of the data to reserve for
#' evaluating the reconstruction error on.
#'
#' @param const numeric. What value to set the test observations to.
#'
#' @details A certain proportion of the entries of X is set to a constant
#' value. We then run singular value decomposition (SVD) on the matrix
#' and for each rank k compute the reconstruction mean squared error between
#' the predicted values \eqn{U_{1:k} D_{1:k} V_{1:k}^T} and the true values at
#' at these missing positions. The rank k with the lowest mean squared error
#' is chosen as the effective rank for whitening.
#'
#' This is basically one step in the imputation method from
#' http://alexhwilliams.info/itsneuronalblog/2018/02/26/crossval
#' 
#' @return a data.frame with number of dimensions and reconstruction error
#' (mse).
#'
#' @seealso \code{\link{whiten}}
#'
#' @examples
#'
#' #######################
#' ## HapMap3 chr1 example
#' data(hm3.chr1)
#' X <- scale2(hm3.chr1$bed)
#'
#' r <- validate.rank(X, maxdim=50)
#' plot(mse ~ dim, data=r)
#'
#' @export
#'
validate.rank <- function(X, maxdim=20, test.prop=0.1, const=0)
{
   X0 <- X1 <- X
   
   mx <- sample(length(X), size=test.prop * length(X))
   
   # training observations
   X0[mx] <- const
   
   # testing observations
   # not strictly necessary to set NAs here, just being pedantic
   X1[-mx] <- NA
   
   n <- nrow(X)
   
   f <- flashpca(X0, ndim=maxdim, do_loadings=TRUE, stand="none", divisor="n1")
   u <- f$vectors
   v <- f$loadings
   d <- sqrt(f$values)
   
   err <- sapply(1:maxdim, function(k) {
      R <- tcrossprod(u[,1:k] %*% diag(d[1:k], k, k), v[,1:k])
      mean((X1[mx] - R[mx])^2)
   })
   
   data.frame(dim=cbind(1:maxdim), mse=err)
}

#' Cross-validation for the FCCA model.
#'
#' @param X An n by p numeric matrix. Note: PLINK data is currently not supported.
#'
#' @param Y An n by k numeric matrix.
#'
#' @param lambda1 Numeric vector. Non-negative L1 penalty on canonical vectors of X.
#'
#' @param lambda2 Numeric vector. Non-negative L1 penalty on canonical vectors of Y.
#'
#' @param gamma1 Numeric vector. Non-negative L2 penalty on X.
#'
#' @param gamma2 Numeric vector. Non-negative L2 penalty on Y.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param nfolds Integer. Number of cross-validation folds.
#'
#' @param folds Integer. The fold identifiers. Overrides the `nfolds'
#' parameter.
#'
#' @param return.models Logical. Whether to return all trained models.
#'
#' @param svd.tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @param verbose Logical. 
#' 
#' @param maxiter Integer. Maximum number of iterations.
#'
#' @param parallal Logical. Whether to use parallelisation. Requires the use
#' of the \code{foreach} package.
#'
#' @return an \code{cv.fcca} object with the following components:
#'
#' \describe{
#'    \item{ndim:} Number of dimensions.
#'    \item{lambda1:} Vector of lambda1 penalties.
#'    \item{lambda2:} Vector of lambda2 penalties.
#'    \item{gamma1:} Vector of gamma1 penalties.
#'    \item{gamma2:} Vector of gamma2 penalties.
#' }
#' 
#' @details
#' Standardisation is done once for X and Y
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
#' r <- cv.fcca(X, Y, ndim=1, nfolds=3,
#'    lambda1=seq(1e-3, 1e-1, length=5),
#'    lambda2=seq(1e-4, 0.5, length=3),
#'    gamma1=10^c(-3, -1),
#'    gamma2=10^c(-3, -1))
#' 
#' @importFrom data.table setDT data.table rbindlist
#' @importFrom foreach `%dopar%` foreach registerDoSEQ
#'
#' @seealso \code{\link{fcca}}
#'
#' @export
#'
cv.fcca <- function(X, Y, lambda1=0, lambda2=0, gamma1=0, gamma2=0, ndim=1,
   nfolds=10, folds=NULL, return.models=FALSE, svd.tol=1e-12,
   verbose=FALSE, maxiter=1000, parallel=FALSE)
{
   if(!is.numeric(X) && !is.null(dim(x))) {
      stop("X must be a numeric matrix")
   }
   if(!is.numeric(Y) && !is.null(dim(Y))) {
      stop("Y must be a numeric matrix")
   }
   if(nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows")
   }

   if(is.null(lambda1) || any(lambda1 < 0)) {
      stop("lambda1 must be a numeric vector >=0") 
   }
   if(is.null(lambda2) || any(lambda2 < 0)) {
      stop("lambda2 must be a numeric vector >=0") 
   }
   if(is.null(gamma1) || any(gamma1 < 0)) {
      stop("gamma1 must be a numeric vector >=0") 
   }
   if(is.null(gamma2) || any(gamma2 < 0)) {
      stop("gamma2 must be a numeric vector >=0") 
   }

   if(!is.numeric(svd.tol) || svd.tol <= .Machine$double.eps
      || length(svd.tol) > 1) {
      stop("svd.tol must be a single number >0")
   }
   
   n <- nrow(X)
   if(is.null(folds)) {
      if(is.null(nfolds) || length(nfolds) > 1 || nfolds < 2 || nfolds >= n) {
	 stop("'nfolds' must be an integer >1 and < sample size")
      }
      folds <- sample(1:nfolds, n, replace=TRUE)
   } else {
      if(length(folds) > n) {
	 stop("Length of 'folds' must match number of rows in X and Y")
      }
      if(any(folds < 1)) {
	 stop("'folds' must be consecutive integers >0")
      }
      if(max(folds) > nfolds) {
	 warning("'folds' will take precedence over 'nfolds' parameter")
      }
      nfolds <- length(unique(folds))
   }

   X <- scale(X)
   Y <- scale(Y)

   # All the SVD + whitening needs to be within the crossval loop..
   run_fold <- function(fold) {
      Xtrn <- X[folds != fold, ]
      Ytrn <- Y[folds != fold, ]
      Xtst <- X[folds == fold, ]
      Ytst <- Y[folds == fold, ]
      n <- nrow(Xtrn)

      s1 <- svd(Xtrn)
      s2 <- svd(Ytrn)
      mx <- s1$d > svd.tol
      my <- s2$d > svd.tol

      d <- lapply(gamma1, function(g1) {
	 d <- lapply(gamma2, function(g2) {
	    if(verbose) {
	       cat("gamma1:", g1, "gamma2: ", g2, "\n")
	    }

	    # Whitened training data
	    Xwtrn <- with(s1,
	       tcrossprod(u[,mx] %*% diag(d[mx] / sqrt(d[mx]^2 + (n - 1) * g1)),
	          v[,mx]))
	    Ywtrn <- with(s2,
	       tcrossprod(u[,my] %*% diag(d[my] / sqrt(d[my]^2 + (n - 1) * g2)),
		  v[,my]))

	    if(verbose) {
	       cat("start scca\n")
	    }
	    st <- system.time({
	       mod <- scca(Xwtrn, Ywtrn, standx="none", standy="none",
	          lambda1=lambda1, lambda2=lambda2, ndim=ndim, verbose=verbose,
		  simplify=FALSE, divisor="none", maxiter=maxiter)
	    })
	    if(verbose) {
	       cat("end scca\n")
	    }

	    Xwtst <- with(s1,
	       tcrossprod(
		  Xtst %*% v[,mx] %*% diag(sqrt((n - 1) / (d[mx]^2 + (n - 1) * g1))),
		  v[,mx]))
	    Ywtst <- with(s2,
	       tcrossprod(
		  Ytst %*% v[,my] %*% diag(sqrt((n - 1) / (d[my]^2 + (n - 1) * g2))),
		  v[,my]))

	    if(verbose) {
	       cat("computing correlations\n")
	    }
	    d <- lapply(seq(along=lambda1), function(i) {
	       d <- lapply(seq(along=lambda2), function(j) {
		  if(verbose) {
		     cat("start")
		  }

		  Za <- crossprod(s1$v[, mx], mod[[i]][[j]]$U)
		  Zb <- crossprod(s2$v[, my], mod[[i]][[j]]$V)
		  mod[[i]][[j]]$a <- with(s1,
		     v[,mx] %*% diag(sqrt((n - 1) / (d[mx]^2 + (n - 1) * g1)))
			%*% Za)
		  mod[[i]][[j]]$b <- with(s2,
		     v[,my] %*% diag(sqrt((n - 1) / (d[my]^2 + (n - 1) * g2)))
			%*% Zb)

		  Px.trn <- Xtrn %*% mod[[i]][[j]]$a
		  Py.trn <- Ytrn %*% mod[[i]][[j]]$b
		  Px.tst <- Xtst %*% mod[[i]][[j]]$a
		  Py.tst <- Ytst %*% mod[[i]][[j]]$b

		  r.trn <- suppressWarnings(diag(cor(Px.trn, Py.trn)))
		  r.tst <- suppressWarnings(diag(cor(Px.tst, Py.tst)))

		  if(verbose) {
		     cat("-end-\n")
		  }
	          dd <- data.table(fold=fold, dim=1:ndim, gamma1=g1, gamma2=g2,
		     lambda1=lambda1[i], lambda2=lambda2[j],
		     r.trn=r.trn, r.tst=r.tst,
		     n.trn=nrow(Px.trn), n.tst=nrow(Px.tst))
		  if(return.models) {
		     list(result=dd, model=mod[[i]][[j]])
		  } else {
		     list(result=dd)
		  }
	       })
	    })
	 })
      })
   }

   res <- if(parallel && requireNamespace("foreach", quietly=TRUE)) {
      if(!foreach::getDoParRegistered()) {
	 foreach::registerDoSEQ()
      }
      foreach::`%dopar%`(foreach::foreach(fold=1:nfolds), { run_fold(fold) })
   } else {
      lapply(1:nfolds, run_fold)
   }

   res.stat <- lapply(seq_along(res), function(i) {
      l <- lapply(seq_along(res[[i]]), function(j) {
	 l <- lapply(seq_along(res[[i]][[j]]), function(k) {
	    l <- lapply(seq_along(res[[i]][[j]][[k]]), function(m) {
	       l <- lapply(seq_along(res[[i]][[j]][[k]][[m]]), function(n) {
		  r <- res[[i]][[j]][[k]][[m]][[n]]$result
		  r[, c("idx.i", "idx.j", "idx.k", "idx.m", "idx.n")
		     := list(i, j, k, m, m)]
		  r
	       })
	       rbindlist(l)
	    })
	    rbindlist(l)
	 })
	 rbindlist(l)
      })
      rbindlist(l)
   })
   res.stat <- rbindlist(res.stat)

   # Average and stderr of correlation, across the dimensions, for all grid
   # points
   res.stat.agg <- res.stat[,
      list(r.trn.mean=mean(r.trn), r.trn.se=sqrt(var(r.trn) / .N),
      r.tst.mean=mean(r.tst), r.tst.se=sqrt(var(r.tst) / .N)),
      by=.(dim, gamma1, gamma2, lambda1, lambda2)]

   # Sum of squared correlations across the dimensions, for each model
   # # on the penalty grid. I.e., the best overall model taking all
   # dimensions into account.
   res.stat.agg.avg <- res.stat.agg[,
      list(avg.sq.cor=mean(r.tst.mean^2, na.rm=TRUE)),
      by=.(gamma1, gamma2, lambda1, lambda2)]
   res.stat.agg.avg.best <- res.stat.agg.avg[which.max(avg.sq.cor), ]
   
   # results are for all folds, need to average over them
   obj <- list(ndim=ndim, lambda1=lambda1, lambda2=lambda2,
      gamma1=gamma1, gamma2=gamma2,
      nfolds=nfolds, folds=folds,
      result.raw=res.stat, 
      result.agg=res.stat.agg,
      result=res.stat.agg.avg,
      result.best=res.stat.agg.avg.best)
   if(return.models) {
      obj$models <- res
   }
   class(obj) <- "cv.fcca"
   obj
}

#' Fit the FCCA model.
#'
#' @param X An n by p numeric matrix. Note: PLINK data is currently not supported.
#'
#' @param Y An n by k numeric matrix.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param lambda1 Numeric vector. Non-negative L1 penalty on canonical vectors of X.
#'
#' @param lambda2 Numeric vector. Non-negative L1 penalty on canonical vectors of Y.
#'
#' @param gamma1 Numeric vector. Non-negative L2 penalty on X.
#'
#' @param gamma2 Numeric vector. Non-negative L2 penalty on Y.
#'
#' @param V Numeric matrix. Initial estimate of the left singular
#' vector of  X'Y.
#'
#' @param standx Character. One of "sd" (empricial standard deviation) or "none".
#'
#' @param standy Character. One of "sd" (empricial standard deviation) or "none".
#'
#' @param svd.tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @param maxiter Integer. Maximum number of iterations.
#'
#' @param verbose Logical. Whether to print more information.
#'
#' @return an \code{fcca} object with the following components:
#'
#' \describe{
#'    \item{ndim:} Number of dimensions.
#'    \item{lambda1:} Vector of lambda1 penalties.
#'    \item{lambda2:} Vector of lambda2 penalties.
#'    \item{gamma1:} Vector of gamma1 penalties.
#'    \item{gamma2:} Vector of gamma2 penalties.
#'    \item{U:} Matrix of left singular vectors
#'    \item{a:} Matrix of canonical vectors for X.
#'    \item{V:} Matrix of right singular vectors
#'    \item{b:} Matrix of canonical vectors for Y.
#'    \item{Px:} Matrix of canonical coordinates, X * a
#'    \item{Py:} Matrix of canonical coordinates, Y * b
#'    \item{r:} Canonical correlations, diag(cor(Px, Py))
#'    \item{converged:} whether the optimisation has converged or not
#' }
#' 
#' @details
#'
#' Standardisation is done once for X and Y
#'
#' @export
#'
fcca <- function(X, Y, ndim=1, lambda1=0, lambda2=0, gamma1=0, gamma2=0,
   V=NULL, standx=c("sd", "none"), standy=c("sd", "none"),
   svd.tol=1e-12, maxiter=1000, verbose=FALSE)
{
   if(mode(X) != "numeric" || any(dim(X) == 0)) {
      stop("X must be a numeric matrix")
   }
   if(mode(Y) != "numeric" || any(dim(Y) == 0)) {
      stop("Y must be a numeric matrix")
   }
   if(nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows")
   }
   if(!is.numeric(svd.tol) || svd.tol <= .Machine$double.eps
      || length(svd.tol) > 1) {
      stop("svd.tol must be a single number >0")
   }
   if(!is.numeric(lambda1) || length(lambda1) > 1 || lambda1 < 0) {
      stop("lambda1 must be a single number >= 0")
   }
   if(!is.numeric(lambda2) || length(lambda1) > 2 || lambda1 < 0) {
      stop("lambda2 must be a single number >= 0")
   }
   if(!is.numeric(gamma1) || length(gamma1) > 1 || gamma1 < 0) {
      stop("gamma1 must be a single number >= 0")
   }
   if(!is.numeric(gamma2) || length(gamma1) > 2 || gamma1 < 0) {
      stop("gamma2 must be a single number >= 0")
   }
   if(!is.numeric(ndim) || length(ndim) > 1 || ndim <= 0) {
      stop("ndim must be a single number > 0")
   }
   standx <- match.arg(standx)
   standy <- match.arg(standy)

   if(verbose) {
      cat("start svd\n") 
   }

   n <- nrow(X)

   if(standx == "sd") {
      X <- scale(X)
   }
   if(standy == "sd") {
      Y <- scale(Y)
   }
   s1 <- svd(X)
   s2 <- svd(Y)
   mx <- s1$d > svd.tol
   my <- s2$d > svd.tol

   if(verbose) {
      cat("end svd\n")
   }

   Xw <- with(s1,
      tcrossprod(u[,mx] %*% diag(d[mx] / sqrt(d[mx]^2 + (n - 1) * gamma1)),
	 v[,mx]))
   Yw <-  with(s2,
      tcrossprod(u[,my] %*% diag(d[my] / sqrt(d[my]^2 + (n - 1) * gamma2)),
	 v[,my]))

   if(is.null(V)) {
      r <- scca(Xw, Yw, ndim=ndim, lambda1=lambda1, lambda2=lambda2,
	 standx="none", standy="none", divisor="none",
	 verbose=verbose, maxiter=maxiter)
   } else {
      r <- scca(Xw, Yw, ndim=ndim, lambda1=lambda1, lambda2=lambda2,
	 standx="none", standy="none", divisor="none", V=V,
	 verbose=verbose, maxiter=maxiter)
   }

   if(verbose) {
      cat("end scca\n")
   }
   Za <- crossprod(s1$v[, mx], r$U)
   Zb <- crossprod(s2$v[, my], r$V)
   a <- with(s1,
      v[,mx] %*% diag(sqrt((n - 1) / (d[mx]^2 + (n - 1) * gamma1))) %*% Za)
   b <- with(s2,
      v[,my] %*% diag(sqrt((n - 1) / (d[my]^2 + (n - 1) * gamma2))) %*% Zb)
      
   rownames(a) <- colnames(X)
   rownames(b) <- colnames(Y)
   Px <- X %*% a
   Py <- Y %*% b
   colnames(Px) <- paste0("Px", 1:ncol(Px))
   colnames(Py) <- paste0("Py", 1:ncol(Py))
   rp <- suppressWarnings(diag(cor(Px, Py)))

   res <- list(ndim=r$ndim, U=r$U, V=r$V, d=r$d,
      a=a, b=b, Px=Px, Py=Py, r=rp,
      converged=r$converged)
   class(res) <- "fcca"
   res
}

#' Prints a cv.fcca object
#'
#' @param x A cv.fcca object to be printed
#' @param ... Ignored
#' @export 
print.cv.fcca <- function(x, ...)
{
   cat(paste0(
      "cv.fcca object; ", x$nfolds,
      "-fold cross-validation\n"))
   invisible(x)
}

#' Prints an fcca object
#'
#' @param x An fcca object to be printed
#' @param ... Ignored
#' @export 
print.fcca <- function(x, ...)
{
   cat("fcca object; ndim=", length(x$d), "\n")
   if(!x$converged) {
      cat("warning: model has not converged\n")
   }
   invisible(x)
}

#' Optimise FCCA via grid search and/or Bayesian optimisation.
#'
#' @param X An n by p numeric matrix. Note: PLINK data is currently not supported.
#'
#' @param Y An n by k numeric matrix.
#'
#' @param ndim Integer. Positive number of canonical vectors to compute.
#'
#' @param nfolds Integer. Number of cross-validation folds.
#'
#' @param folds Integer. The fold identifiers. Overrides the `nfolds'
#' parameter.
#'
#' @param lambda1.grid Numeric vector. Non-negative L1 penalty on canonical
#' vectors of X for use in the grid search.
#'
#' @param lambda2.grid Numeric vector. Non-negative L1 penalty on canonical vectors of Y
#' for use in the grid search.
#'
#' @param gamma1.grid Numeric vector. Non-negative L2 penalty on X for
#' use in the grid search.
#'
#' @param gamma2.grid Numeric vector. Non-negative L2 penalty on Y for
#' use in the grid search.
#'
#' @param lambda1.bopt Numeric vector. Non-negative L1 penalty on canonical
#' vectors of X for use in the Bayesian optimisation.
#'
#' @param lambda2.bopt Numeric vector. Non-negative L1 penalty on canonical vectors of Y
#' for use in the Bayesian optimisation.
#'
#' @param gamma1.bopt Numeric vector. Non-negative L2 penalty on X for
#' use in the Bayesian optimisation.
#'
#' @param gamma2.bopt Numeric vector. Non-negative L2 penalty on Y for
#' use in the Bayesian optimisation.
#' 
#' @param method Character. Either "bapt" (Bayesian optimisation) or "grid"
#' (grid search).
#'
#' @param parallal Logical. Whether to use parallelisation. Requires the use
#' of the \code{foreach} package.
#' 
#' @param final.model Logical. Whether to fit and return the final model trained on
#' all the data using the optimal hyperparameters.
#' 
#' @param final.model.cv Logical. Whether to fit a model in cross-validation 
#' using the optimal hyperparameters, and return it.
#'
#' @param final.model.reorder Logical. Whether to reorder the canonical
#' coordinates based on the cross-validation correlations (in decreasing
#' magnitude).
#'
#' @param maxiter Integer. Maximum number of iterations for FCCA fitting.
#'
#' @param svd.tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @param verbose Logical. 
#'
#' @return an \code{optim.cv.fcca} object with the following components:
#'
#' \describe{
#'    \item{ndim:} Number of dimensions.
#'    \item{folds:} The vector of folds.
#'    \item{nfolds:} The number of folds.
#'    \item{grid.path:} The results for the grid search.
#'    \item{bopt:} The mlrMBO result (if using Bayesian optimisation).
#'    \item{bopt.path:} The Bayesian optimisation results (if using Bayesian
#' optimisation).
#'    \item{opt.param:} The optimal hyperparameters .
#'    \item{final.model:} An FCCA model trained on all the data using the
#'	 optimal hyperparameters.
#'    \item{final.model.cv:} An FCCA model trained in cross-validation using
#'	 the optimal hyperparameters.
#'    \item{final.model.cv.Px:} The canonical coordinates for X from
#'	 the cross-validated FCCA model.
#'    \item{final.model.cv.Py:} The canonical coordinates for Y from
#'	 the cross-validated FCCA model.
#'    \item{final.model.reordered:} Whether the canonical coordinates for the 
#'	 final model were reordered by the cross-validated canonical
#'	 correlations.
#' }
#' 
#' @examples
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
#' r <- optim.cv.fcca(X, Y, ndim=1, nfolds=3,
#'    lambda1.grid=seq(1e-3, 1e-1, length=5),
#'    lambda2.grid=seq(1e-4, 0.5, length=3),
#'    gamma1.grid=10^c(-3, -1),
#'    gamma2.grid=10^c(-3, -1), method="grid",
#'    final.model.cv=TRUE)
#'
#' # The optimal hyperparameters
#' r$opt.param
#'
#' # The cross-validated canonical correlations
#' diag(cor(r$final.model.cv.Px, r$final.model.cv.Py))
#'
#' @importFrom mlrMBO makeMBOControl
#' @importFrom data.table copy setcolorder setnames as.data.table
#' @importFrom foreach foreach `%dopar%` registerDoSEQ
#' @export
optim.cv.fcca <- function(X, Y, ndim=1, nfolds=5, folds=NULL,
   lambda1.grid=c(0, seq(1e-4, 0.002, length=4)),
   lambda2.grid=c(0, seq(1e-4, 0.002, length=4)),
   gamma1.grid=10^(-2:0), gamma2.grid=10^(-2:0),
   lambda1.bopt=c(0, 1), lambda2.bopt=c(0, 1),
   gamma1.bopt=10^c(-3, 4), gamma2.bopt=10^c(-3, 4),
   method=c("bopt", "grid"), parallel=FALSE,
   final.model=TRUE, final.model.cv=FALSE, final.model.reorder=FALSE,
   maxiter=1e3, svd.tol=1e-12, verbose=FALSE)
{
   method <- match.arg(method)
   if(final.model.reorder && !final.model.cv) {
      stop("final.model.cv must be TRUE if final.model.reorder is TRUE")
   }

   if(method == "bopt"
      && (!requireNamespace("mlrMBO", quietly=TRUE)
	 || !requireNamespace("smoof", quietly=TRUE))) {
      stop("Bayesian optimisation requires the packages: mlrMBO, smoof")
   }

   if(!is.null(folds)) {
      nfolds <- max(folds)
   } else {
      folds <- sample(nfolds, nrow(X), replace=TRUE)
   }

   if(method == "bopt") {
      if(is.null(lambda1.bopt)) {
         stop("lambda1.bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(is.null(lambda2.bopt)) {
         stop("lambda2.bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(is.null(gamma1.bopt)) {
         stop("gamma1.bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(is.null(gamma2.bopt)) {{
         stop("gamma2.bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(all(sapply(
	 list(lambda1.bopt, lambda2.bopt, gamma1.bopt, gamma2.bopt),
	    length) == 1))
	    stop("Can't do Bayesian optimisation if lambda1.bopt, ",
	       "lambda2.bopt, gamma1.bopt, gamma2.bopt are all constants")
      }
   }

   if(is.null(lambda1.grid) || any(lambda1.grid < 0)) {
      stop("lambda1.grid must be a numeric vector of values >= 0")
   }
   if(is.null(lambda2.grid) || any(lambda2.grid < 0)) {
      stop("lambda2.grid must be a numeric vector of values >=0")
   }
   if(is.null(gamma1.grid) || any(gamma1.grid < 0)) {
      stop("gamma1.grid must be a numeric vector of values >=0")
   }
   if(is.null(gamma2.grid) || any(gamma2.grid < 0)) {
      stop("gamma2.grid must be a numeric vector of values >=0")
   }
   
   # Setup initial 'warmup' results for mlrMBO
   des.fcca.cv <- cv.fcca(X, Y, ndim=ndim,
      lambda1=lambda1.grid, lambda2=lambda2.grid,
      gamma1=gamma1.grid, gamma2=gamma2.grid,
      folds=folds, parallel=parallel,
      svd.tol=svd.tol, verbose=verbose)
   
   des.fcca <- copy(des.fcca.cv$result)
   setcolorder(des.fcca, c(3, 4, 1, 2))
   des.fcca[, log10_gamma1 := log10(gamma1)]
   des.fcca[, log10_gamma2 := log10(gamma2)]

   if(method == "bopt") {
      des.fcca.d <- as.data.frame(des.fcca[, 
	 .(lambda1, lambda2, gamma1, gamma2, avg.sq.cor)])
      colnames(des.fcca.d)[5] <- "y"

      # Super ugly hack to get around the fact that mlrMBO cannot handle
      # constant hyperpameters but needs them as extra arguments
      par.set.l <- vector("list", 4)
      more.args <- list()

      if(length(lambda1.bopt) > 1) {
	 par.set.l[[1]] <- ParamHelpers::makeNumericParam("lambda1",
	    lower=lambda1.bopt[1], upper=lambda1.bopt[2])
      } else {
	 more.args <- c(more.args, c("lambda1"=lambda1.bopt))
      }

      if(length(lambda2.bopt) > 1) {
	 par.set.l[[2]] <- ParamHelpers::makeNumericParam("lambda2",
	    lower=lambda2.bopt[1], upper=lambda2.bopt[2])
      } else {
	 more.args <- c(more.args, c("lambda2"=lambda2.bopt))
      }
      
      if(length(gamma1.bopt) > 1) {
	 par.set.l[[3]] <- ParamHelpers::makeNumericParam("gamma1",
	    lower=log10(gamma1.bopt[1]), upper=log10(gamma1.bopt[2]),
	    trafo=function(x) 10^x)
      } else {
	 more.args <- c(more.args, c("gamma1"=gamma1.bopt))
      }

      if(length(gamma2.bopt) > 1) {
	 par.set.l[[4]] <- ParamHelpers::makeNumericParam("gamma2",
	    lower=log10(gamma2.bopt[1]), upper=log10(gamma2.bopt[2]),
	    trafo=function(x) 10^x)
      } else {
	 more.args <- c(more.args, c("gamma2"=gamma2.bopt))
      }

      par.set.l <- par.set.l[!sapply(par.set.l, is.null)]
      par.set <- do.call(ParamHelpers::makeParamSet, par.set.l)

      obj.fun.fcca <- smoof::makeSingleObjectiveFunction(
         name="fcca",
	 fn=function(x, ...) {
	    more.args <- list(...)
	    if("lambda1" %in% names(x)) {
	       lambda1 <- x["lambda1"]
	    } else if("lambda1" %in% names(more.args)) {
	       lambda1 <- more.args$lambda1
	    }
	    if("lambda2" %in% names(x)) {
	       lambda2 <- x["lambda2"]
	    } else if("lambda2" %in% names(more.args)) {
	       lambda2 <- more.args$lambda2
	    }
	    if("gamma1" %in% names(x)) {
	       gamma1 <- x["gamma1"]
	    } else if("gamma1" %in% names(more.args)) {
	       gamma1 <- more.args$gamma1
	    }
	    if("gamma2" %in% names(x)) {
	       gamma2 <- x["gamma2"]
	    } else if("gamma2" %in% names(more.args)) {
	       gamma2 <- more.args$gamma2
	    }

	    r <- cv.fcca(X, Y, ndim=ndim,
	       lambda1=lambda1, lambda2=lambda2, gamma1=gamma1, gamma2=gamma2,
	       folds=folds, parallel=parallel, svd.tol=svd.tol,
	       maxiter=maxiter, verbose=FALSE)
	    m <- r$result$avg.sq.cor
	    if(verbose) {
	       cat("x:", x, "m:", m, "\n")
	    }
	    #ifelse(is.nan(m) || !is.finite(m), runif(1, 0, 1e-2), m)
	    ifelse(is.nan(m), 0, m)
	 },
	 par.set=par.set,
	 minimize=FALSE,
	 noisy=TRUE
      )

      ctrl <- mlrMBO::makeMBOControl()
      ctrl <- mlrMBO::setMBOControlTermination(ctrl, iters=50L)
      ctrl <- mlrMBO::setMBOControlInfill(ctrl,
	 crit=mlrMBO::makeMBOInfillCritEI())
      surr.km.fcca <- mlrMBO::makeMBOLearner(control=ctrl, fun=obj.fun.fcca,
	 config=list(show.learner.output=verbose,
	 on.learner.warning=ifelse(verbose, "warn", "quiet")))

      des.fcca.d <- subset(des.fcca.d, !is.nan(y))
      des.fcca.d <- des.fcca.d[, c(names(par.set$pars), "y")]

      run.fcca <- mlrMBO::mbo(obj.fun.fcca, design=des.fcca.d,
         learner=surr.km.fcca, control=ctrl, show.info=verbose,
	 more.args=more.args)
      run.fcca.path <- as.data.table(run.fcca$opt.path)

      # Note that gamma1, gamma2 are retured on log10 scale
      opt.param <- c(
	 lambda1=run.fcca$x$lambda1, lambda2=run.fcca$x$lambda2,
	    gamma1=10^(run.fcca$x$gamma1), gamma2=10^(run.fcca$x$gamma2))
      if(length(more.args) > 0) {
	 opt.param <- c(opt.param, unlist(more.args))
	 opt.param <- opt.param[c("lambda1", "lambda2", "gamma1", "gamma2")]
      }
   } else {
      dmax <- des.fcca[which.max(avg.sq.cor), ]
      opt.param <- c(lambda1=dmax$lambda1, lambda2=dmax$lambda2,
	 gamma1=10^(dmax$log10_gamma1), gamma2=10^(dmax$log10_gamma2))
   }

   mod.fcca <- mod.fcca.cv <- NULL
   Px.cv <- Py.cv <- NULL
   reordered <- FALSE

   # A final model on all the data
   if(final.model) {
      mod.fcca <- fcca(X, Y, ndim=ndim,
	 lambda1=opt.param["lambda1"], lambda2=opt.param["lambda2"],
	 gamma1=opt.param["gamma1"], gamma2=opt.param["gamma2"],
	 verbose=verbose, svd.tol=svd.tol, maxiter=maxiter)

      # It's difficult to extract the cross-validation results
      # for the final model from the optimisation process so we do it again
      # here for the optimal hyperparameters
      if(final.model.cv) {
	 mod.fcca.cv <- cv.fcca(X, Y, ndim=ndim, folds=folds,
	    lambda1=opt.param["lambda1"], lambda2=opt.param["lambda2"],
	    gamma1=opt.param["gamma1"], gamma2=opt.param["gamma2"],
	    verbose=verbose, return.models=TRUE, svd.tol=svd.tol,
	    maxiter=maxiter)

	 Px.cv <- Py.cv <- matrix(0, nrow(X), ndim)
	 colnames(Px.cv) <- paste0("Px", 1:ndim)
	 colnames(Py.cv) <- paste0("Py", 1:ndim)
	 for(fold in 1:nfolds) {
	    mod <- mod.fcca.cv$models[[fold]][[1]][[1]][[1]][[1]]$model
	    Px.cv[folds == fold,] <- X[folds == fold,] %*% mod$a
	    Py.cv[folds == fold,] <- Y[folds == fold,] %*% mod$b
	 }

	 # Re-order the final model's canonical vectors based on the
	 # cross-validated canonical correlations (decreasing magnitude)
	 if(final.model.reorder && ndim > 1) {
	    ord <- mod.fcca.cv$result.agg[,
	       order(r.tst.mean, decreasing=TRUE)]
	    if(any(ord != seq_along(ord))) {
	       reordered <- TRUE
	    }
	    mod.fcca$U <- mod.fcca$U[, ord]   
	    mod.fcca$V <- mod.fcca$V[, ord]   
	    mod.fcca$a <- mod.fcca$a[, ord]   
	    mod.fcca$b <- mod.fcca$b[, ord]   
	    mod.fcca$Px <- mod.fcca$Px[, ord]   
	    mod.fcca$Py <- mod.fcca$Py[, ord]   
	    mod.fcca$r <- mod.fcca$r[ord]   
	 }
      }
   }

   if(method == "grid") {
      res <- list(
	 ndim=ndim, folds=folds, nfolds=nfolds, grid.path=des.fcca,
	 opt.param=opt.param, final.model=mod.fcca, 
	 final.model.cv=mod.fcca.cv, final.model.cv.Px=Px.cv,
	 final.model.cv.Py=Py.cv, final.model.reordered=reordered)
   } else {
      res <- list(
	 ndim=ndim, folds=folds, nfolds=nfolds,
	 grid.path=des.fcca, bopt=run.fcca,
	 bopt.path=run.fcca.path, opt.param=opt.param,
	 final.model=mod.fcca, final.model.cv=mod.fcca.cv,
	 final.model.cv.Px=Px.cv, final.model.cv.Py=Py.cv,
	 final.model.reordered=reordered)
   }
   class(res) <- "optim.cv.fcca"
   res
}

#' Prints an optim.cv.fcca object
#'
#' @param x An optim.cv.fcca object to be printed
#' @param ... Ignored
#' @export 
print.optim.cv.fcca <- function(x, ...)
{
   cat(paste0(
      "optim.cv.fcca object; ", x$nfolds,
      "-fold cross-validation (re-ordered: ",
      x$final.model.reordered, ")\n"))
   invisible(x)
}

