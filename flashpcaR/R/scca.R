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
	 x <- res[[1]][[1]]
	 class(x) <- "scca"
	 x
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
   cat("scca object; ndim=", length(x$d), "\n")
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

#' @export
cv.scca.ridge <- function(X, Y, lambda1, lambda2,
   gamma1=0, gamma2=0, nfolds=10, svd.tol=1e-12, ...)
{
   if(!is.numeric(X) && !is.null(dim(x))) {
      stop("X must be a numeric matrix")
   }
   if(!is.numeric(Y) && !is.null(dim(Y))) {
      stop("Y must be a numeric matrix")
   }
   
   n <- nrow(X)
   folds <- sample(1:nfolds, n, replace=TRUE)

   X <- scale(X)
   Y <- scale(Y)

   # All the SVD + whitening needs to be within the crossval loop..
   d <- foreach(fold=1:nfolds) %dopar% {
      Xtrn <- X[folds != fold, ]
      Ytrn <- Y[folds != fold, ]
      n <- nrow(Xtrn)
      s1 <- svd(Xtrn)
      s2 <- svd(Ytrn)
      mx <- s1$d > svd.tol
      my <- s2$d > svd.tol

      d <- lapply(gamma1, function(g1) {
	 d <- lapply(gamma2, function(g2) {
	    cat("gamma1:", g1, "gamma2: ", g2, "\n")

	    # Whitened training data
	    Xwtrn <- sqrt(n - 1) * with(s1,
	       tcrossprod(u[,mx] %*% diag(d[mx] / sqrt(d[mx]^2 + (n - 1) * g1)),
	          v[,mx]))
	    Ywtrn <- sqrt(n - 1) * with(s2,
	       tcrossprod(u[,my] %*% diag(d[my] / sqrt(d[my]^2 + (n - 1) * g2)),
	          v[,my]))

	    cat("start scca\n")
	    print(system.time({
	       res <- scca(Xwtrn, Ywtrn, standx="none", standy="none",
	          lambda1=lambda1, lambda2=lambda2, ndim=1)
	    }))
	    cat("end scca\n")

	    cat("foo\n")
	    # (X'X / (n-1))^{-1/2}
	    sx.invsqrt <- with(s1,
	       sqrt(n - 1) * tcrossprod(
	          v[,mx] %*% diag(1 / sqrt(d[mx]^2 + (n - 1) * g1)), v[,mx]))

	    cat("bar\n")
	    # (Y'Y / (n-1))^{-1/2}
	    sy.invsqrt <- with(s2,
	       sqrt(n - 1) * tcrossprod(
	          v[,my] %*% diag(1 / sqrt(d[my]^2 + (n - 1) * g2)), v[,my]))

	    d <- lapply(seq(along=lambda1), function(i) {
	       d <- lapply(seq(along=lambda2), function(j) {
	          cat(fold, i, j, "\n")
	          a <- sx.invsqrt %*% res[[i]][[j]]$U
	          b <- sy.invsqrt %*% res[[i]][[j]]$V
	          Px <- X[folds == fold,] %*% a
	          Py <- Y[folds == fold,] %*% b
	          r <- suppressWarnings(diag(cor(Px, Py)))
	          data.frame(gamma1=g1, gamma2=g2,
		     lambda1=lambda1[i], lambda2=lambda2[j], r=r,
		     n=nrow(Px))
	       })
	       do.call(rbind, d)
	    })
	    do.call(rbind, d)
	 })
	 do.call(rbind, d)
      })
      do.call(rbind, d)
   }
   res <- do.call(rbind, d)
   res
}

#' @export
scca.ridge <- function(X, Y, gamma1, gamma2, svd.tol=1e-12, ...)
{
   n <- nrow(X)
   s1 <- svd(X)
   s2 <- svd(Y)
   mx <- s1$d > svd.tol
   my <- s2$d > svd.tol

   Xw <- sqrt(n - 1) * with(s1,
      tcrossprod(u[,mx] %*% diag(d[mx] / sqrt(d[mx]^2 + (n - 1) * gamma1)),
	 v[,mx]))
   Yw <- sqrt(n - 1) * with(s2,
      tcrossprod(u[,my] %*% diag(d[my] / sqrt(d[my]^2 + (n - 1) * gamma2)),
	 v[,my]))

   # (X'X / (n-1))^{-1/2}
   sx.invsqrt <- with(s1,
      sqrt(n - 1) * tcrossprod(
         v[,mx] %*% diag(1 / sqrt(d[mx]^2 + (n - 1) * gamma1)), v[,mx]))

   # (Y'Y / (n-1))^{-1/2}
   sy.invsqrt <- with(s2,
      sqrt(n - 1) * tcrossprod(
         v[,my] %*% diag(1 / sqrt(d[my]^2 + (n - 1) * gamma2)), v[,my]))

   r <- scca(Xw, Yw, standx="none", standy="none", ...)
   a <- sx.invsqrt %*% r$U
   b <- sy.invsqrt %*% r$V
   Px <- X %*% a
   Py <- Y %*% b
   colnames(Px) <- paste0("Px", 1:ncol(Px))
   colnames(Py) <- paste0("Py", 1:ncol(Py))
   rp <- diag(cor(Px, Py))

   res <- list(U=r$U, V=r$V, a=a, b=b, Px=Px, Py=Py, r=rp)
   class(res) <- "scca.ridge"
   res
}

