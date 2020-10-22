
# Import so that testthat loads data.table
#' @import data.table
NULL

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
#' @param check_sign Logical. Whether to check and correct the sign
#' of the singular vectors to prevent sign flipping in cross-validation.
#'
#' @param return_models Logical. Whether to return all trained models.
#'
#' @param svd_tol Numeric. Tolerance under which to truncate singular values of X and Y.
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
#' @importFrom data.table data.table copy setcolorder setnames as.data.table rbindlist
#' @importFrom foreach `%dopar%` foreach registerDoSEQ
#'
#' @seealso \code{\link{fcca}}
#'
#' @export
#'
cv.fcca <- function(X, Y, lambda1=0, lambda2=0, gamma1=0, gamma2=0, ndim=1,
   nfolds=10, folds=NULL, check_sign=TRUE, return_models=FALSE,
   svd_tol=1e-12, verbose=FALSE, maxiter=1000, parallel=FALSE)
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

   if(!is.numeric(svd_tol) || svd_tol <= .Machine$double.eps
      || length(svd_tol) > 1) {
      stop("svd_tol must be a single number >0")
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
      mx <- s1$d > svd_tol
      my <- s2$d > svd_tol

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
		  mod[[i]][[j]]$A <- with(s1,
		     v[,mx] %*% diag(sqrt((n - 1) / (d[mx]^2 + (n - 1) * g1)))
			%*% Za)
		  mod[[i]][[j]]$B <- with(s2,
		     v[,my] %*% diag(sqrt((n - 1) / (d[my]^2 + (n - 1) * g2)))
			%*% Zb)

		  if(check_sign) {
		     sc <- check.eig.sign(mod[[i]][[j]]$A, mod[[i]][[j]]$B, Xtrn, Ytrn)
		     mod[[i]][[j]]$A <- sc$A
		     mod[[i]][[j]]$B <- sc$B
		  }
		  Px.trn <- Xtrn %*% mod[[i]][[j]]$A 
		  Py.trn <- Ytrn %*% mod[[i]][[j]]$B
		  Px.tst <- Xtst %*% mod[[i]][[j]]$A
		  Py.tst <- Ytst %*% mod[[i]][[j]]$B

		  r.trn <- suppressWarnings(diag(cor(Px.trn, Py.trn)))
		  r.tst <- suppressWarnings(diag(cor(Px.tst, Py.tst)))

		  if(verbose) {
		     cat("-end-\n")
		  }
	          dd <- data.table(fold=fold, dim=1:ndim, gamma1=g1, gamma2=g2,
		     lambda1=lambda1[i], lambda2=lambda2[j],
		     r.trn=r.trn, r.tst=r.tst,
		     n.trn=nrow(Px.trn), n.tst=nrow(Px.tst))
		  if(return_models) {
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
      result_raw=res.stat, 
      result_agg=res.stat.agg,
      result=res.stat.agg.avg,
      result_best=res.stat.agg.avg.best)
   if(return_models) {
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
#' @param svd_tol Numeric. Tolerance under which to truncate singular values of X and Y.
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
   check_sign=FALSE, svd_tol=1e-12, maxiter=1000, verbose=FALSE)
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
   if(!is.numeric(svd_tol) || svd_tol <= .Machine$double.eps
      || length(svd_tol) > 1) {
      stop("svd_tol must be a single number >0")
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
   mx <- s1$d > svd_tol
   my <- s2$d > svd_tol

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
   A <- with(s1,
      v[,mx] %*% diag(sqrt((n - 1) / (d[mx]^2 + (n - 1) * gamma1))) %*% Za)
   B <- with(s2,
      v[,my] %*% diag(sqrt((n - 1) / (d[my]^2 + (n - 1) * gamma2))) %*% Zb)
      
   rownames(A) <- colnames(X)
   rownames(B) <- colnames(Y)

   if(check_sign) {
      sc <- check.eig.sign(A, B, X, Y)
      A <- sc$A
      B <- sc$B
   }

   Px <- X %*% A
   Py <- Y %*% B
   colnames(Px) <- paste0("Px", 1:ncol(Px))
   colnames(Py) <- paste0("Py", 1:ncol(Py))
   rp <- suppressWarnings(diag(cor(Px, Py)))

   res <- list(ndim=r$ndim, U=r$U, V=r$V, d=r$d,
      A=A, B=B, Px=Px, Py=Py, r=rp,
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
#' @param lambda1_grid Numeric vector. Non-negative L1 penalty on canonical
#' vectors of X for use in the grid search.
#'
#' @param lambda2_grid Numeric vector. Non-negative L1 penalty on canonical vectors of Y
#' for use in the grid search.
#'
#' @param gamma1_grid Numeric vector. Non-negative L2 penalty on X for
#' use in the grid search.
#'
#' @param gamma2_grid Numeric vector. Non-negative L2 penalty on Y for
#' use in the grid search.
#'
#' @param lambda1_bopt Numeric vector. Non-negative L1 penalty on canonical
#' vectors of X for use in the Bayesian optimisation.
#'
#' @param lambda2_bopt Numeric vector. Non-negative L1 penalty on canonical vectors of Y
#' for use in the Bayesian optimisation.
#'
#' @param gamma1_bopt Numeric vector. Non-negative L2 penalty on X for
#' use in the Bayesian optimisation.
#'
#' @param gamma2_bopt Numeric vector. Non-negative L2 penalty on Y for
#' use in the Bayesian optimisation.
#' 
#' @param method Character. Either "bapt" (Bayesian optimisation) or "grid"
#' (grid search).
#'
#' @param parallal Logical. Whether to use parallelisation. Requires the use
#' of the \code{foreach} package.
#' 
#' @param final_model Logical. Whether to fit and return the final model trained on
#' all the data using the optimal hyperparameters.
#' 
#' @param final_model_cv Logical. Whether to fit a model in cross-validation 
#' using the optimal hyperparameters, and return it.
#'
#' @param final_model_reorder Logical. Whether to reorder the canonical
#' coordinates based on the cross-validation correlations (in decreasing
#' magnitude).
#'
#' @param check_sign Logical. Whether to check and correct the sign
#' of the singular vectors to prevent sign flipping in cross-validation.
#'
#' @param maxiter Integer. Maximum number of iterations for FCCA fitting.
#'
#' @param svd_tol Numeric. Tolerance under which to truncate singular values of X and Y.
#'
#' @param verbose Logical. 
#'
#' @return an \code{optim.cv.fcca} object with the following components:
#'
#' @details
#' Performs optimisation of the FCCA hyperparmeters {lambda1, lambda2, gamma1,
#' gamma2}, initially via grid search and optionally further tuned
#' by Bayesian optimisation.
#' 
#' The objective function is the cross-validated average r^2 between the 
#' predictions for X and Y (averaged over all k dimensions).
#'
#' \describe{
#'    \item{ndim:} Number of dimensions.
#'    \item{folds:} The vector of folds.
#'    \item{nfolds:} The number of folds.
#'    \item{grid_path:} The results for the grid search.
#'    \item{bopt:} The mlrMBO result (if using Bayesian optimisation).
#'    \item{bopt_path:} The Bayesian optimisation results (if using Bayesian
#' optimisation).
#'    \item{opt_param:} The optimal hyperparameters .
#'    \item{final_model:} An FCCA model trained on all the data using the
#'	 optimal hyperparameters.
#'    \item{final_model_cv:} An FCCA model trained in cross-validation using
#'	 the optimal hyperparameters.
#'    \item{final_model_cv_Px:} The canonical coordinates for X from
#'	 the cross-validated FCCA model.
#'    \item{final_model_cv_Py:} The canonical coordinates for Y from
#'	 the cross-validated FCCA model.
#'    \item{final_model_reordered:} Whether the canonical coordinates for the 
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
#'    lambda1_grid=seq(1e-3, 1e-1, length=5),
#'    lambda2_grid=seq(1e-4, 0.5, length=3),
#'    gamma1_grid=10^c(-3, -1),
#'    gamma2_grid=10^c(-3, -1), method="grid",
#'    final_model.cv=TRUE)
#'
#' # The optimal hyperparameters
#' r$opt_param
#'
#' # The cross-validated canonical correlations
#' diag(cor(r$final_model_cv_Px, r$final_model_cv_Py))
#'
#' @importFrom mlrMBO makeMBOControl
#' @importFrom data.table copy setcolorder setnames as.data.table
#' @importFrom foreach foreach `%dopar%` registerDoSEQ
#' @export
optim.cv.fcca <- function(X, Y, ndim=1, nfolds=5, folds=NULL,
   lambda1_grid=c(0, seq(1e-4, 0.002, length=4)),
   lambda2_grid=c(0, seq(1e-4, 0.002, length=4)),
   gamma1_grid=10^(-2:0), gamma2_grid=10^(-2:0),
   lambda1_bopt=c(0, 1), lambda2_bopt=c(0, 1),
   gamma1_bopt=10^c(-3, 4), gamma2_bopt=10^c(-3, 4),
   method=c("bopt", "grid"), parallel=FALSE,
   final_model=TRUE, final_model_cv=FALSE, final_model_reorder=FALSE,
   check_sign=TRUE, maxiter=1e3, svd_tol=1e-12, verbose=FALSE)
{
   method <- match.arg(method)
   if(final_model_reorder && !final_model_cv) {
      stop("final_model_cv must be TRUE if final_model_reorder is TRUE")
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
      if(is.null(lambda1_bopt)) {
         stop("lambda1_bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(is.null(lambda2_bopt)) {
         stop("lambda2_bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(is.null(gamma1_bopt)) {
         stop("gamma1_bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(is.null(gamma2_bopt)) {{
         stop("gamma2_bopt must be a vector c(min, max), or a numeric >=0")
      }

      if(all(sapply(
	 list(lambda1_bopt, lambda2_bopt, gamma1_bopt, gamma2_bopt),
	    length) == 1))
	    stop("Can't do Bayesian optimisation if lambda1_bopt, ",
	       "lambda2_bopt, gamma1_bopt, gamma2_bopt are all constants")
      }
   }

   if(is.null(lambda1_grid) || any(lambda1_grid < 0)) {
      stop("lambda1_grid must be a numeric vector of values >= 0")
   }
   if(is.null(lambda2_grid) || any(lambda2_grid < 0)) {
      stop("lambda2_grid must be a numeric vector of values >=0")
   }
   if(is.null(gamma1_grid) || any(gamma1_grid < 0)) {
      stop("gamma1_grid must be a numeric vector of values >=0")
   }
   if(is.null(gamma2_grid) || any(gamma2_grid < 0)) {
      stop("gamma2_grid must be a numeric vector of values >=0")
   }
   
   # Setup initial 'warmup' results for mlrMBO
   des.fcca.cv <- cv.fcca(X, Y, ndim=ndim,
      lambda1=lambda1_grid, lambda2=lambda2_grid,
      gamma1=gamma1_grid, gamma2=gamma2_grid,
      folds=folds, parallel=parallel, check_sign=check_sign,
      svd_tol=svd_tol, verbose=verbose)
   
   des.fcca <- copy(des.fcca.cv$result)
   setcolorder(des.fcca, c(3, 4, 1, 2))
   des.fcca[, log10_gamma1 := log10(gamma1)]
   des.fcca[, log10_gamma2 := log10(gamma2)]

   if(method == "bopt") {
      des.fcca.d <- as.data.frame(des.fcca[, 
	 .(lambda1, lambda2, gamma1, gamma2, avg.sq.cor)])
      colnames(des.fcca.d)[5] <- "y"

      # Super ugly hacksto get around the fact that mlrMBO cannot handle
      # constant hyperpameters but needs them as extra arguments
      par.set.l <- vector("list", 4)
      more.args <- list()

      if(length(lambda1_bopt) > 1) {
	 par.set.l[[1]] <- ParamHelpers::makeNumericParam("lambda1",
	    lower=lambda1_bopt[1], upper=lambda1_bopt[2])
      } else {
	 more.args <- c(more.args, c("lambda1"=lambda1_bopt))
      }

      if(length(lambda2_bopt) > 1) {
	 par.set.l[[2]] <- ParamHelpers::makeNumericParam("lambda2",
	    lower=lambda2_bopt[1], upper=lambda2_bopt[2])
      } else {
	 more.args <- c(more.args, c("lambda2"=lambda2_bopt))
      }
      
      if(length(gamma1_bopt) > 1) {
	 par.set.l[[3]] <- ParamHelpers::makeNumericParam("gamma1",
	    lower=log10(gamma1_bopt[1]), upper=log10(gamma1_bopt[2]),
	    trafo=function(x) 10^x)
      } else {
	 more.args <- c(more.args, c("gamma1"=gamma1_bopt))
      }

      if(length(gamma2_bopt) > 1) {
	 par.set.l[[4]] <- ParamHelpers::makeNumericParam("gamma2",
	    lower=log10(gamma2_bopt[1]), upper=log10(gamma2_bopt[2]),
	    trafo=function(x) 10^x)
      } else {
	 more.args <- c(more.args, c("gamma2"=gamma2_bopt))
      }

      par.set.l <- par.set.l[!sapply(par.set.l, is.null)]
      par.set <- do.call(ParamHelpers::makeParamSet, par.set.l)

      # Objective function for maximisation
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
	       folds=folds, parallel=parallel, check_sign=check_sign,
	       svd_tol=svd_tol, maxiter=maxiter, verbose=FALSE)
	    m <- r$result$avg.sq.cor
	    if(verbose) {
	       cat("x:", x, "m:", m, "\n")
	    }
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
   if(final_model) {
      mod.fcca <- fcca(X, Y, ndim=ndim,
	 lambda1=opt.param["lambda1"], lambda2=opt.param["lambda2"],
	 gamma1=opt.param["gamma1"], gamma2=opt.param["gamma2"],
	 check_sign=check_sign, verbose=verbose,
	 svd_tol=svd_tol, maxiter=maxiter)

      # It's difficult to extract the cross-validation results
      # for the final model from the optimisation process so we do it again
      # here for the optimal hyperparameters
      if(final_model_cv) {
	 mod.fcca.cv <- cv.fcca(X, Y, ndim=ndim, folds=folds,
	    lambda1=opt.param["lambda1"], lambda2=opt.param["lambda2"],
	    gamma1=opt.param["gamma1"], gamma2=opt.param["gamma2"],
	    verbose=verbose, return_models=TRUE, check_sign=check_sign,
	    svd_tol=svd_tol, maxiter=maxiter)

	 Px.cv <- Py.cv <- matrix(0, nrow(X), ndim)
	 colnames(Px.cv) <- paste0("Px", 1:ndim)
	 colnames(Py.cv) <- paste0("Py", 1:ndim)
	 for(fold in 1:nfolds) {
	    mod <- mod.fcca.cv$models[[fold]][[1]][[1]][[1]][[1]]$model
	    Px.cv[folds == fold,] <- X[folds == fold,] %*% mod$A
	    Py.cv[folds == fold,] <- Y[folds == fold,] %*% mod$B
	 }

	 # Re-order the final model's canonical vectors based on the
	 # cross-validated canonical correlations (decreasing magnitude)
	 if(final_model_reorder && ndim > 1) {
	    ord <- mod.fcca.cv$result_agg[,
	       order(r.tst.mean, decreasing=TRUE)]
	    if(any(ord != seq_along(ord))) {
	       reordered <- TRUE
	    }

	    # reorder the final model
	    mod.fcca$U <- mod.fcca$U[, ord]   
	    mod.fcca$V <- mod.fcca$V[, ord]   
	    mod.fcca$a <- mod.fcca$a[, ord]   
	    mod.fcca$b <- mod.fcca$b[, ord]   
	    mod.fcca$d <- mod.fcca$d[ord]   
	    mod.fcca$r <- mod.fcca$r[ord]   
	    mod.fcca$Px <- mod.fcca$Px[, ord]   
	    mod.fcca$Py <- mod.fcca$Py[, ord]   
	    colnames(mod.fcca$Px) <- paste0("Px", 1:ndim)
	    colnames(mod.fcca$Py) <- paste0("Py", 1:ndim)

	    # reorder the cross-validated predictions too
	    Px.cv <- Px.cv[, ord]
	    Py.cv <- Py.cv[, ord]
	    colnames(Px.cv) <- paste0("Px", 1:ndim)
	    colnames(Py.cv) <- paste0("Py", 1:ndim)
	 }
      }
   }

   if(method == "grid") {
      res <- list(
	 ndim=ndim, folds=folds, nfolds=nfolds, grid_path=des.fcca,
	 opt_param=opt.param, final_model=mod.fcca, 
	 final_model_cv=mod.fcca.cv, final_model_cv_Px=Px.cv,
	 final_model_cv_Py=Py.cv, final_model_reordered=reordered)
   } else {
      res <- list(
	 ndim=ndim, folds=folds, nfolds=nfolds,
	 grid_path=des.fcca, bopt=run.fcca,
	 bopt_path=run.fcca.path, opt_param=opt.param,
	 final_model=mod.fcca, final_model_cv=mod.fcca.cv,
	 final_model_cv_Px=Px.cv, final_model_cv_Py=Py.cv,
	 final_model_reordered=reordered)
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
      x$final_model_reordered, ")\n"))
   invisible(x)
}
