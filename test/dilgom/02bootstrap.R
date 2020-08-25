

library(devtools)
library(data.table)
library(doMC)
library(ggplot2)
library(GGally)
library(boot)
#library(mlrMBO)

registerDoMC(cores=20)

load_all("~/Code/flashpca/flashpcaR")

#load("dilgom_data.RData")
load("dilgom_mbo_run.RData")
load("rna.RData")
load("metab.RData")

X <- scale(rna)
Y <- scale(metab)

ndim <- 3

s3 <- scca.ridge(X, Y, ndim=ndim,
  lambda1=run$x$x1, lambda2=run$x$x2,
  gamma1=10^(run$x$x3), gamma2=10^(run$x$x4))

boot.fun <- function(original, idx)
{
   s <- scca.ridge(X[idx,], Y[idx,], ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
      V=s3$V, maxiter=2000)

   # important to check convergence, otherwise model coefficients are not
   # valid
   if(!s$converged) {
      cat("model didn't converge\n")
      return(NULL)
   }

   cat(".")

   # boot() can't handle returning of matrices, so we squeeze the matrices
   # as a single vector. Avoids having to re-run bootstrap for each
   # of a and b and their dimensions.
   # Note that as.numeric operates _by_column_
   #c(as.numeric(s$a), as.numeric(s$b))
   as.numeric(s$a)
}

boot2 <- function(nreps, t0)
{
   rb <- foreach(rep=1:nreps, .combine="rbind",
         .export=c("run", "X", "Y", "ndim")) %dopar% {
      idx <- sample(nrow(X), replace=TRUE)
      boot.fun(NULL, idx)
   }
   if(nrow(rb) < nreps) {
      warning(nreps - nrow(rb), " models did not converge")
   }
   obj <- list(t0=t0, t=rb, sim="ordinary", R=nrow(rb))
   class(obj) <- "boot"
   attr(obj, "boot_type") <- "boot"
   obj
}

# Take a bootstrap object and add more replicates
boot2.more <- function(obj, nreps)
{
   rb <- foreach(rep=1:nreps, .combine="rbind",
         .export=c("run", "X", "Y", "ndim")) %dopar% {
      idx <- sample(nrow(X), replace=TRUE)
      boot.fun(NULL, idx)
   }
   if(nrow(rb) < nreps) {
      warning(nreps - nrow(rb), " models did not converge")
   }

   obj <- list(t0=obj$t0, t=rbind(obj$t, rb),
      sim="ordinary", R=obj$R + nrow(rb))
   class(obj) <- "boot"
   attr(obj, "boot_type") <- "boot"
   obj
}

# If this is too big, the results will be VERY large
#nreps <- 10000
#res.boot <- boot(data=data.frame(idx=1:nrow(X)),
#   statistic=boot.fun, R=nreps, parallel="multicore", ncpus=20)
#save(res.boot, file="dilgom_scca_ridge_boot.RData")

nreps <- 1e3
system.time({
   res.boot <- boot2(nreps, t0=s3$a)
})
save(res.boot, file="dilgom_scca_ridge_boot_1.RData")

res.boot <- boot2.more(res.boot, nreps=9000)

ncolx <- ncol(X)
ncoly <- ncol(Y)
#var.types <- rep(c("a", "b"), c(ncolx * ndim, ncoly * ndim))
#var.names <- rep(c(colnames(X), colnames(Y)), each=ndim)
var.types <- rep("a", ncolx * ndim)
var.names <- rep(colnames(X), ndim)

res.boot.ci <- foreach(j=seq(along=res.boot$t0), .combine="rbind") %dopar% {
   bc <- try(boot.ci(res.boot, index=j, conf=0.95, type=c("norm", "perc")))
   #if(var.types[j] == "a") {
   #   dim <- floor((j - 1) / ncolx) + 1
   #} else {
   #   dim <- floor((j - ncolx * ndim - 1) / ncoly) + 1
   #}
   dim <- floor((j - 1) / ncolx) + 1
   if(is(bc, "try-error")) {
      r <- data.table(nreps=nreps, idx=j, var.type=var.types[j],
         var.name=var.names[j], dim=dim,
         boot.method=NA, t0=res.boot$t0[j],
	  t.boot.mean=mean(res.boot$t[,j]),
	  t.boot.median=mean(res.boot$t[,j]),
         lower=NA, upper=NA)
      return(r)
   }
   m <- t(sapply(bc[-(1:3)], function(x) tail(c(x), n=2)))
   colnames(m) <- c("lower", "upper")
   r <- data.table(nreps=nreps, idx=j, var.type=var.types[j],
      var.name=var.names[j], dim=dim,
      boot.method=rownames(m), t0=res.boot$t0[j],
      t.boot.mean=mean(res.boot$t[,j]), t.boot.median=mean(res.boot$t[,j]), m)
   r
}

save(res.boot.ci, file="dilgom_scca_ridge_boot_ci.RData")


