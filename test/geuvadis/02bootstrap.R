
rm(list=ls())
graphics.off()

options(error=dump.frames)

#library(flashpcaR)
library(mixOmics)
library(PMA)
library(data.table)
library(doMC)
library(mlrMBO)
library(ggplot2)
library(boot)
library(biomaRt)
library(GGally)

library(devtools)

load_all("~/Code/flashpca/flashpcaR")

registerDoMC(cores=4)

set.seed(2721)

load("exp_snps.rda")
load("geuvadis_data.RData")
load("geuvadis_mbo_run.RData")
load("geuvadis_pcca.RData")
load("geuvadis_fcca_final_model.RData")

ndim <- 3


################################################################################
# Bootstrap the optimal model

# This function is compatible with boot::boot in case anyone wants to use
# boot::boot 
boot.fun <- function(original, idx)
{
   s <- fcca(X[idx,], Y[idx,], ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
      V=s3$V)

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
   as.numeric(s$b)
}

# The default boot::boot function calls parallel::parLapply which seems
# slower than doMC/foreach. We construct a minimal 'boot' object which can be
# passed to boot::boot.ci
boot2 <- function(nreps)
{
   rb <- foreach(1:nreps, .combine="rbind",
	 .export=c("run", "X", "Y", "ndim")) %dopar% {
      idx <- sample(nrow(X), replace=TRUE)
      boot.fun(NULL, idx)
   }
   if(nrow(rb) < nreps) {
      warning(nreps - nrow(rb), " models did not converge")
   }
   obj <- list(t0=as.numeric(s3$b), t=rb, sim="ordinary", R=nrow(rb))
   class(obj) <- "boot"
   attr(obj, "boot_type") <- "boot"
   obj
}

# Take a bootstrap object and add more replicates
boot2.more <- function(obj, nreps)
{
   rb <- foreach(1:nreps, .combine="rbind",
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
nreps <- 100
res.boot <- boot2(nreps)
save(res.boot, file="geuvadis_fcca_boot_1.RData")

res.boot <- boot2.more(res.boot, nreps=1e3)
save(res.boot, file="geuvadis_fcca_boot.RData")

ncolx <- ncol(X)
ncoly <- ncol(Y)
#var.types <- rep(c("a", "b"), c(ncolx * ndim, ncoly * ndim))
#var.names <- rep(c(colnames(X), colnames(Y)), each=ndim)
var.types <- rep("b", ncoly * ndim)
var.names <- rep(colnames(Y), ndim)

# boot.ci can only handle one result at a time
res.boot.ci <- foreach(j=seq(along=res.boot$t0), .combine="rbind") %dopar% {
   bc <- try(boot.ci(res.boot, index=j, conf=0.95, type=c("norm", "perc")))
   #if(var.types[j] == "a") {
   #   dim <- floor((j - 1) / ncolx) + 1
   #} else {
   #   dim <- floor((j - ncolx * ndim - 1) / ncoly) + 1
   #}
   dim <- floor((j - 1) / ncoly) + 1
   if(is(bc, "try-error")) {
      r <- data.table(nreps=res.boot$R, idx=j, var.type=var.types[j],
	 var.name=var.names[j], dim=dim,
	 boot.method=NA, t0=res.boot$t0[j],
	 t.boot.mean=mean(res.boot$t[,j]), t.boot.median=mean(res.boot$t[,j]),
	 lower=NA, upper=NA)
      return(r)
   }
   m <- t(sapply(bc[-(1:3)], function(x) tail(c(x), n=2)))
   colnames(m) <- c("lower", "upper")
   r <- data.table(nreps=res.boot$R, idx=j, var.type=var.types[j],
      var.name=var.names[j], dim=dim,
      boot.method=rownames(m), t0=res.boot$t0[j], 
      t.boot.mean=mean(res.boot$t[,j]), t.boot.median=mean(res.boot$t[,j]),
      m)
   r
}

save(res.boot.ci, file="geuvadis_fcca_boot_ci.RData")

# check what the bootstrap samples look like
r <- res.boot.ci[order(abs(t.boot.median), decreasing=TRUE),
   .SD[boot.method == "normal"]][1:20, ]

png("geuvadis_expr_bootstrap_checks.png", width=900, height=900)
par(mfrow=c(4, 5))
for(j in seq(along=r$idx)) {
   x <- res.boot$t[,r$idx[j]]
   plot(density(x), main=r$var.name[j])
   abline(v=0, lty=3)
   rug(x)
}
dev.off()


