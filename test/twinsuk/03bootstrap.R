
library(devtools)
library(data.table)
library(doMC)
library(ggplot2)
library(GGally)
library(boot)
library(gplots)
library(gridExtra)

registerDoMC(cores=20)

load_all("~/Code/flashpca/flashpcaR")

load("twinsuk_fat_meth_info.RData")
load("twinsuk_LL_data.RData")
load("twinsuk_LL_mbo_run.RData")
load("twinsuk_LL_final_model.RData")

################################################################################
# Bootstrap analysis
#
boot.fun <- function(original, idx)
{
   s <- fcca(X.LL[idx,], Y.LL[idx,], ndim=s3$ndim,
      lambda1=run.LL$x$x1, lambda2=run.LL$x$x2,
      gamma1=10^(run.LL$x$x3), gamma2=10^(run.LL$x$x4),
      V=s3$V, verbose=TRUE)

   # important to check convergence, otherwise model coefficients
   # are not valid
   if(!s$converged) {
      cat("model didn't converge\n")
      return(NULL)
   }

   cat(".")

   # boot() can't handle returning of matrices, so we squeeze the matrices
   # as a single vector. Avoids having to re-run bootstrap for each
   # of a and b and their dimensions.
   # Note that as.numeric operates _by_column_
   c(as.numeric(s$a), as.numeric(s$b))
}

boot2 <- function(nreps, t0)
{
   rb <- foreach(rep=1:nreps, .combine="rbind",
         .export=c("run.LL", "X.LL", "Y.LL", "ndim")) %dopar% {
      idx <- sample(nrow(X.LL), replace=TRUE)
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
	 .export=c("run.LL", "X.LL", "Y.LL", "ndim")) %dopar% {
      idx <- sample(nrow(X.LL), replace=TRUE)
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
nreps <- 1e2
system.time({
   res.boot <- boot2(nreps, t0=c(s3$a, s3$b))
})
save(res.boot, file="twinsuk_LL_fcca_boot_1.RData")

res.boot <- boot2.more(res.boot, nreps=1000)
save(res.boot, file="twinsuk_LL_fcca_boot.RData")

ncolx <- ncol(X.LL)
ncoly <- ncol(Y.LL)
var.types <- rep(c("a", "b"), c(ncolx * ndim, ncoly * ndim))
var.names <- rep(c(colnames(X.LL), colnames(Y.LL)), each=ndim)

################################################################################
# Bootstrap confidence intervals
#
res.boot.ci <- foreach(j=seq(along=res.boot$t0), .combine="rbind") %dopar% {
   bc <- try(boot.ci(res.boot, index=j, conf=0.95, type=c("norm", "perc")))
   if(var.types[j] == "a") {
      dim <- floor((j - 1) / ncolx) + 1
   } else {
      dim <- floor((j - ncolx * ndim - 1) / ncoly) + 1
   }
   if(is(bc, "try-error")) {
      r <- data.table(nreps=res.boot$R, idx=j, var.type=var.types[j],
         var.name=var.names[j], dim=dim,
         boot.method=NA, t0=res.boot$t0[j],
	 t.boot.mean=mean(res.boot$t[,j]),
	 t.boot.median=mean(res.boot$t[,j]),
         lower=NA, upper=NA)
      return(r)
   }
   m <- t(sapply(bc[-(1:3)], function(x) tail(c(x), n=2)))
   colnames(m) <- c("lower", "upper")
   r <- data.table(nreps=res.boot$R, idx=j,
      var.type=var.types[j],
      var.name=var.names[j], dim=dim,
      boot.method=rownames(m), t0=res.boot$t0[j],
      t.boot.mean=mean(res.boot$t[,j]),
      t.boot.median=mean(res.boot$t[,j]), m)
   r
}

save(res.boot.ci, file="twinsuk_LL_fcca_boot_ci.RData")

# Look at the bootstrap distributions for the top expression probes
r <- res.boot.ci[order(abs(t.boot.median), decreasing=TRUE),
   .SD[boot.method == "normal" & var.type == "a"]]

png("twinsuk_LL_expr_bootstrap_checks.png", width=900, height=900)
par(mfrow=c(3, 4))
for(j in seq(along=r$idx)) {
   x <- res.boot$t[,r$idx[j]]
   plot(density(x), main=r$var.name[j])
   abline(v=0, lty=3)
   rug(x)
}
dev.off()

# Look at the bootstrap distributions for the top methylation probes
r <- res.boot.ci[order(abs(t.boot.median), decreasing=TRUE),
   .SD[boot.method == "normal" & var.type == "b"]][1:20,]

png("twinsuk_LL_meth_bootstrap_checks.png", width=900, height=900)
par(mfrow=c(4, 5))
for(j in seq(along=r$idx)) {
   x <- res.boot$t[,r$idx[j]]
   plot(density(x), main=r$var.name[j])
   abline(v=0, lty=3)
   rug(x)
}
dev.off()

