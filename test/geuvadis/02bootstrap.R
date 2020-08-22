
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

registerDoMC(cores=20)

set.seed(2721)

load("geuvadis_data.RData")

load("geuvadis_mbo_run.RData")

ndim <- 3

################################################################################
# Run scca over thw whole data, using the optimal penalties
system.time({
   s3 <- scca.ridge(X, Y, ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
   gamma1=10^(run$x$x3), gamma2=10^(run$x$x4), verbose=TRUE)
})

################################################################################
# Plot the gene expression by populations
pdf("scca_ridge_bayesopt_final_by_pop.pdf", width=8, height=8)
pairs(cbind(s3$Px[,1:3], s3$Py[, 1:3]),
   col=factor(sample_info$pop), main="SCCA bayes optim", gap=0)
pairs(cbind(r4$Px[,1:3], r4$Py[, 1:3]),
   col=factor(sample_info$pop), main="PCCA", gap=0)
dev.off()

df1 <- data.table(
   pop=sample_info$pop, method="FCCA", s3$Py[,1:3])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", r4$Py[,1:3])
df <- rbind(df1, df2)
df[, method_label := factor(method, levels=c("PCCA", "FCCA"))]

g1 <- ggplot(df, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")
g1 <- g1 + guides(colour=guide_legend(override.aes=list(alpha=1)))

ggsave(g1, file="scca_ridge_bayesopt_final_by_pop_v2.pdf",
   width=7, height=3.5)

g2 <- ggscatmat(df[method == "FCCA",], columns=3:5, color="pop")
ggsave(g2, file="scca_ridge_bayesopt_final_by_pop_alldims_v2.pdf")

################################################################################
# Bootstrap the optimal model

# This function is compatible with boot::boot in case anyone wants to use
# boot::boot 
boot.fun <- function(original, idx)
{
   s <- scca.ridge(X[idx,], Y[idx,], ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
      V=s3$V)

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
boot2 <- function(nreps) {
   rb <- foreach(1:nreps, .combine="rbind",
	 .export=c("run", "X", "Y", "ndim")) %dopar% {
      idx <- sample(nrow(X), replace=TRUE)
      boot.fun(NULL, idx)
   }
   obj <- list(t0=as.numeric(s3$b), t=rb, sim="ordinary", R=nreps)
   class(obj) <- "boot"
   attr(obj, "boot_type") <- "boot"
   obj
}

# If this is too big, the results will be VERY large
nreps <- 5000
res.boot <- boot2(nreps)
save(res.boot, file="geuvadis_scca_ridge_boot.RData")

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
      r <- data.table(idx=j, var.type=var.types[j],
	 var.name=var.names[j], dim=dim,
	 boot.method=NA, t0=res.boot$t0[j],
	 lower=NA, upper=NA)
      return(r)
   }
   m <- t(sapply(bc[-(1:3)], function(x) tail(c(x), n=2)))
   colnames(m) <- c("lower", "upper")
   r <- data.table(idx=j, var.type=var.types[j],
      var.name=var.names[j], dim=dim,
      boot.method=rownames(m), t0=res.boot$t0[j], m)
   r
}

save(res.boot.ci, file="geuvadis_scca_ridge_boot_ci.RData")

# check whether bootstrap samples are approximately normal
png("geuvadis_expr_bootstrap_checks.png", width=900, height=900)
par(mfrow=c(5, 5))
for(j in sample(ncol(res.boot$t), 25)) {
   qqnorm(res.boot$t[,j])
   qqline(res.boot$t[,j])
}
dev.off()

