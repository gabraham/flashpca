

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
#library(DiceKriging)

library(devtools)

load_all("~/Code/flashpca/flashpcaR")

registerDoMC(cores=10)

set.seed(2721)

load("exp_snps.rda")

# These SNPs have been LD pruned already, and filtered by MAF
maf <- colMeans(snps) / 2
summary(maf)

X <- scale(snps)
Y <- scale(exp_gene)

X <- X[, sample(ncol(X), 1000)]
Y <- Y[, sample(ncol(Y), 2000)]

################################################################################
# PCCA

n <- nrow(X)
f1 <- svd(X)
f2 <- svd(Y)
mx <- f1$d > 1e-10
my <- f2$d > 1e-10
kx <- 5 # number of SNP PCs
ky <- 30 # number of gene expression PCs

summary(colMeans(f1$u))
summary(colMeans(f2$u))
r4 <- cancor(f1$u[, 1:kx], f2$u[, 1:ky], xcenter=FALSE, ycenter=FALSE)
r4$Px <- f1$u[, 1:kx] %*% r4$xcoef
r4$Py <- f2$u[, 1:ky] %*% r4$ycoef
colnames(r4$Py) <- paste0("Py", 1:ncol(r4$Py))
colnames(r4$Px) <- paste0("Px", 1:ncol(r4$Px))

################################################################################
# FCCA, cross-validation
nfolds <- 5
folds <- sample(nfolds, nrow(X), replace=TRUE)

ndim <- 3

obj.fun <- smoof::makeSingleObjectiveFunction(
   name="scca",
   fn=function(x) {
      r <- cv.scca.ridge(X, Y, ndim=ndim, lambda1=x[1], lambda2=x[2],
	 gamma1=x[3], gamma2=x[4], folds=folds, verbose=FALSE)
      m <- r$result$avg.sq.cor
      cat("x:", x, "m:", m, "\n")
      ifelse(is.nan(m) || !is.finite(m), runif(1, 0, 1e-2), m)
   },
   par.set=makeParamSet(
      makeNumericParam("x1", lower=1e-4, upper=1e-1),
      makeNumericParam("x2", lower=1e-4, upper=1e-1),
      makeNumericParam("x3", lower=-3, upper=4, trafo=function(x) 10^x),
      makeNumericParam("x4", lower=-3, upper=4, trafo=function(x) 10^x)
   ),
   minimize=FALSE,
   noisy=TRUE
)

# Setup initial 'warmup' results for mlrMBO
run.time1 <- system.time({
   des.cv <- cv.scca.ridge(X, Y, ndim=ndim,
      lambda1=seq(1e-4, 0.002, length=3),
      lambda2=seq(1e-4, 0.002, length=3),
      gamma1=10^(-2:0), gamma2=10^(-2:0), folds=folds)
})
des <- copy(des.cv$result)
setcolorder(des, c(3, 4, 1, 2))
des[, gamma1 := log10(gamma1)]
des[, gamma2 := log10(gamma2)]
setnames(des, c("x1", "x2", "x3", "x4", "y"))
des <- as.data.frame(des)

ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters=50L)
ctrl <- setMBOControlInfill(ctrl, crit=makeMBOInfillCritEI())
surr.km <- makeMBOLearner(control=ctrl, obj.fun)

run.time2 <- system.time({
   run <- mbo(obj.fun, design=des, learner=surr.km,
      control=ctrl, show.info=TRUE)
})
run.path <- as.data.table(run$opt.path)
setnames(run.path, c("x1", "x2", "x3", "x4"),
   c("lambda1", "lambda2", "log10.gamma1", "log10.gamma2"))

save(run, file="geuvadis_mbo_run.RData")

# sanity check
s1 <- cv.scca.ridge(X, Y, ndim=ndim,
  lambda1=run$x$x1, lambda2=run$x$x2,
  gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
  folds=folds, verbose=FALSE)
c(s1$result$avg.sq.cor, run$y)

pdf("scca_ridge_bayesopt_runs.pdf")
plot(run)
dev.off()

s3 <- scca.ridge(X, Y, ndim=ndim,
  lambda1=run$x$x1, lambda2=run$x$x2,
  gamma1=10^(run$x$x3), gamma2=10^(run$x$x4))

################################################################################
# Bootstrap results

library(boot)
boot.fun <- function(original, idx)
{
   s <- scca.ridge(X[idx,], Y[idx,], ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4))

   # boot() can't handle returning of matrices, so we squeeze the matrices
   # as a single vector. Avoids having to re-run bootstrap for each of a and b and their
   # dimensions.
   # Note that as.numeric operates _by_column_
   c(as.numeric(s$a), as.numeric(s$b))
}

# If this is too big, the results will be VERY large
nreps <- 1000
res.boot <- boot(data=data.frame(idx=1:nrow(X)),
   statistic=boot.fun, R=nreps, parallel="multicore", ncpus=20)

ncolx <- ncol(X)
ncoly <- ncol(Y)
var.types <- rep(c("a", "b"), c(ncolx * ndim, ncoly * ndim))
var.names <- rep(c(colnames(X), colnames(Y)), each=ndim)

res.boot.ci <- foreach(j=seq(along=res.boot$t0), .combine="rbind") %dopar% {
   bc <- try(boot.ci(res.boot, index=j, conf=0.95,
      type=c("norm", "perc")))
   if(var.types[j] == "a") {
      dim <- floor((j - 1) / ncolx) + 1
   } else {
      dim <- floor((j - ncolx * ndim - 1) / ncoly) + 1
   }
   if(is(bc, "try-error")) {
      r <- data.table(idx=j, var.type=var.types[j], var.name=var.names[j],
	 dim=dim,
	 boot.method=NA, t0=res.boot$t0[j], lower=NA, upper=NA)
      return(r)
   }
   m <- t(sapply(bc[-(1:3)], function(x) tail(c(x), n=2)))
   colnames(m) <- c("lower", "upper")
   r <- data.table(idx=j, var.type=var.types[j], var.name=var.names[j],
      dim=dim,
      boot.method=rownames(m), t0=res.boot$t0[j], m)
   r
}

# Which SNPs have CIs that don't cross 0
res.boot.a.signif <- res.boot.ci[sign(lower) == sign(upper) 
   & boot.method == "normal" & var.type == "a",]

# Which expression probes have CIs that don't cross 0
res.boot.b.signif <- res.boot.ci[sign(lower) == sign(upper) 
   & boot.method == "normal" & var.type == "b",]

################################################################################


pdf("scca_ridge_bayesopt_final_by_pop.pdf", width=8, height=8)
pairs(cbind(s3$Px[,1:3], s3$Py[, 1:3]),
   col=factor(sample_info$pop), main="SCCA bayes optim", gap=0)
pairs(cbind(r4$Px[,1:3], r4$Py[, 1:3]),
   col=factor(sample_info$pop), main="PCCA", gap=0)
dev.off()

df1 <- data.table(
   pop=sample_info$pop, method="FCCA", s3$Py[,1:2])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", r4$Py[,1:2])
df <- rbind(df1, df2)
df[, method_label := factor(method, levels=c("PCCA", "FCCA"))]


g1 <- ggplot(df, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
cols <- c("GBR"="blue", "FIN"="orange", "YRI"="olivedrab", "TSI"="red")
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")

ggsave(g1, file="scca_ridge_bayesopt_final_by_pop_v2.pdf",
   width=7, height=3.5)

warning("test-set R^2 are not decreasing with dimension, necessarily")
warning("do we need to sort the canonical directions somehow?")



