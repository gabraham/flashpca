
library(devtools)
library(data.table)
library(doMC)
library(ggplot2)
library(GGally)
library(boot)
library(mlrMBO)

registerDoMC(cores=5)

load_all("~/Code/flashpca/flashpcaR")

# Gene expression data
load("twinsuk_fat_expr.RData")

# Methylation data
load("twinsuk_fat_meth.RData")
load("twinsuk_fat_meth_info.RData")

load("annot.rda")

set.seed(23984)

ndim <- 1
nfolds <- 5
folds <- sample(nfolds, nrow(X), replace=TRUE)

LL.mod <- data.table(
   Symbol=c("CPA3", "ENPP3", "FCER1A", "GATA2",
      "HDC", "HS.132563", "MS4A2", "MS4A3", "MS4A3", "SLC45A3", "SPRYD5"),
   Probe_Id=c("ILMN_1766551", "ILMN_1749131", "ILMN_1688423",
      "ILMN_2102670", "ILMN_1792323", "ILMN_1899034", "ILMN_1806721",
      "ILMN_1695530", "ILMN_1751625", "ILMN_1726114", "ILMN_1753648"))

LL.mod[annot, on="Probe_Id", Chr := i.Chromosome]

# Restrict CpG to those on same chromosomes as the LL module
X.LL <- scale(expr[, LL.mod$Probe_Id])
prb <- d.meth.info[CHR %in% unique(LL.mod$Chr), IlmnID]
Y.LL <- scale(meth.m.imp[, colnames(meth.m.imp) %in% prb])

save(X.LL, Y.LL, file="twinsuk_LL_data.RData")

folds <- sample(nfolds, nrow(X.LL), replace=TRUE)

# Setup initial 'warmup' results for mlrMBO
des.LL.cv <- cv.fcca(X.LL, Y.LL, ndim=ndim,
   lambda1=seq(1e-4, 0.002, length=3),
   lambda2=seq(1e-4, 0.002, length=3),
   gamma1=10^(-2:0), gamma2=10^(-2:0),
   folds=folds, verbose=TRUE)
des.LL <- copy(des.LL.cv$result)

setcolorder(des.LL, c(3, 4, 1, 2))
des.LL[, gamma1 := log10(gamma1)]
des.LL[, gamma2 := log10(gamma2)]
setnames(des.LL, c("x1", "x2", "x3", "x4", "y"))
des.LL <- as.data.frame(des.LL)

obj.fun.LL <- smoof::makeSingleObjectiveFunction(
   name="scca",
   fn=function(x) {
      r <- cv.fcca(X.LL, Y.LL, ndim=ndim, lambda1=x[1], lambda2=x[2],
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

max.time.sec <- 120
ctrl <- makeMBOControl()
ctrl <- setMBOControlTermination(ctrl, iters=50L, time.budget=max.time.sec)
ctrl <- setMBOControlInfill(ctrl, crit=makeMBOInfillCritEI())
surr.km.LL <- makeMBOLearner(control=ctrl, obj.fun.LL)
run.LL <- mbo(obj.fun.LL, design=des.LL, learner=surr.km.LL,
   control=ctrl, show.info=TRUE)

run.LL.d <- as.data.table(run.LL$opt.path)
save(run.LL, run.LL.d, file="twinsuk_mbo_run_LL.RData")

################################################################################
# Model on all data
s3 <- fcca(X.LL, Y.LL, ndim=ndim,
   lambda1=run.LL$x$x1, lambda2=run.LL$x$x2,
   gamma1=10^(run.LL$x$x3), gamma2=10^(run.LL$x$x4),
   verbose=TRUE)
save(s3, file="twinsuk_LL_final_model.RData")

################################################################################
# Get test set predictions
folds <- sample(1:nfolds, size=nrow(X.LL), replace=TRUE)
Px <- Py <- matrix(0, nrow(X.LL), ndim)
colnames(Px) <- paste0("Px", 1:ndim)
colnames(Py) <- paste0("Py", 1:ndim)
for(fold in 1:nfolds)
{
   s3a <- fcca(X.LL[folds != fold,], Y.LL[folds != fold,],
      ndim=ndim, lambda1=run.LL$x$x1, lambda2=run.LL$x$x2,
      gamma1=10^(run.LL$x$x3), gamma2=10^(run.LL$x$x4))
   Px[folds == fold, ] <- X.LL[folds == fold,] %*% s3a$a
   Py[folds == fold, ] <- Y.LL[folds == fold,] %*% s3a$b
}

P <- cbind(Px, Py)
png("twinsuk_LL_fcca_scatter.png", width=700, height=700)
pairs(P, gap=0)
dev.off()

r.ll.sec <- cv.fcca(X.LL, Y.LL,
   lambda1=seq(1e-4, 0.002, length=10),
   lambda2=seq(1e-4, 0.002, length=10),
   gamma1=c(1e6, 10^run.LL$x$x3),
   gamma2=c(1e6, 10^run.LL$x$x4),
   nfolds=5, ndim=3)

# What's the maximum average squared cross-validated correlation
r.ll.sec$result[, max(avg.sq.cor)]

# What's the maximum average squared cross-validated correlation,
# for the high L2 penalty (approx same as assuming diagonal covariance)
r.ll.sec$result[gamma1 == 1e6 & gamma2 == 1e6, max(avg.sq.cor)]

