
library(devtools)
library(data.table)
library(doMC)
library(ggplot2)
library(GGally)
library(boot)
library(mlrMBO)

registerDoMC(cores=10)

load_all("~/Code/flashpca/flashpcaR")

load("rna.RData")
load("metab.RData")

X <- scale(rna)
Y <- scale(metab)

set.seed(23984)

dim(X)
dim(Y)
summary(colMeans(X))
summary(colMeans(Y))

save(X, Y, file="dilgom_data.RData")

#lambda1 <- c(0, seq(1e-4, 1e-2, length=12))
#lambda2 <- c(0, seq(1e-4, 1e-2, length=12))
#gamma1 <- 10^(-4:5)
#gamma2 <- 10^(-4:5)

#r0 <- cv.scca(X, Y, lambda1=lambda1, lambda2=lambda2, nfolds=10,
#   verbose=FALSE, parallel=TRUE, standx="none", standy="none")
#r0b <- scca(X, Y, lambda1=r0$best.lambda1, lambda2=r0$best.lambda2,
#   standx="none", standy="none")
#diag(cor(r0b$Px, r0b$Py))

fout <- "dilgom_cv_scca_ridge.RData"

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
   des.cv <- cv.fcca(X, Y, ndim=ndim,
      lambda1=seq(1e-4, 0.002, length=3),
      lambda2=seq(1e-4, 0.002, length=3),
      gamma1=10^(-2:0), gamma2=10^(-2:0), folds=folds)
})
des <- copy(des.cv$result)

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

save(run, file="dilgom_mbo_run.RData")

