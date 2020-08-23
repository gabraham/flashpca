

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

lambda1 <- c(0, seq(1e-4, 1e-2, length=12))
lambda2 <- c(0, seq(1e-4, 1e-2, length=12))
gamma1 <- 10^(-4:5)
gamma2 <- 10^(-4:5)

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
   des.cv <- cv.scca.ridge(X, Y, ndim=ndim,
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

run.path <- as.data.table(run$opt.path)
setnames(run.path, c("x1", "x2", "x3", "x4"),
   c("lambda1", "lambda2", "log10.gamma1", "log10.gamma2"))

s3 <- scca.ridge(X, Y, ndim=ndim,
  lambda1=run$x$x1, lambda2=run$x$x2,
  gamma1=10^(run$x$x3), gamma2=10^(run$x$x4))


boot.fun <- function(original, idx)
{
   s <- scca.ridge(X[idx,], Y[idx,], ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4))

   # boot() can't handle returning of matrices, so we squeeze the matrices
   # as a single vector. Avoids having to re-run bootstrap for each
   # of a and b and their dimensions.
   # Note that as.numeric operates _by_column_
   #c(as.numeric(s$a), as.numeric(s$b))
   as.numeric(s$a)
}

# If this is too big, the results will be VERY large
nreps <- 10000
res.boot <- boot(data=data.frame(idx=1:nrow(X)),
   statistic=boot.fun, R=nreps, parallel="multicore", ncpus=20)
save(res.boot, file="dilgom_scca_ridge_boot.RData")

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

#url <- "https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/v3/humanht-12_v3_0_r3_11283641_a_txt.zip"
#tmpdir <- tempdir()
#f <- paste0(tmpdir, basename(url))
#if(!file.exists(f)) {
#   download.file(url, f)
#}
#annot <- fread(unzip(f, "HumanHT-12_V3_0_R3_11283641_A.txt"),
#   skip=8, fill=TRUE)
#save(annot, file="annot.rda")
load("annot.rda")

setnames(res.boot.ci, "var.name", "Probe_Id")
res.boot.ci[annot, on="Probe_Id", Symbol := i.Symbol]


# Select the genes with CI not overlapping with zero
res.boot.b.signif <- res.boot.ci[sign(lower) == sign(upper)
    & boot.method == "normal" & var.type == "a",]

res.boot.b.signif[, t0_abs := abs(t0)]
setorder(res.boot.b.signif, -t0_abs)


#pdf("dilgom_scca_ridge_bayesopt_final.pdf", width=8, height=8)
#pairs(cbind(s3$Px[,1:3], s3$Py[, 1:3]), main="SCCA bayes optim", gap=0)
#dev.off()
#
#dx1 <- data.table(Probe_Id=rownames(s3$a), s3$a)
#dy1 <- data.table(Metabolite=rownames(s3$b), s3$b)
#
#
#
#dx1[annot, on="Probe_Id", Symbol := i.Symbol]

LL.mod <- data.table(
   Symbol=c("CPA3", "ENPP3", "FCER1A", "GATA2",
      "HDC", "HS.132563", "MS4A2", "MS4A3", "MS4A3", "SLC45A3", "SPRYD5"),
   Probe_Id=c("ILMN_1766551", "ILMN_1749131", "ILMN_1688423",
   "ILMN_2102670", "ILMN_1792323", "ILMN_1899034", "ILMN_1806721",
   "ILMN_1695530", "ILMN_1751625", "ILMN_1726114", "ILMN_1753648"))

LL.mod[, in.data := Probe_Id %in% colnames(rna)]

res.boot.b.signif[, in.LL.module := FALSE]
res.boot.b.signif[LL.mod, on="Probe_Id", in.LL.module := TRUE]

#setnames(dx1, paste0("V", 1:3), paste0("Component ", 1:3))
#setnames(dy1, paste0("V", 1:3), paste0("Component ", 1:3))
#
#dx1[order(abs(`Component 1`), decreasing=TRUE)[1:20], ]
#dx1[order(abs(`Component 2`), decreasing=TRUE)[1:20], ]
#dx1[order(abs(`Component 3`), decreasing=TRUE)[1:20], ]
#
#dy1[order(abs(`Component 1`), decreasing=TRUE)[1:20], ]
#dy1[order(abs(`Component 2`), decreasing=TRUE)[1:20], ]
#dy1[order(abs(`Component 3`), decreasing=TRUE)[1:20], ]

an.mod <- fread("Module_Genes_DILGOM.csv", skip=2)[, -5]
setnames(an.mod, c("Probe_Id", "Symbol", "Replicated.Module",
   "CoreModuleGenesDILGOM.YFS"))

res.boot.b.signif[an.mod, on="Probe_Id", AN.module := Replicated.Module]

#dx1d <- melt(dx1, measure.vars=grep("^Component", colnames(dx1), value=TRUE))
#dx1d[, pos := 1:.N, by=variable]
#dx1d[, type := "gene"]
#dx1d[, name := Probe_Id]
#dx1d[, col := ifelse(in.LL.module, "red", "black"), by=variable]
#
#g1 <- ggplot(dx1d, aes(sample=value))
## Need to hack the colour so as to not split the plot into group
#g1 <- g1 + stat_qq(
#   colour=dx1d$col[dx1d[, order(value), by=variable]$V1]) + stat_qq_line()
#g1 <- g1 + facet_grid(. ~ variable)
#g1 <- g1 + theme_bw()
#g1 <- g1 + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)
#ggsave(g1, file="dilgom_scca_ridge_expr.png", width=10, height=4)
#
#setnames(dy1, paste0("V", 1:3), paste0("Component ", 1:3))
#dy1d <- melt(dy1, measure.vars=grep("^Component", colnames(dy1), value=TRUE))
#g2 <- ggplot(dy1d, aes(sample=value))
## Need to hack the colour so as to not split the plot into group
#g2 <- g2 + stat_qq()
##   colour=dx1d$col[dx1d[, order(value), by=variable]$V1]) 
#g2 <- g2 + stat_qq_line()
#g2 <- g2 + facet_grid(. ~ variable)
#g2 <- g2 + theme_bw()
#g2 <- g2 + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)
#ggsave(g2, file="dilgom_scca_ridge_metab.png", width=10, height=4)
#
#stop()
#
#Pxy <- cbind(r2$Px, r2$Py)
#
#warning("these correlations are way too high, overfitting?")
#g2 <- ggpairs(data.frame(Pxy), mapping=ggplot2::aes(alpha=0.5))
#g2 <- g2 + theme_bw()
#ggsave(g2, file="dilgom_cca_Pxy.png", width=9, height=9)
#

