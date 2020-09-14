
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

load("exp_snps.rda")

# These SNPs have been LD pruned already, and filtered by MAF
warning("check where these SNPs came from")
warning("some of these SNPs are close to monomorphic")
#maf <- colMeans(snps) / 2
#summary(maf)
vs <- apply(snps, 2, var)
snps <- snps[, vs > 0.1]

X <- scale(snps[, sample(ncol(snps), 5000)])
#Y <- scale(exp_gene[, sample(ncol(exp_gene), 7000)])
Y <- scale(exp_gene)

save(X, Y, file="geuvadis_data.RData")


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

save(r4, file="geuvadis_pcca.RData")

################################################################################
# FCCA, cross-validation
nfolds <- 5
folds <- sample(nfolds, nrow(X), replace=TRUE)

ndim <- 3

obj.fun <- smoof::makeSingleObjectiveFunction(
   name="scca",
   fn=function(x) {
      r <- cv.fcca(X, Y, ndim=ndim, lambda1=x[1], lambda2=x[2],
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
s1 <- cv.fcca(X, Y, ndim=ndim,
  lambda1=run$x$x1, lambda2=run$x$x2,
  gamma1=10^(run$x$x3), gamma2=10^(run$x$x4),
  folds=folds, verbose=FALSE)
c(s1$result$avg.sq.cor, run$y)

pdf("fcca_bayesopt_runs.pdf")
plot(run)
dev.off()


################################################################################
# Get test set predictions
folds <- sample(1:5, size=nrow(X), replace=TRUE)
Px <- Py <- matrix(0, nrow(X), ndim)
colnames(Px) <- paste0("Px", 1:ndim)
colnames(Py) <- paste0("Py", 1:ndim)
for(fold in 1:5)
{
   s3a <- fcca(X[folds != fold,], Y[folds != fold,],
      ndim=ndim, lambda1=run$x$x1, lambda2=run$x$x2,
      gamma1=10^(run$x$x3), gamma2=10^(run$x$x4))
   Px[folds == fold, ] <- X[folds == fold,] %*% s3a$a
   Py[folds == fold, ] <- Y[folds == fold,] %*% s3a$b
}

P <- cbind(Px, Py)
png("geuvadis_fcca_scatter.png", width=700, height=700)
pairs(P, gap=0, col=factor(sample_info$pop))
dev.off()

################################################################################
# Run scca over thw whole data, using the optimal penalties
system.time({
   s3 <- fcca(X, Y, ndim=ndim,
      lambda1=run$x$x1, lambda2=run$x$x2,
   gamma1=10^(run$x$x3), gamma2=10^(run$x$x4), verbose=TRUE)
})

save(s3, file="geuvadis_fcca_final_model.RData")

################################################################################
# Plot the gene expression by populations
#pdf("fcca_bayesopt_final_by_pop.pdf", width=8, height=8)
#pairs(cbind(s3$Px[,1:3], s3$Py[, 1:3]),
#   col=factor(sample_info$pop), main="SCCA bayes optim", gap=0)
#pairs(cbind(r4$Px[,1:3], r4$Py[, 1:3]),
#   col=factor(sample_info$pop), main="PCCA", gap=0)
#dev.off()

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

ggsave(g1, file="fcca_bayesopt_final_by_pop_v2.pdf",
   width=7, height=3.5)

g2 <- ggscatmat(df[method == "FCCA",], columns=3:5, color="pop")
ggsave(g2, file="fcca_bayesopt_final_by_pop_alldims_v2.pdf")


