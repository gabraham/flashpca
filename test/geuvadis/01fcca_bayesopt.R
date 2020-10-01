
rm(list=ls())
graphics.off()

options(error=dump.frames)

#library(flashpcaR)
library(PMA)
library(data.table)
library(doMC)
library(mlrMBO)
library(ggplot2)
library(GGally)
library(pROC)

library(devtools)

load_all("~/Code/flashpca/flashpcaR")

registerDoMC(cores=20)

set.seed(2729)

load("exp_snps.rda")

# These SNPs have been LD pruned already, and filtered by MAF
warning("check where these SNPs came from")
warning("some of these SNPs are close to monomorphic")
vs <- apply(snps, 2, var)
snps <- snps[, vs > 0.1]

#X <- scale(snps[, sample(ncol(snps), 15000)])
Y <- scale(exp_gene)
X <- scale(snps[, sample(ncol(snps), 3500)])
#Y <- scale(exp_gene[, sample(ncol(exp_gene), 10000)])

save(X, Y, file="geuvadis_data.RData")

nfolds <- 5
folds <- sample(nfolds, nrow(X), replace=TRUE)
ndim <- 3

################################################################################
# PCCA
run.pcca <- cv.pcca(X, Y, ndim=ndim, kx=5, ky=30,
   folds=folds, final.model=TRUE)

png("geuvadis_pcca_cv_pairs_Py.png")
pairs(run.pcca$final.model.cv.Py, gap=0,
   col=factor(sample_info$pop))
dev.off()

save(run.pcca, file="geuvadis_pcca_model.RData")

#################################################################################
# FCCA
system.time({
   run.fcca <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
      verbose=FALSE, final.model.cv=TRUE, final.model.reorder=TRUE,
      parallel=TRUE)
})

png("geuvadis_fcca_cv_pairs.png")
pairs(
   cbind(run.fcca$final.model.cv.Px, run.fcca$final.model.cv.Py), gap=0,
   col=factor(sample_info$pop))
dev.off()

png("geuvadis_fcca_cv_pairs_Py.png")
pairs(run.fcca$final.model.cv.Py, gap=0,
   col=factor(sample_info$pop))
dev.off()

save(run.fcca, file="geuvadis_fcca_model.RData")

#################################################################################
## SCCA (L1 only, no L2)
run.scca <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
   gamma1.grid=0, gamma2.grid=0, 
   gamma1.bopt=0, gamma2.bopt=0,
   verbose=FALSE, final.model.cv=TRUE, final.model.reorder=TRUE,
   parallel=TRUE)

png("geuvadis_scca_cv_pairs.png")
pairs(
   cbind(run.scca$final.model.cv.Px, run.scca$final.model.cv.Py), gap=0,
   col=factor(sample_info$pop))
dev.off()

png("geuvadis_scca_cv_pairs_Py.png")
pairs(run.scca$final.model.cv.Py, gap=0,
   col=factor(sample_info$pop))
dev.off()

save(run.scca, file="geuvadis_scca_model.RData")

################################################################################
# Ridge CCA (L2 only, no L1)

run.rcca <- optim.cv.fcca(X, Y, ndim=ndim, folds=folds,
   lambda1.grid=0, lambda2.grid=0, 
   lambda1.bopt=0, lambda2.bopt=0,
   verbose=FALSE, final.model.cv=TRUE, final.model.reorder=TRUE,
   parallel=TRUE)

png("geuvadis_rcca_cv_pairs.png")
pairs(
   cbind(run.rcca$final.model.cv.Px, run.rcca$final.model.cv.Py), gap=0,
   col=factor(sample_info$pop))
dev.off()

png("geuvadis_rcca_cv_pairs_Py.png")
pairs(run.rcca$final.model.cv.Py, gap=0,
   col=factor(sample_info$pop))
dev.off()

save(run.rcca, file="geuvadis_rcca_model.RData")

#################################################################################
## PMA (sparse CCA)
#nperms <- 100
#niter <- 50
#pma.grid.len <- 12
#pen <- expand.grid(seq(1e-3, 0.999, length=pma.grid.len),
#   seq(1e-3, 0.999, length=pma.grid.len))
#p1 <- CCA.permute(X, Y, nperms=nperms, niter=niter,
#   penaltyxs=pen$Var1, penaltyzs=pen$Var2)
#p2 <- CCA(X, Y, penaltyx=p1$bestpenaltyx, penaltyz=p1$bestpenaltyz,
#   v=p1$v.init, K=ndim)
#p2.Px <- X %*% p2$u
#p2.Py <- Y %*% p2$v
#colnames(p2.Py) <- paste0("Py", 1:ncol(p2.Py))
#pen$Z <- p1$zstats
#save(p1, p2, p2.Px, p2.Py,
#   file="geuvadis_pma_results.RData")
#
#png("geuvadis_pma_pairs.png", width=900, height=900)
#pairs(cbind(p2.Px, p2.Py), gap=0, col=factor(sample_info$pop))
#dev.off()

################################################################################
# Predictions on entire data, not cross-validated
df1 <- data.table(
   pop=sample_info$pop, method="FCCA", run.fcca$final.model$Py[,1:3])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", run.pcca$final.model$Py[,1:3])
df3 <- data.table(
   pop=sample_info$pop, method="RCCA", run.rcca$final.model$Py[,1:3])
df4 <- data.table(
   pop=sample_info$pop, method="SCCA", run.scca$final.model$Py[,1:3])
df.final <- rbind(df1, df2, df3, df4)
df.final[, method_label :=
   factor(method, levels=c("PCCA", "RCCA", "SCCA", "FCCA"))]

g1 <- ggplot(df.final, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")
g1 <- g1 + guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(g1, file="geuvadis_all_models_final_by_pop.pdf", width=7, height=6)

# Cross-validated predictions
df1 <- data.table(
   pop=sample_info$pop, method="FCCA", run.fcca$final.model.cv.Py[,1:3])
df2 <- data.table(
   pop=sample_info$pop, method="PCCA", run.pcca$final.model.cv.Py[,1:3])
df3 <- data.table(
   pop=sample_info$pop, method="RCCA", run.rcca$final.model.cv.Py[,1:3])
df4 <- data.table(
   pop=sample_info$pop, method="SCCA", run.scca$final.model.cv.Py[,1:3])
df.cv <- rbind(df1, df2, df3, df4)
df.cv[, method_label :=
   factor(method, levels=c("PCCA", "RCCA", "SCCA", "FCCA"))]

g1 <- ggplot(df.cv, aes(x=Py1, y=Py2, colour=pop))
g1 <- g1 + geom_point(alpha=0.5)
g1 <- g1 + facet_wrap(. ~ method_label, scales="free")
g1 <- g1 + theme_bw()
g1 <- g1 + scale_colour_viridis_d(name="Population")
g1 <- g1 + scale_x_continuous("Canonical coordinate 1")
g1 <- g1 + scale_y_continuous("Canonical coordinate 2")
g1 <- g1 + guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(g1, file="geuvadis_all_models_cv_by_pop.pdf", width=7, height=6)

# Measure how well FIN is separated from the rest in Py2, for each method
df.final[, auc(factor(pop == "FIN") ~ Py2, quiet=TRUE), by=method]
df.final[, auc(factor(pop == "FIN") ~ Py3, quiet=TRUE), by=method]
df.cv[, auc(factor(pop == "FIN") ~ Py2, quiet=TRUE), by=method]
df.cv[, auc(factor(pop == "FIN") ~ Py3, quiet=TRUE), by=method]

