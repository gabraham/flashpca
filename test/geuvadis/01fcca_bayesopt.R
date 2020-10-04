
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
X <- scale(snps[, sample(ncol(snps), 5000)])
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


warning("Strongly suspect singular vector sign flipping in CV")

pdf("fcca_check_folds.pdf")
for(fold in 1:nfolds) {
   A <- run.fcca$final.model.cv$models[[fold]][[1]][[1]][[1]][[1]]$model$a
   B <- run.fcca$final.model.cv$models[[fold]][[1]][[1]][[1]][[1]]$model$b
   Px.tmp <- X[folds == fold,] %*% A
   Py.tmp <- Y[folds == fold,] %*% B
   colnames(Px.tmp) <- paste0("Py", 1:ncol(Py.tmp))
   pairs(Py.tmp, col=factor(sample_info$pop[folds == fold]), gap=0, pch=19,
      main=paste0("Fold ", fold),
      panel=function(x, y, ...) {
	 points(x, y, ...)
	 abline(h=0, v=0, lty=3)
      })
}
dev.off()

B3 <- sapply(1:nfolds, function(fold)
   run.fcca$final.model.cv$models[[fold]][[1]][[1]][[1]][[1]]$model$b[,3])
sapply(1:nfolds, function(fold) {
    b <- run.fcca$final.model.cv$models[[fold]][[1]][[1]][[1]][[1]]$model$b[,3]
    #cy <- cor(b, t(Y[folds == fold,]))
    cy <- cor(b, t(Y[folds != fold,]))
    #cy <- cor(b, t(Y))
    rowSums(sign(cy) * cy^2)
})


Px <- Py <- matrix(0, nrow(X), ndim)
for(fold in 1:nfolds) {
   A <- run.fcca$final.model.cv$models[[fold]][[1]][[1]][[1]][[1]]$model$a
   B <- run.fcca$final.model.cv$models[[fold]][[1]][[1]][[1]][[1]]$model$b
   Px.tmp <- X[folds == fold,] %*% A
   Py.tmp <- Y[folds == fold,] %*% B

   mx <- cor(A, t(X[folds != fold,]))
   sx <- rowSums(sign(mx) * mx^2)
   my <- cor(B, t(Y[folds != fold,]))
   sy <- rowSums(sign(my) * my^2)
   
   for(j in 1:ndim) {
      cat("dim: ", j, "\n")
      cat("sx: ", sx[j], ", sy: ", sy[j], "\n")
      if(sign(sx[j]) != sign(sy[j])) {
         if(sx[j] < sy[j]) {
            cat("flipx\n")
            sx[j] <- -sx[j]
         } else {
            cat("flipy\n")
            sy[j] <- -sy[j]
         }
      }
      Px.tmp[,j] <- sign(sx[j]) * Px.tmp[,j]
      Py.tmp[,j] <- sign(sy[j]) * Py.tmp[,j]
   }

   #Px.tmp <- sweep(Px.tmp, 2, sign(sx), FUN="*")
   #Py.tmp <- sweep(Py.tmp, 2, sign(sy), FUN="*")
   Px[folds == fold, ] <- Px.tmp
   Py[folds == fold, ] <- Py.tmp
}
   
   #mx <- crossprod(Px.tmp, X[folds == fold,])
   #s.left <- rowSums(sign(mx) * mx^2)

   #my <- crossprod(Py.tmp, Y[folds == fold,])
   #s.right <- rowSums(sign(my) * my^2)

   #A.flip <- A
   #A.flip[] <- 0
   #B.flip <- B
   #B.flip[] <- 0

   #for(j in 1:ndim) {
   #   cat("dim: ", j, "\n")
   #   cat("s.left: ", s.left[j], ", s.right: ", s.right[j], "\n")
   #   if(sign(s.left[j]) != sign(s.right[j])) {
   #      if(s.left[j] < s.right[j]) {
   #         cat("flip left\n")
   #         s.left[j] <- -s.left[j]
   #      } else {
   #         cat("flip right\n")
   #         s.right[j] <- -s.right[j]
   #      }
   #   }
   #   A.flip[,j] <- sign(s.left[j]) * A[,j]
   #   B.flip[,j] <- sign(s.right[j]) * B[,j]
   #}

#   Px[folds == fold,] <- X[folds == fold, ] %*% A.flip
#   Py[folds == fold,] <- Y[folds == fold, ] %*% B.flip
#}

diag(cor(Px, run.fcca$final.model.cv.Px))
diag(cor(Py, run.fcca$final.model.cv.Py))

png("fcca_sign_flip.png", width=900, height=300)
par(mfrow=c(1, 3))
plot(run.fcca$final.model.cv.Py[,c(1,3)], col=factor(sample_info$pop),
   main="cv", pch=19)
abline(h=0, v=0, lty=2)
plot(Py[, c(1, 3)], col=factor(sample_info$pop), main="cv flipped", pch=19)
abline(h=0, v=0, lty=2)
plot(run.fcca$final.model$Py[, c(1, 2)], # been reordered
   col=factor(sample_info$pop), main="Final model", pch=19)
abline(h=0, v=0, lty=2)
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

