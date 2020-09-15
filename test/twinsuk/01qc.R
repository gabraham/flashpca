
#library(flashpcaR)
library(devtools)
library(data.table)
library(ggplot2)

set.seed(21982)

load_all("~/Code/flashpca/flashpcaR")

# Gene expression data
load("twinsuk_fat_expr.RData")

# Methylation data
load("twinsuk_fat_meth.RData")
load("twinsuk_fat_meth_info.RData")

d.info <- fread("d.info.txt.gz")

################################################################################
mean(is.na(expr))
mean(is.na(meth.m.imp))

################################################################################


f1 <- flashpca(expr, ndim=10, stand="sd", check_geno=FALSE)
colnames(f1$projection) <- paste0("PC", 1:ncol(f1$projection))

png("twinsuk_fat_expr_pca.png")
pairs(f1$projection[, 1:5], gap=0)
dev.off()

d.info[, ID := as.integer(gsub("^TWPID", "", indiv.id))]

d.info2 <- copy(d.info)
d.info2 <- cbind(d.info2, f1$projection)
setorder(d.info2, ID)

d.info2.m <- melt(d.info2, measure.vars=patterns("^PC"))
g1 <- ggplot(d.info2.m, aes(x=ID, y=value))
g1 <- g1 + geom_point()
g1 <- g1 + facet_wrap(variable ~ ., ncol=4, nrow=3)
g1 <- g1 + theme_bw()

ggsave(g1, file="twinsuk_fat_expr_check.png")

f2 <- flashpca(meth.m.imp, ndim=10, stand="sd", check_geno=FALSE)
colnames(f2$projection) <- paste0("PC", 1:ncol(f2$projection))

png("twinsuk_fat_meth_pca.png")
pairs(f2$projection[, 1:5], gap=0)
dev.off()

w <- sample(ncol(meth.m.imp), 1e4)
meth.m.imp.thin <- meth.m.imp[, w]
meth.chr <- d.meth.info[
   data.table(IlmnID=colnames(meth.m.imp.thin)),
   on="IlmnID", Chromosome_36]

f3 <- flashpca(meth.m.imp.thin, ndim=10, stand="sd",
   check_geno=FALSE, do_loadings=TRUE)
colnames(f3$projection) <- paste0("PC", 1:ncol(f3$projection))

png("twinsuk_fat_meth_thinned_pca.png")
pairs(f3$projection[, 1:10], gap=0)
dev.off()


d.info3 <- copy(d.info)
d.info3 <- cbind(d.info3, f3$projection)
setorder(d.info3, ID)

d.info3.m <- melt(d.info3, measure.vars=patterns("^PC"))
g2 <- ggplot(d.info3.m, aes(x=ID, y=value))
g2 <- g2 + geom_point()
g2 <- g2 + facet_wrap(variable ~ ., ncol=4, nrow=3)
g2 <- g2 + theme_bw()

ggsave(g2, file="twinsuk_fat_meth_check.png")


d.info4 <- copy(d.info)
d.info4 <- cbind(d.info4, f3$projection)
setorder(d.info4, ID)

d.info4.m <- melt(d.info4, measure.vars=patterns("^PC"))
g2 <- ggplot(d.info4.m, aes(x=ID, y=value))
g2 <- g2 + geom_point()
g2 <- g2 + facet_wrap(variable ~ ., ncol=4, nrow=3)
g2 <- g2 + theme_bw()

ggsave(g2, file="twinsuk_fat_meth_check.png")


d.info4[, ID.grp := Hmisc::cut2(ID, g=10)]
g3 <- ggplot(d.info4, aes(x=PC1, y=PC5, colour=ID.grp))
g3 <- g3 + geom_point()
g3 <- g3 + theme_bw()

ggsave(g3, file="twinsuk_fat_meth_check_PCs.png")

colnames(f3$loadings) <- paste0("V", 1:ncol(f3$loadings))
d5 <- data.table(x=1:nrow(f3$loadings), f3$loadings)
d5m <- melt(d5, id.vars="x")
g5 <- ggplot(d5m, aes(x=x, y=value))
g5 <- g5 + geom_point(alpha=0.1)
g5 <- g5 + facet_wrap(. ~ variable)
g5 <- g5 + theme_bw()

ggsave(g5, file="twinsuk_fat_meth_check_loadings.png")

png("twinsuk_fat_meth_check_loadings_pairs.png",
   width=480 * 2, height=480 * 2)
pairs(f3$loadings, gap=0, col=meth.chr)
dev.off()

