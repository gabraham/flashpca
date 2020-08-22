
rm(list=ls())
graphics.off()

options(error=dump.frames)

#library(flashpcaR)
library(data.table)
library(ggplot2)
library(biomaRt)
library(GGally)

library(devtools)

load_all("~/Code/flashpca/flashpcaR")

load("geuvadis_scca_ridge_boot_ci.RData")

res.boot.b.signif <- res.boot.ci[boot.method == "percent" & sign(lower) == sign(upper),]
res.boot.b.signif[, gene.ens := gsub("\\.[0-9]*$", "", var.name)]

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

bm.res <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=res.boot.b.signif$gene.ens,
  uniqueRows=TRUE)
setDT(bm.res)
res.boot.b.signif[bm.res,
   on=c("gene.ens"="ensembl_gene_id"),
   symbol := i.external_gene_name]
res.boot.b.signif[, t0_abs := abs(t0)]
setorder(res.boot.b.signif, -t0_abs)

top.k <- 5
res.boot.b.signif.top <- res.boot.b.signif[,
   head(.SD[, .(var.name, symbol, t0)], n=top.k), by=dim]
exp_gene_top <- data.table(pop=sample_info$pop,
   exp_gene[, res.boot.b.signif.top$var.name])
setnames(exp_gene_top, c("pop", res.boot.b.signif.top$symbol))
exp_gene_top.m <- melt(exp_gene_top, id.vars="pop")

# check whether a gene appears in more than one dimension
table(colSums(res.boot.b.signif.top[, table(dim, symbol)]))

exp_gene_top.m[res.boot.b.signif.top, on=c("variable"="symbol"),
   dim := i.dim]
exp_gene_top.m[, dim_label := paste0("Dimension ", dim)]

g2 <- ggplot(exp_gene_top.m[dim %in% 1:2,], aes(x=pop, y=value, colour=pop))
g2 <- g2 + geom_violin()
g2 <- g2 + geom_point(alpha=0.1)
g2 <- g2 + facet_wrap(dim_label ~ variable, scales="free_y", ncol=top.k)
g2 <- g2 + theme_bw()
g2 <- g2 + scale_x_discrete("Population")
g2 <- g2 + scale_y_continuous("Gene expression level")
g2 <- g2 + scale_colour_viridis_d(name="Population")
g2 <- g2 + theme(legend.position="none")

ggsave(g2, file="geuvadis_scca_ridge_bayesopt_top_genes_by_pop.png", width=8)

