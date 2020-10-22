
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

load("exp_snps.rda")
load("geuvadis_fcca_boot_ci.RData")

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

# Symbols can be duplicated or missing
res.boot.b.signif[is.na(symbol), symbol := gene.ens] 
res.boot.b.signif[, dup := .N > 1, by=symbol]
res.boot.b.signif[dup == FALSE, symbol_v2 := symbol]
res.boot.b.signif[dup == TRUE,
   symbol_v2 := paste0(symbol, " (#", seq(.N), ")"), by=symbol]

res.boot.b.signif[, t0_abs := abs(t0)]
#setorder(res.boot.b.signif, -t0_abs)
res.boot.b.signif[, t.boot.median_abs := abs(t.boot.median)]
setorder(res.boot.b.signif, -t.boot.median_abs)
res.boot.b.signif[, symbol_fac := factor(symbol_v2,
   levels=res.boot.b.signif[order(t.boot.median), symbol_v2])]

g1 <- ggplot(res.boot.b.signif[t.boot.median_abs > 0.002,],
   aes(x=symbol_fac, y=t.boot.median))
g1 <- g1 + geom_pointrange(aes(ymin=lower, ymax=upper))
g1 <- g1 + theme_bw()
g1 <- g1 + coord_flip()
g1 <- g1 + geom_hline(yintercept=0, linetype="dashed", alpha=0.5)
g1 <- g1 + scale_x_discrete("Gene")
g1 <- g1 + scale_y_continuous("Bootstrap effect size (median, 95% CIs)")
ggsave(g1, file="geuvadis_fcca_bayesopt_top_gene_bootstrap.png",
   height=15, width=5)

top.k <- 5
res.boot.b.signif.top <- res.boot.b.signif[,
   head(.SD[, .(var.name, symbol, t.boot.median_abs)], n=top.k), by=dim]
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

ggsave(g2, file="geuvadis_fcca_bayesopt_top_genes_by_pop.png", width=8)

g3 <- ggplot(exp_gene_top.m[dim %in% 1:2,], aes(x=pop, y=value, colour=pop))
g3 <- g3 + geom_violin()
g3 <- g3 + geom_point(alpha=0.1)
g3 <- g3 + facet_wrap( ~ variable, ncol=top.k, scales="free_y")
g3 <- g3 + theme_bw()
g3 <- g3 + scale_x_discrete("Population")
g3 <- g3 + scale_y_continuous("Gene expression level")
g3 <- g3 + scale_colour_viridis_d(name="Population")
g3 <- g3 + theme(legend.position="none")

ggsave(g3, file="geuvadis_fcca_bayesopt_top_genes_by_pop_dim1.png",
   width=8, height=4)

