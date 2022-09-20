library(ape)
library(phytools)

##### randsample_5682tree #####
load("../results/randsample_5682tree/randsample_5682tree_ltt95.RData")
source("../code/FYP/phytools_ltt95.R")

pdf("../results/FYP_figs/randsample_1477taxa_ltt95.pdf", 
    width = 7, height = 7.3)
plot(ultra_subset_1477taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4)
text(x=-300, y=1000, "Malaysia",font=2, cex=1.4, col="forestgreen")
ltt.lines(ultra_aa_mitogenome_Malaysia_tree, lwd = 2,
          backward=T, col="forestgreen")
dev.off()

pdf("../results/FYP_figs/randsample_1359taxa_ltt95.pdf", 
    width = 7, height = 7.3)
plot(ultra_subset_1359taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4)
text(x=-300, y=950, "Panama",font=2, cex=1.4, col="orange")
ltt.lines(ultra_aa_mitogenome_Panama_tree, lwd = 2,
          backward=T, col="orange")
dev.off()

pdf("../results/FYP_figs/randsample_both_ltt95.pdf", 
    width = 7, height = 14)
par(mfrow=c(2,1))
plot(ultra_subset_1477taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 0.85)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4, padj = 0.5)
text(x=-300, y=1000, "Malaysia",font=2, cex=1.4, col="forestgreen")
ltt.lines(ultra_aa_mitogenome_Malaysia_tree, lwd = 2,
          backward=T, col="forestgreen")
plot(ultra_subset_1359taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 0.85)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4, padj = 0.5)
text(x=-300, y=950, "Panama",font=2, cex=1.4, col="orange")
ltt.lines(ultra_aa_mitogenome_Panama_tree, lwd = 2,
          backward=T, col="orange")
dev.off()

##### stochastic_5682tree #####
load( "../results/stochastic_5682tree/treesim_phytools_ltt95_only.RData")
source("../code/FYP/phytools_ltt95.R")

pdf("../results/FYP_figs/stochastic_1477taxa_ltt95.pdf", 
    width = 7, height = 7.3)
plot(treesim_1477taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4)
text(x=-200, y=1000, "Malaysia",font=2, cex=1.4, col="forestgreen")
ltt.lines(ultra_aa_mitogenome_Malaysia_tree, lwd = 2,
          backward=T, col="forestgreen")
dev.off()

pdf("../results/FYP_figs/stochastic_1359taxa_ltt95.pdf", 
    width = 7, height = 7.3)
plot(treesim_1359taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4)
text(x=-240, y=950, "Panama",font=2, cex=1.4, col="orange")
ltt.lines(ultra_aa_mitogenome_Panama_tree, lwd = 2,
          backward=T, col="orange")
dev.off()


pdf("../results/FYP_figs/stochastic_both_ltt95.pdf", 
    width = 7, height = 14)
par(mfrow=c(2,1))
plot(treesim_1477taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 0.85)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4, padj = 0.5)
text(x=-200, y=1000, "Malaysia",font=2, cex=1.4, col="forestgreen")
ltt.lines(ultra_aa_mitogenome_Malaysia_tree, lwd = 2,
          backward=T, col="forestgreen")
plot(treesim_1359taxa_trees_ltt95, log = T, 
     xaxis="negative", shaded = T, bg = "lightgrey",
     method = "lineages", mode = "median")
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 0.85)
mtext(side=2, line=3, "ln(Number of lineages)", font=2, cex=1.4, padj = 0.5)
text(x=-240, y=950, "Panama",font=2, cex=1.4, col="orange")
ltt.lines(ultra_aa_mitogenome_Panama_tree, lwd = 2,
          backward=T, col="orange")
dev.off()
