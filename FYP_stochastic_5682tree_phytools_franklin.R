rm(list=ls())
library(ape)
library(phytools)
# Simulate pure-birth/birth-death stochastic trees
Ntimes <- 2500
speciation = 0.1355742
extinction = 0.09847907

# phytools::pbtree
allmito_5682taxa_randsim_pbtrees <- pbtree(b=speciation, d=extinction,
                                           n = 10000,
                                           nsim = Ntimes)
pb_subset_1477taxa_trees <- 
  phytools::drop.tip.multiPhylo(allmito_5682taxa_randsim_pbtrees,
                                tip = sample(allmito_5682taxa_randsim_pbtrees[[1]]$tip.label,
                                             8523)) # 10000 - 1477 = 8523
save.image(file='pbtree_ltt95.RData')
pb_subset_1477taxa_trees_ltt95 <- ltt95(pb_subset_1477taxa_trees, mode = "median",
                                        log = T, method = "lineages", shaded = T)
save.image(file='pbtree_ltt95.RData')

pb_subset_1359taxa_trees <- 
  phytools::drop.tip.multiPhylo(allmito_5682taxa_randsim_pbtrees,
                                tip = sample(allmito_5682taxa_randsim_pbtrees[[1]]$tip.label,
                                             8641)) # 10000 - 1359 = 8641
save.image(file='pbtree_ltt95.RData')
pb_subset_1359taxa_trees_ltt95 <- ltt95(pb_subset_1359taxa_trees, mode = "median",
                                        log = T, method = "lineages", shaded = T)
save.image(file='pbtree_ltt95.RData')
