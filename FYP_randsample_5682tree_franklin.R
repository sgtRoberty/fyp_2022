rm(list=ls())
library(ape)
library(phytools)
# Randomly sample global tree N times
ultra_allmito_5682taxa_tree <- 
  ape::read.tree("ultra_allmito_5682taxa_tree.tre")

Ntimes <- 2500

ultra_subset_1477taxa_trees <- rep(ultra_allmito_5682taxa_tree, Ntimes)
ultra_subset_1359taxa_trees <- rep(ultra_allmito_5682taxa_tree, Ntimes)

for (i in 1:Ntimes){
  ultra_subset_1477taxa_trees[[i]] <- ape::keep.tip(ultra_subset_1477taxa_trees[[i]], 
                                                    tip = sample(ultra_allmito_5682taxa_tree$tip.label,
                                                                 1477))
}

for (i in 1:Ntimes){
  ultra_subset_1359taxa_trees[[i]] <- ape::keep.tip(ultra_subset_1359taxa_trees[[i]], 
                                                    tip = sample(ultra_allmito_5682taxa_tree$tip.label,
                                                                 1359))
}

ultra_subset_1477taxa_trees_ltt95 <- ltt95(ultra_subset_1477taxa_trees, mode = "median",
                                           log = T, method = "lineages", shaded = T)
ultra_subset_1359taxa_trees_ltt95 <- ltt95(ultra_subset_1359taxa_trees, mode = "median",
                                           log = T, method = "lineages", shaded = T)
#Save the environment data
rm(ultra_subset_1477taxa_trees, ultra_subset_1359taxa_trees)
save.image(file='randsample_5682tree_ltt95.RData')
