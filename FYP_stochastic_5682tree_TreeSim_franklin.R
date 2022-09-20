# TreeSim
rm(list=ls())
library(ape)
library(phytools)
library(TreeSim)
Ntimes <- 2500
speciation = 0.1355742
extinction = 0.09847907

treesim_1477taxa_trees <- TreeSim::sim.bd.taxa.age(n = 1477, numbsim = Ntimes,
                                                   lambda = speciation, 
                                                   mu = extinction, 
                                                   frac = 1477/10000,
                                                   age = 268.4729)
save.image(file='treesim_ltt95.RData')

treesim_1359taxa_trees <- TreeSim::sim.bd.taxa.age(n = 1359, numbsim = Ntimes,
                                                   lambda = speciation, 
                                                   mu = extinction, 
                                                   frac = 1359/10000,
                                                   age = 321.8095)
save.image(file='treesim_ltt95.RData')

class(treesim_1477taxa_trees) <- "multiPhylo"
class(treesim_1359taxa_trees) <- "multiPhylo"
treesim_1477taxa_trees_ltt95 <- ltt95(treesim_1477taxa_trees, 
                                      log = T, xaxis="negative", 
                                      shaded = T,
                                      method = "lineages", 
                                      mode = "median")
save.image(file='treesim_ltt95.RData')
treesim_1359taxa_trees_ltt95 <- ltt95(treesim_1359taxa_trees, 
                                      log = T, xaxis="negative", 
                                      shaded = T,
                                      method = "lineages", 
                                      mode = "median")
save.image(file='treesim_ltt95.RData')
