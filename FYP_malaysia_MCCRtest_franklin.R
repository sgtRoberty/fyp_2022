rm(list=ls())
library(ape)
library(phytools)
library(laser)
Ntotal <- 10000

ultra_aa_mitogenome_Malaysia_tree <- 
  ape::read.tree("ultra_aa_mitogenome_Malaysia_tree.tre")

pdf("malaysia_ltt_tree_plot.pdf")
malaysia_ltt <- ltt(force.ultrametric(ultra_aa_mitogenome_Malaysia_tree,
                                      method = "extend"), plot = T,
                    show.tree = T, lwd=2, main="ultra_aa_mitogenome_Malaysia_tree")
dev.off()
malaysia_ltt

malaysia_mccr_1000sim <- mccr(malaysia_ltt,
                              rho=Ntip(ultra_aa_mitogenome_Malaysia_tree)/Ntotal,
                              nsim=1000)
malaysia_mccr_1000sim
pdf("malaysia_mccr_1000sim.pdf")
plot(malaysia_mccr_1000sim, main = expression(paste("malaysia_mccr_1000sim: ", 
                                                    "null distribution of ",
                                                    gamma)))
dev.off()

# LASER
malaysia_mccrTest_1000sim <- mccrTest(CladeSize = Ntotal,
                                      NumberMissing = Ntotal - Ntip(ultra_aa_mitogenome_Malaysia_tree),
                                      NumberOfReps = 1000,
                                      ObservedGamma = malaysia_ltt$gamma)

##### Save the environment data #####
save.image(file='malaysia_ltt_mccr_output.RData')

##### 5000 sims #####
malaysia_mccr_5000sim <- mccr(malaysia_ltt,
                              rho=Ntip(ultra_aa_mitogenome_Malaysia_tree)/Ntotal,
                              nsim=5000)
malaysia_mccr_5000sim
pdf("malaysia_mccr_5000sim.pdf")
plot(malaysia_mccr_5000sim, main = expression(paste("malaysia_mccr_5000sim: ", 
                                                    "null distribution of ",
                                                    gamma)))
dev.off()

malaysia_mccrTest_5000sim <- mccrTest(CladeSize = Ntotal,
                                      NumberMissing = Ntotal - Ntip(ultra_aa_mitogenome_Malaysia_tree),
                                      NumberOfReps = 5000,
                                      ObservedGamma = malaysia_ltt$gamma)

##### Save the environment data #####
save.image(file='malaysia_ltt_mccr_output.RData')
