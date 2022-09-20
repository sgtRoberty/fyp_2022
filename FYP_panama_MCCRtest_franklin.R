rm(list=ls())
library(ape)
library(phytools)
library(laser)
Ntotal <- 10000

ultra_aa_mitogenome_Panama_tree <- 
  ape::read.tree("ultra_aa_mitogenome_Panama_tree.tre")

pdf("panama_ltt_tree_plot.pdf")
panama_ltt <- ltt(force.ultrametric(ultra_aa_mitogenome_Panama_tree,
                                      method = "extend"), plot = T,
                    show.tree = T, lwd=2, main="ultra_aa_mitogenome_Panama_tree")
dev.off()
panama_ltt

panama_mccr_1000sim <- mccr(panama_ltt,
                           rho=Ntip(ultra_aa_mitogenome_Panama_tree)/Ntotal,
                           nsim=1000)
panama_mccr_1000sim
pdf("panama_mccr_1000sim.pdf")
plot(panama_mccr_1000sim, main = expression(paste("panama_mccr_1000sim: ", 
                                                  "null distribution of ",
                                                  gamma)))
dev.off()

# LASER
panama_mccrTest_1000sim <- mccrTest(CladeSize = Ntotal,
                                      NumberMissing = Ntotal - Ntip(ultra_aa_mitogenome_Panama_tree),
                                      NumberOfReps = 1000,
                                      ObservedGamma = panama_ltt$gamma)
##### Save the environment data #####
save.image(file='panama_ltt_mccr_output.RData')

##### 5000 sims #####
panama_mccr_5000sim <- mccr(panama_ltt,
                            rho=Ntip(ultra_aa_mitogenome_Panama_tree)/Ntotal,
                            nsim=5000)
panama_mccr_5000sim
pdf("panama_mccr_5000sim.pdf")
plot(panama_mccr_5000sim, main = expression(paste("panama_mccr_5000sim: ", "null distribution of ",
                                                  gamma)))
dev.off()

panama_mccrTest_5000sim <- mccrTest(CladeSize = Ntotal,
                                      NumberMissing = Ntotal - Ntip(ultra_aa_mitogenome_Panama_tree),
                                      NumberOfReps = 5000,
                                      ObservedGamma = panama_ltt$gamma)

##### Save the environment data #####
save.image(file='panama_ltt_mccr_output.RData')
