##### Load packages #####
library(vegan)
library(ips)
library(spdep)
library(adephylo)
library(phytools)
library(ggplot2)
library(ggtree)
library(ape)
library(castor)
library(paleotree)
library(tidyr)
library(phylo)
library(TreeTools)

?pathd8()
?ltt.plot()

##### Set-up #####
dev.off()
rm(list=ls())
#setwd("")

##### Dating #####
allmito_5682taxa_tree <- read.newick("../data/allmito_nostarttree_constr_part.raxml.lastTree.TMP")
Ntip(allmito_5682taxa_tree)

is.rooted(allmito_5682taxa_tree)
is.binary(allmito_5682taxa_tree)

# Rooting at Megaloptera + Neuroptera
# allmito_5682taxa_tree_rooted <- root(allmito_5682taxa_tree, 
#                                     outgroup = c("GBDL01738",
#                                                  "SRAA00097",
#                                                  "SRAA00098",
#                                                  "SRAA00099",
#                                                  "GBDL01734"),
#                                     resolve.root = T)

# Alternative rooting at Raphidioptera
allmito_5682taxa_tree_rooted <- root(allmito_5682taxa_tree, 
                                     outgroup = c("SRAA00101", "SRAA00102"),
                                     resolve.root = T)
is.rooted(allmito_5682taxa_tree_rooted)
is.binary(allmito_5682taxa_tree_rooted)

# Calibrating tree by assigning dates from Zhang et al 2018 to internal nodes 
calib_minmaxage <- rbind(cal1=c("GBDL01466","BGLP004", "minage", 205.6), # Haliplidae + Carabidae (fossil)
                         cal2=c("BIOD00570","GBDL00060", "minage", 155.7), # Ptiliidae + Staphylinidae (fossil)
                         cal3=c("GBDL01338","GBDL01326", "minage", 155.7), # Trogidae + Scarabaeidae (fossil)
                         cal4=c("SPSO00040","BIOD01854", "minage", 150.8), # Histeridaes + Hydrophilidae (fossil)
                         cal5=c("CERN00072","GBDL00619", "minage", 122.5), # Cerambycidae + Chrysomelidae (fossil)
                         cal6=c("GBDL00997","BIOD01856", "minage", 155.7), # Laemophloeidae + Erotylidae (fossil)
                         cal7=c("SRAA00015","BIOD00777", "minage", 155.7), # Anthribidae + Curculionidae (fossil)
                         cal8=c("GBDL01172","BIOD02122", "minage", 155.7), # Trogossitidae + Melyridae (fossil)
                         cal9=c("BIOD00705","GBDL00752", "minage", 125.5), # Ripiphoridae + Tenebrionidae (fossil)
                         cal10=c("GBDL01413","GBDL01473", "minage", 164.7), # Byrrhidae + Buprestidae (fossil)
                         cal11=c("BIOD01666","GBDL01148", "minage", 109.0), # Dryopidae + Limnichidae (fossil)
                         cal12=c("BIOD00699","GBDL01373", "minage", 155.7), # Lampyridae + Omethidae (fossil)
                         cal13=c("SRAA00099", "BGLP001", "minage", 295), # Neuropterida/Myrmeleontidae + Coleoptera/Carabidae (fossil)
                         cal14=c("SRAA00099", "BGLP001", "maxage", 323.2),
                         cal15=c("SRAA00052","SRAA00047", "minage", 155.7), # Derodontidae + Clambidae (fossil)
                         cal16=c("GBDL00308","BIOD00861", "fixage", 96.5), # Lampyridae + Rhagophthalmidae
                         cal17=c("SPSO00092","SPSO00075","fixage", 136.64), # Silvanidae + Cryptophagidae
                         cal18=c("GBDL01373","BIOD02005", "fixage", 145.62), # Omethidae + Artematopodidae
                         cal19=c("BIOD00338", "BIOD02122", "fixage", 111.91), # Prionoceridae + Melyridae
                         cal20=c("BIOD02378","GBDL01078","fixage", 80.01), # Biphyllidae + Byturidae
                         cal21=c("BIOD00196","GBDL00422","fixage", 127.55), # Nitidulidae + Kateretidae
                         cal22=c("BIOD01722","BIOD00777","fixage", 123.69), # Brentidae + Curculionidae
                         cal23=c("BIOD00193","GBDL01121","fixage", 148.04) # Ptiliidae + Hydraenidae
)
colnames(calib_minmaxage) <- c("tax1", "tax2", "age_type", "age")

ultra_allmito_5682taxa_tree <- pathd8(allmito_5682taxa_tree_rooted, 
                                      exec ="../../../pathd8/PATHd8.exe", 
                                      calibration=calib_minmaxage, seql = 6372)

# LTT plots
png("../results/5682taxa_aa_mitogenome_ltt.png", 
    width = 650, height = 650)
ltt.plot(ultra_allmito_5682taxa_tree$path_tree, 
         backward = T, log="y",
         xlab="Time before present (Ma)", 
         ylab="log(Number of lineages)",
         main=paste("5682-taxon aa mitogenome tree \n(",
                    Ntip(ultra_allmito_5682taxa_tree$path_tree),
                    " taxa)", 
                    sep = ""))
dev.off()

# Write out 5682-taxon, no zero branch length
is.ultrametric(ultra_allmito_5682taxa_tree$path_tree)
is.binary(ultra_allmito_5682taxa_tree$path_tree)
ultra_allmito_extd <- ultra_allmito_5682taxa_tree$path_tree
ultra_allmito_extd$edge.length[ultra_allmito_extd$edge.length == 0] <- 1e-8
ultra_allmito_extd <- phytools::force.ultrametric(ultra_allmito_extd,
                                                  method = "extend")
is.ultrametric(ultra_allmito_extd)

write.tree(ultra_allmito_extd,
           "../results/randsample_5682tree/ultra_allmito_5682taxa_tree.tre")

##### Split the tree #####
table(aa_supermatrix_metadata$family)
table(aa_supermatrix_metadata$superfamily)

# Malaysia tree
ultra_aa_mitogenome_Malaysia_tree <- keep.tip(ultra_allmito_5682taxa_tree$path_tree, 
                                              Malaysia_metadata$db_id)
ltt.plot(ultra_aa_mitogenome_Malaysia_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main=paste("Malaysia aa mitogenome LTT\n", 
                    Ntip(ultra_aa_mitogenome_Malaysia_tree),
                    "taxa"), 
         col="forestgreen")
dev.off()
max(branching.times(ultra_aa_mitogenome_Malaysia_tree))

# Panama tree
ultra_aa_mitogenome_Panama_tree <- keep.tip(ultra_allmito_5682taxa_tree$path_tree, 
                                            Panama_metadata$db_id)
ltt.plot(ultra_aa_mitogenome_Panama_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main=paste("Panama aa mitogenome LTT\n", 
                    Ntip(ultra_aa_mitogenome_Panama_tree),
                    "taxa"), 
         col="red")
dev.off()
max(branching.times(ultra_aa_mitogenome_Panama_tree))

# Write out time-calibrated trees
write.tree(ultra_aa_mitogenome_Malaysia_tree,
           "../results/mccr_test/ultra_aa_mitogenome_Malaysia_tree.tre")
write.tree(ultra_aa_mitogenome_Panama_tree,
           "../results/mccr_test/ultra_aa_mitogenome_Panama_tree.tre")

##### LTT plots #####
# Malaysia and Panama LTT plots
par(mfrow = c(1, 2))
ltt.plot(ultra_aa_mitogenome_Malaysia_tree, backward=TRUE, log="y", ylab="Lineages",
         xlab="Time before present (Millions of years)", main=paste("Malaysia aa mitogenome LTT \n(",
                                                                    Ntip(ultra_aa_mitogenome_Malaysia_tree), " taxa)", 
                                                                    sep = ""),
         col="forestgreen", lwd=1)
ltt.plot(ultra_aa_mitogenome_Panama_tree, backward=TRUE, log="y", ylab="Lineages",
         xlab="Time before present (Millions of years)", main=paste("Panama aa mitogenome LTT \n(",
                                                                    Ntip(ultra_aa_mitogenome_Panama_tree), " taxa)", 
                                                                    sep = ""), 
         col="red")
dev.off()

par(mfrow = c(1, 1))
ltt.plot(ultra_allmito_5682taxa_tree$path_tree, backward=TRUE, 
         log="y", ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)",
         main = "Global aa mitogenome tree \n and two local trees, Malaysia and Panama",
         col="black")
ltt.lines(ultra_aa_mitogenome_Malaysia_tree, backward=T, col="forestgreen")
ltt.lines(ultra_aa_mitogenome_Panama_tree, backward=T, col="red")
dev.off()

# Randomly sampled global tree and Malaysia tree
ltt.plot(keep.tip(ultra_allmito_5682taxa_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_allmito_5682taxa_tree$path_tree),
                               size = Ntip(ultra_aa_mitogenome_Malaysia_tree), 
                               replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global aa mitogenome tree\nand Malaysia tree (",
                                 Ntip(ultra_aa_mitogenome_Malaysia_tree)," taxa)", sep = ""))
ltt.lines(ultra_aa_mitogenome_Malaysia_tree, backward=TRUE, col="forestgreen")
dev.off()

# Randomly sampled global tree and Panama tree
ltt.plot(keep.tip(ultra_allmito_5682taxa_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_allmito_5682taxa_tree$path_tree),
                               size = Ntip(ultra_aa_mitogenome_Panama_tree), 
                               replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global aa mitogenome tree\nand Panama tree (",
                                 Ntip(ultra_aa_mitogenome_Panama_tree)," taxa)", sep = ""))
ltt.lines(ultra_aa_mitogenome_Panama_tree, backward=TRUE, col="red")
dev.off()

##### Gamma statistics #####
library(castor)
?gamma_statistic()
gamma_statistic(ultra_allmito_5682taxa_tree$path_tree)
gamma_statistic(ultra_aa_mitogenome_Malaysia_tree)
gamma_statistic(ultra_aa_mitogenome_Panama_tree)

library(ape)
?gammaStat()
gammaStat(ultra_allmito_5682taxa_tree$path_tree)
gammaStat(ultra_aa_mitogenome_Malaysia_tree)
gammaStat(ultra_aa_mitogenome_Panama_tree)

##### laser (Rabosky, 2006) superseded by more recently developed programmes #####
library(devtools)
devtools::install_github(repo="cran/laser", 
                         INSTALL_opts="--byte-compile", 
                         dependencies=FALSE)
library(laser)
gamStat(branching.times(ultra_allmito_5682taxa_tree$path_tree))
gamStat(branching.times(ultra_aa_mitogenome_Malaysia_tree))
gamStat(branching.times(ultra_aa_mitogenome_Panama_tree))

##### MCCR test#####
# FYP_malaysia_MCCRtest_franklin.R
# FYP_panama_MCCRtest_franklin.R
library(phytools)
?ltt

library(laser)
?mccrTest

##### check ultrametric trees for BAMM #####
is.ultrametric(ultra_aa_mitogenome_Malaysia_tree)
min(ultra_aa_mitogenome_Malaysia_tree$edge.length)
head(table(ultra_aa_mitogenome_Malaysia_tree$edge.length), 20)

# nnls
?phytools::force.ultrametric
tree_nnls <- phytools::force.ultrametric(ultra_aa_mitogenome_Malaysia_tree,
                                         method = "nnls")
is.binary(tree_nnls)
is.ultrametric(tree_nnls)
min(tree_nnls$edge.length)
head(table(tree_nnls$edge.length), 20)

# extend
tree_extd <- ultra_aa_mitogenome_Malaysia_tree
tree_extd$edge.length[tree_extd$edge.length == 0] <- 1e-8
head(table(tree_extd$edge.length), 20)
tree_extd <- phytools::force.ultrametric(tree_extd,
                                         method = "extend")
is.binary(tree_extd)
is.ultrametric(tree_extd)
max(tree_extd$edge.length)
min(tree_extd$edge.length)
head(table(tree_extd$edge.length), 20)
