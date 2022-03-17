install.packages(c("seqinr", "ape"))
library(seqinr)
library(ape)

#Convert fasta to nexus
aa_supermatrix = read.fasta("../data/5_aa_supermatrix.fasta")
write.nexus.data(aa_supermatrix, file="../data/5_aa_supermatrix.nex",
                   format="protein")

#Convert newick to nexus
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
ape::write.nexus(backbone_mtaa, file='../data/5_backbone_mtaa.nex')

backbone_nuclear <- ape::read.tree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre")
ape::write.nexus(backbone_nuclear, file='../data/5_backbone_nuclear.nex')

#Subset aa_supermatrix data for analysis
typeof(aa_supermatrix)
typeof(backbone_mtaa)
backbone_mtaa$tip.label

subsetted <- aa_supermatrix[backbone_mtaa$tip.label]
excluded <- aa_supermatrix
excluded[backbone_mtaa$tip.label] <- NULL
random <- sample(excluded, 500)
new_subset <- c(subsetted, random)

write.nexus.data(new_subset, file="../data/5_aa_supermatrix_newsubset.nex",
                 format="protein")
