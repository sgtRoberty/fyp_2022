install.packages(c("seqinr", "ape"))
library(seqinr)
library(ape)

######Convert fasta to nexus#####
convt <- function(filename, datatype) {
  data = read.fasta(filename)
  filename = substr(filename, 1, nchar(filename)-5)
  write.nexus.data(data, file=paste(filename, "nex", sep=""),
                   format=datatype)
}

convt("../data/5_aa_supermatrix.fasta", datatype = "protein")

#####Convert newick to nexus#####
convttree <- function(filename){
  data <- ape::read.tree(filename)
  filename = substr(filename, 1, nchar(filename)-3)
  ape::write.nexus(data, file=paste(filename, "nex", sep=""))
}

convttree("../data/5_backbone_mtaa.tre")
convttree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre")

######Subset aa_supermatrix data for analysis#####
aa_supermatrix <- read.fasta("../data/5_aa_supermatrix.fasta")
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_nuclear <- ape::read.tree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre")
typeof(aa_supermatrix)
typeof(backbone_mtaa)

included <- aa_supermatrix[backbone_mtaa$tip.label]
excluded <- aa_supermatrix
excluded[backbone_mtaa$tip.label] <- NULL
random <- sample(excluded, 507)
subsetted <- c(included, random)

write.nexus.data(subsetted, file="../data/5_aa_supermatrix_1000subset.nex",
                 format="protein")
