install.packages(c("seqinr", "ape", "paleotree"))
library(seqinr)
library(ape)

######Convert fasta to nexus#####
convt <- function(filename, datatype) {
  data = read.fasta(filename)
  for (x in 1:length(data)){#Remove stop codon (*) that cannot be read by MrBayes
    for (y in 1:length(data[[x]])){
      if (data[[x]][y] == "*"){
        data[[x]][y] = "-"
      }
    }
  } 
  filename = substr(filename, 1, nchar(filename)-5)
  write.nexus.data(data, file=paste(filename, "nex", sep=""),
                   format=datatype)
}

convt("../data/5_aa_supermatrix.fasta", datatype = "protein")

#####Remove branch labels#####
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL

#####Convert newick to nexus#####
convttree <- function(filename){
  data <- ape::read.tree(filename)
  data$node.label <- NULL #Remove branch labels that cannot be read by MrBayes
  filename = substr(filename, 1, nchar(filename)-3)
  ape::write.nexus(data, file=paste(filename, "nex", sep=""))
}

convttree("../data/5_backbone_mtaa.tre")
convttree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre")

######Subset aa_supermatrix data for analysis#####
read.fasta.protein <- function(filename){ #read.fasta modified to omit stop codon (*)
  data <- read.fasta(filename)
  for (x in 1:length(data)){
    for (y in 1:length(data[[x]])){
      if (data[[x]][y] == "*"){
        data[[x]][y] = "-"
      }
    }
  }
  return(data)
}

aa_supermatrix <- read.fasta.protein("../data/5_aa_supermatrix.fasta")
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

#####Create backbone constraints for MrBayes#####
library(paleotree)
?createMrBayesConstraints
createMrBayesConstraints(backbone_nuclear, partial = TRUE,
                         file = "../data/5_backbone_nuclear_constraints.nex")
