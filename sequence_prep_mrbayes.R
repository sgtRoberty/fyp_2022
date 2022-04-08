install.packages(c("seqinr", "ape", "paleotree", "BiocManager"))
BiocManager::install("ggtree")

######Convert fasta to nexus#####
library(seqinr)
library(ape)

convt <- function(filename, datatype) {
  data = read.fasta(filename)
  for (x in 1:length(data)){#Remove stop codon (*) that cannot be read by MrBayes
    for (y in 1:length(data[[x]])){
      if (data[[x]][y] == "*"){
        data[[x]][y] = "-"
      }
      else if (data[[x]][y] == "j"){
        data[[x]][y] = "(i,l)"
      }
      else if (data[[x]][y] == "b"){
        data[[x]][y] = "(d,n)"
      }
      else if (data[[x]][y] == "z"){
        data[[x]][y] = "(e,q)"
      }
    }
  } 
  filename = substr(filename, 1, nchar(filename)-5)
  write.nexus.data(data, file=paste(filename, "nex", sep=""),
                   format=datatype)
}

convt("../data/5_aa_supermatrix.fasta", datatype = "protein")

#####Convert newick to nexus#####
.getTreesFromDotdotdot <- function(...)
{
  obj <- list(...)
  if (length(obj) == 1 && !inherits(obj[[1]], "phylo")) obj <- obj[[1]]
  obj
}
write.nexus.tree <- function(..., file = "", name = "", translate = TRUE)
{
  obj <- .getTreesFromDotdotdot(...)
  ntree <- length(obj)
  cat("#NEXUS\n", file = file)
  cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),
      file = file, append = TRUE)
  
  N <- length(obj[[1]]$tip.label)
  
  cat("BEGIN TREES;\n", file = file, append = TRUE)
  if (translate) {
    cat("\tTRANSLATE\n", file = file, append = TRUE)
    obj <- .compressTipLabel(obj)
    X <- paste("\t\t", 1:N, "\t", attr(obj, "TipLabel"), ",", sep = "")
    ## We remove the last comma:
    X[length(X)] <- gsub(",", "", X[length(X)])
    cat(X, file = file, append = TRUE, sep = "\n")
    cat("\t;\n", file = file, append = TRUE)
    class(obj) <- NULL
    for (i in 1:ntree)
      obj[[i]]$tip.label <- as.character(1:N)
  } else {
    if (is.null(attr(obj, "TipLabel"))) {
      for (i in 1:ntree)
        obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
    } else {
      attr(obj, "TipLabel") <- checkLabel(attr(obj, "TipLabel"))
      obj <- .uncompressTipLabel(obj)
    }
  }
  
  if (name == ""){
    title <- names(obj)
  }
  else if (name != ""){
    title <- name
  }
  if (is.null(title))
    title <- rep("UNTITLED", ntree)
  else {
    if (any(s <- title == "")) title[s] <- "UNTITLED"
  }
  
  for (i in 1:ntree) {
    if (class(obj[[i]]) != "phylo") next
    root.tag <- if (is.rooted(obj[[i]])) "= [&R] " else "= [&U] "
    cat("\tTREE", title[i], root.tag, file = file, append = TRUE)
    cat(write.tree(obj[[i]], file = ""),
        "\n", sep = "", file = file, append = TRUE)
  }
  
  cat("END;\n", file = file, append = TRUE)
}

convttree <- function(filename, name = "", translate){
  data <- ape::read.tree(filename)
  data$node.label <- NULL #Remove branch labels that cannot be read by MrBayes
  filename = substr(filename, 1, nchar(filename)-3)
  write.nexus.tree(data, file = paste(filename, "nex", sep=""),
                   name = name, translate = translate)
}

convttree("../data/5_backbone_mtaa.tre", "mtaa", translate = T)
convttree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre",
          name = "nuclear_consensus", translate = T)

######Subset aa_supermatrix data for analysis#####
read.fasta.protein <- function(filename){ #read.fasta modified to omit stop codon (*)
  data <- read.fasta(filename)
  for (x in 1:length(data)){
    for (y in 1:length(data[[x]])){
      if (data[[x]][y] == "*"){
        data[[x]][y] = "-"
      }
      else if (data[[x]][y] == "j"){
        data[[x]][y] = "(i,l)"
      }
      else if (data[[x]][y] == "b"){
        data[[x]][y] = "(d,n)"
      }
      else if (data[[x]][y] == "z"){
        data[[x]][y] = "(e,q)"
      }
    }
  }
  return(data)
}

aa_supermatrix <- read.fasta.protein("../data/5_aa_supermatrix.fasta")
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL #Remove branch labels
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

#####Add tips at random to the tree#####
library(phytools)
mtaa_1000 <- add.random(backbone_mtaa, tips = names(random))
ggtree(mtaa_1000) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
  geom_point2(aes(subset=!isTip), color="red", size=1)

mtaa_5682 <- add.random(backbone_mtaa, tips = names(excluded))
ggtree(mtaa_5682) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
  geom_point2(aes(subset=!isTip), color="red", size=1)

write.nexus.tree(mtaa_1000, file = "../data/mtaa_1000_random.nex", 
                 name = "mtaa_1000_random")
write.nexus.tree(mtaa_5682, file = "../data/mtaa_5000.nex", 
                 name = "mtaa_5000")

#####Create backbone constraints for MrBayes#####
library(paleotree)
?createMrBayesConstraints
createMrBayesConstraints(backbone_nuclear, partial = TRUE,
                         file = "../data/5_backbone_nuclear_constraints.nex")

#####Tree manipulation#####
library(ggtree)
#Tree visualisation
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL#Remove branch labels

ggtree(backbone_mtaa) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
  geom_point2(aes(subset=!isTip), color="red", size=1)

#Tree subsetting
subset.tree <- function(filename, number, name){
  tree <- ape::read.tree(filename)
  tree$node.label <- NULL
  subset <- keep.tip(tree, sample(tree$tip.label, number))
  graph <- ggtree(subset) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
    geom_point2(aes(subset=!isTip), color="red", size=1)
  filename = substr(filename, 1, nchar(filename)-4)
  write.nexus.tree(subset, file = paste(filename, "_", number, "subset.nex", sep=""),
                   name = name)
  return(graph)
}#Subset trees in Newick format

subset.tree("../data/5_backbone_mtaa.tre", 200, name = "mtaa_200")
subset.tree("../data/5_backbone_mtaa.tre", 100, "mtaa_100")
subset.tree("../data/5_backbone_mtaa.tre", 50, "mtaa_50")

subset.tree.nex <- function(filename, number, name){
  tree <- ape::read.nexus(filename)
  tree$node.label <- NULL
  subset <- keep.tip(tree, sample(tree$tip.label, number))
  graph <- ggtree(subset) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
    geom_point2(aes(subset=!isTip), color="red", size=1)
  filename = substr(filename, 1, nchar(filename)-4)
  write.nexus.tree(subset, file = paste(filename, "_", number, "subset.nex", sep=""),
                   name = name)
  return(graph)
}#Subset trees in Nexus format

subset.tree.nex("../data/cynmix.nex.con.tre", 15, "cynmix_15subset")
