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

#####Remove branch labels#####
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL

#####Convert newick to nexus#####
.getTreesFromDotdotdot <- function(...)
{
  obj <- list(...)
  if (length(obj) == 1 && !inherits(obj[[1]], "phylo")) obj <- obj[[1]]
  obj
}

write.nexus.tree <- function(..., file = "", translate = TRUE)
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
  
  title <- names(obj)
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

convttree <- function(filename, translate){
  data <- ape::read.tree(filename)
  data$node.label <- NULL #Remove branch labels that cannot be read by MrBayes
  filename = substr(filename, 1, nchar(filename)-3)
  write.nexus.tree(data, file = paste(filename, "nex", sep=""),
                   translate = translate)
}

convttree("../data/5_backbone_mtaa.tre", translate = T)
convttree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre",
          translate = T)

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

#####Tree manipulation#####
library(ggtree)
#Tree visualisation
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL

ggtree(backbone_mtaa) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
  geom_point2(aes(subset=!isTip), color="red", size=1)

#Tree subsetting
subset.tree <- function(filename, number){
  tree <- ape::read.tree(filename)
  tree$node.label <- NULL
  subset <- keep.tip(tree, sample(tree$tip.label, number))
  graph <- ggtree(subset) + geom_text2(aes(label=label, subset=!isTip), hjust=-.2) +
    geom_point2(aes(subset=!isTip), color="red", size=1)
  filename = substr(filename, 1, nchar(filename)-4)
  write.nexus.tree(subset, file = paste(filename, "_", number, "subset.nex", sep=""))
  return(graph)
}


subset.tree("../data/5_backbone_mtaa.tre", 200)
subset.tree("../data/5_backbone_mtaa.tre", 100)
subset.tree("../data/5_backbone_mtaa.tre", 50)

cynmix <- ape::read.nexus("../data/cynmix.nex.con.tre")
cynmix <- keep.tip(cynmix, sample(cynmix$tip.label, 15))
write.nexus.tree(cynmix, file = "../data/cynmix_15subset.nex")
