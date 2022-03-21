install.packages(c("seqinr", "ape", "paleotree"))
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
  data$node.label <- NULL #Remove branch labels that cannot be read by MrBayes
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

#####Remove branch labels#####
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL

#####Splits function#####
listdescendants <- function(tree, n, nodes = T, tips = T, inc.n = F){
  if(nodes == F & tips == F){
    stop("Nothing to return!")
  }
  if(n < 1 | n > Ntip(tree) + tree$Nnode){
    stop("Argument to n is not a valid node or tip")
  }
  # if(inc.n == T){
  #   if(nodes == F & ! n <= Ntip(tree)) warning("inc.n = T ignored because n is an internal node and nodes = F")
  # }
  
  chs <- tree$edge[tree$edge[,1] == n, 2]
  
  out <- unlist(lapply(chs, function(ch){
    if(ch <= Ntip(tree)){
      if(tips == T){
        return(ch)
      } else if(tips == "labels"){
        return(tree$tip.label[ch])
      }
    } else {
      ot <- listdescendants(tree, ch, nodes, tips)
      if(nodes == T){
        ot <- c(ch,ot)
      }
      return(ot)
    }
  }))
  
  if(inc.n == T & (nodes == T | n <= Ntip(tree))){
    out <- c(n, out)
  }
  return(out)
}

splits <- function(phy, rettype = c("all", "split"), label.tips = F, split.detailed = F, rename = F){
  
  rettype = match.arg(rettype)
  if( rettype == "all" & split.detailed ){
    message("Warning: split.detailed = TRUE is not applicable if rettype = \"all\"")
  }
  
  # get all node numbers
  nodes <- Ntip(phy) + 1:phy$Nnode
  nodelabels <- if(rename | (! rename & !"node.label" %in% names(phy))) nodes else phy$node.label
  
  # get children below each node
  children <- c(
    # tips
    as.list(phy$tip.label),
    # internal nodes
    lapply(nodes, function(n) listdescendants(phy, n, nodes = F, 
                                              tips = if(label.tips) "labels" else T)) )
  # make sure tip names are sorted alphabetically for future comparison
  children <- lapply(children, sort)
  
  if(rettype == "all"){
    return(setNames(children[-c(1:Ntip(phy))], nodelabels))
  } else {
    # convert edge table into list of node id splits
    splits <- lapply(nodes, function(n) phy$edge[phy$edge[, 1] == n, 2])
    
    # set names of splits
    names(splits) <- nodelabels
    
    # retrieve tips below each node
    splits <- lapply(splits, function(s){
      # Get children
      ch <- children[s]
      # ensure splits are sorted alphabetically for future comparison  
      so <- order(unlist(lapply(ch, '[[', 1)))
      # output
      if(split.detailed) return(list(nodes = s[so], tips = ch[so])) else return(ch[so])
    })
    return(splits)
  }
}

backbone_nuclear_split <- splits(backbone_nuclear)

#####Create backbone constraints for MrBayes#####
library(paleotree)
?createMrBayesConstraints
createMrBayesConstraints(backbone_nuclear, partial = TRUE,
                         file = "../data/5_backbone_nuclear_constraints.nex")

