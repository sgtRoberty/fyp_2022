install.packages(c("seqinr", "ape", "paleotree", "BiocManager"))
BiocManager::install("ggtree")

###### Convert fasta to nexus #####
library(seqinr)
library(ape)

convt.fasta.nex <- function(filename, format) {
  if (format == "protein"){
    seqtype <- "AA"
  }
  else if (format == "dna"){
    seqtype <- "DNA"
  }
  data = read.fasta(filename, seqtype)
  # Remove stop codon (*) and ambiguous characters that cannot be read by MrBayes
  for (x in 1:length(data)){
    data[[x]] <- gsub("\\*", "-", 
                      gsub("J", "(I,L)", 
                           gsub("B", "(D,N)", 
                                gsub("Z", "(E,Q)", data[[x]]))))
  }
  filename = substr(filename, 1, nchar(filename)-5)
  write.nexus.data(data, file=paste(filename, "nex", sep=""),
                   format)
}

convt.fasta.nex("../data/5_aa_supermatrix.fasta", format = "protein")
convt.fasta.nex("../data/representatives.fasta", format = "protein")

##### Convert tree formats(Nexus/Newick) #####
.getTreesFromDotdotdot <- function(...)
{
  obj <- list(...)
  if (length(obj) == 1 && !inherits(obj[[1]], "phylo")) obj <- obj[[1]]
  obj
}
write.nexus.tree <- function(..., file = "", names = "", translate = TRUE)
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
  
  if (names == ""){
    title <- names(obj)
  }
  else if (names != ""){
    title <- names
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

convttree <- function(filename, name = "", 
                      translate, input = c("newick", "nexus")){
  if (input == "newick"){
    tree_obj <- ape::read.tree(filename)
    tree_obj$node.label <- NULL #Remove branch labels that cannot be read by MrBayes
    filename = substr(filename, 1, nchar(filename)-3)
    if (length(phytools::as.multiPhylo(tree_obj)) == 1){
      write.nexus.tree(tree_obj, file = paste(filename, "nex", sep=""),
                       names = name, translate)
    } else {
      write.nexus.tree(tree_obj, file = paste(filename, "nex", sep=""),
                       names = names(tree_obj), translate)
    }
  }
  else if (input == "nexus"){
    tree_obj <- ape::read.nexus(filename)
    filename = substr(filename, 1, nchar(filename)-3)
    if (length(phytools::as.multiPhylo(tree_obj)) == 1){
      ape::write.tree(tree_obj, file = paste(filename, "tre", sep=""),
                      tree.names = name)
    } else {
      ape::write.tree(tree_obj, file = paste(filename, "tre", sep=""),
                      tree.names = names(tree_obj))
    }
  }
}

convttree("../data/5_backbone_mtaa.tre", "mtaa", translate = T, input = "newick")
convttree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre",
          name = "nuclear_consensus", translate = T, input = "newick")

convttree("../data/mtaa_1000_random.nex", name = "mtaa", input = "nexus")
convttree("../data/1000subset_1000mtaa_test.run1.nex", input = "nexus")
convttree("../data/1000subset_1000mtaa_5dtest.run1.nex", input = "nexus")

###### Subset aa_supermatrix data for analysis #####
# read.fasta modified to omit stop codon (*) and ambiguous characters
read.fasta.protein <- function(filename, seqtype){ 
  data <- read.fasta(filename, seqtype)
  for (x in 1:length(data)){
    data[[x]] <- gsub("\\*", "-", 
                      gsub("J", "(I,L)", 
                           gsub("B", "(D,N)", 
                                gsub("Z", "(E,Q)", data[[x]]))))
  }
  return(data)
}

aa_supermatrix <- read.fasta.protein("../data/5_aa_supermatrix.fasta", 
                                     seqtype = "AA")
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL # Remove branch labels
backbone_nuclear <- ape::read.tree("../data/5_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID.tre")

paste(backbone_nuclear$tip.label[! backbone_nuclear$tip.label %in% AllMitoMetadata$db_id], 
      collapse = ",")

included <- aa_supermatrix[backbone_mtaa$tip.label]
excluded <- aa_supermatrix
excluded[backbone_mtaa$tip.label] <- NULL
random <- sample(excluded, 507)
subsetted <- c(included, random)
names_1000subset <- names(subsetted)

# Nexus format of supermatrix, stop codon and ambiguous characters converted
write.nexus.data(subsetted, file="../data/5_aa_supermatrix_1000subset.nex",
                 format="protein")

# Read fasta without converting stop codon and ambiguous characters
aa_supermatrix_raw <- read.fasta("../data/5_aa_supermatrix.fasta", 
                                 seqtype = "AA") # 5682 taxa, seqlen: 6372
included_raw <- aa_supermatrix_raw[backbone_mtaa$tip.label]
subsetted_raw <- aa_supermatrix_raw[names(subsetted)]

# Nexus format of 1000 subset data
write.nexus.data(subsetted_raw, file="../data/5_aa_supermatrix_1000subset_raw.nex",
                 format="protein")
write.nexus.data(included_raw, file="../data/5_aa_supermatrix_493subset_raw.nex",
                 format="protein")

# Fasta format of 1000 subset data
subsetted_raw <- read.nexus.data(file="../data/5_aa_supermatrix_1000subset_raw.nex")
subsetted_raw <- aa_supermatrix_raw[names(subsetted_raw)]
write.fasta(subsetted_raw, names = names(subsetted_raw),
          file.out = "../data/5_aa_supermatrix_1000subset_raw.fasta")

# Fasta format of 1000 subset data for RAxML
subsetted_raxml <- subsetted_raw
for (x in 1:length(subsetted_raw)){
  subsetted_raxml[[x]] <- gsub("J", "X", gsub("B", "X", gsub("Z", "X", 
                                                           subsetted_raw[[x]])))
}
write.fasta(subsetted_raxml, names = names(subsetted_raxml),
            file.out = "../data/5_aa_supermatrix_1000subset_raxml.fasta") 

# Fasta format of supermatrix for RAxML
aa_supermatrix_raxml <- aa_supermatrix_raw
for (x in 1:length(aa_supermatrix_raw)){
  aa_supermatrix_raxml[[x]] <- gsub("J", "X", gsub("B", "X", gsub("Z", "X", 
                                                                aa_supermatrix_raw[[x]])))
}
write.fasta(aa_supermatrix_raxml, names = names(aa_supermatrix_raxml),
            file.out = "../data/5_aa_supermatrix_raxml.fasta") 

##### Mitogenome metadata #####
AllMitoMetadata <- readr::read_csv("../data/SITE-100_Mitogenome_Metadata_2022-05-17.csv",
                                   col_types = cols(subgenus = "c", 
                                                    notes = "c", 
                                                    authority = "c"))
aa_supermatrix_raw <- read.fasta("../data/5_aa_supermatrix.fasta", 
                                 seqtype = "AA") # 5682 taxa, seqlen: 6372

aa_supermatrix_metadata <- filter(AllMitoMetadata, 
                                  db_id %in% names(aa_supermatrix_raw)) # 5609 taxa
length(names(aa_supermatrix_raw)[! names(aa_supermatrix_raw) %in% AllMitoMetadata$db_id])
# 73 taxa in the 5_aa_supermatrix don't have a match in AllMitoMetadata
paste(names(aa_supermatrix_raw)[! names(aa_supermatrix_raw) %in% AllMitoMetadata$db_id], 
      collapse = ',')

# Non-Coleoptera, outgroup taxa
table(aa_supermatrix_metadata$order)
non_coleoptera_metadata <- filter(aa_supermatrix_metadata, 
                                  order != "Coleoptera")
non_coleoptera_metadata <- split(non_coleoptera_metadata, 
                                 f = non_coleoptera_metadata$order)

# Taxa in four countries
table(aa_supermatrix_metadata$country)
Malaysia_metadata <- filter(aa_supermatrix_metadata, country == "Malaysia")
Panama_metadata <- filter(aa_supermatrix_metadata, country == "Panama")

table(Malaysia_metadata$family)
table(Panama_metadata$family)
table(Malaysia_metadata$superfamily)
table(Panama_metadata$superfamily)

##### Add tips at random to the tree #####
library(phytools)
library(ggtree)
quick.tree.plot <- function(treename){
  ggtree(treename) +
    geom_point2(aes(subset=!isTip), color="red", size=1)
}

mtaa_1000 <- add.random(backbone_mtaa, tips = names(random))
quick.tree.plot(mtaa_1000)

mtaa_5682 <- add.random(backbone_mtaa, tips = names(excluded))
quick.tree.plot(mtaa_5682)

write.nexus.tree(mtaa_1000, file = "../data/mtaa_1000_random.nex", 
                 name = "mtaa_1000_random")
write.nexus.tree(mtaa_5682, file = "../data/mtaa_5000.nex", 
                 name = "mtaa_5000")

##### Create backbone constraints for MrBayes #####
library(paleotree)
?createMrBayesConstraints
createMrBayesConstraints(backbone_nuclear, partial = TRUE,
                         file = "../data/5_backbone_nuclear_constraints.nex")

##### Tree manipulation #####
library(ggtree)
library(ggplot2)
## Tree visualisation
# mtaa backbone tree
backbone_mtaa <- ape::read.tree("../data/5_backbone_mtaa.tre")
backbone_mtaa$node.label <- NULL # Remove branch labels
ggtree(backbone_mtaa) +
  geom_point2(aes(subset=!isTip), color="red", size=1)

# RapidNJ tree
rapidnj_tree <- ape::read.tree("../data/rapidnj_test_5000taxa_bs.txt")
rapidnj_tree$tip.label <- gsub("['\"]", "", rapidnj_tree$tip.label)
rapidnj_tree$node.label <- NULL
rapidnj_tree_plot <- ggtree(rapidnj_tree, branch.length='none', layout="circular") %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily)) + geom_rootedge(rootedge = 30) + 
  geom_tippoint(aes(color=superfamily), alpha=0.45, size=0.4) + 
  theme(legend.position="right")
ggsave("../results/rapidnj_test_5000taxa_bs_superfamilytree.pdf", plot = rapidnj_tree_plot, 
       width = 30, height = 30)

# RapidNJ tree, no negative branch lengths
rapidnj_tree_noneg <- read.newick("../data/rapidnj_test_5000taxa_noneg.txt")
rapidnj_tree_noneg$tip.label <- gsub("['\"]", "", rapidnj_tree_noneg$tip.label)
rapidnj_tree_noneg_rooted <- root(rapidnj_tree_noneg, 
                                  outgroup = non_coleoptera_metadata$Megaloptera$db_id,
                                  resolve.root = T)
is.rooted(rapidnj_tree_noneg_rooted)

rapidnj_tree_noneg_rooted_plot <- ggtree(rapidnj_tree_noneg_rooted, size = 0.2) %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily), size = 0.2) + 
  geom_tippoint(aes(color=superfamily), alpha=0.45, size=0.2) + 
  labs(color = "Superfamily") +
  theme(legend.position="right")
ggsave("../results/rapidnj_tree_noneg_rooted_superfamily.pdf", 
       plot = rapidnj_tree_noneg_rooted_plot, 
       width = 30, height = 40)

rapidnj_tree_noneg_rooted_plot_bw <- ggtree(rapidnj_tree_noneg_rooted, 
                                            branch.length = 'none',
                                         layout="circular") %<+% 
  AllMitoMetadata + geom_rootedge(rootedge = 20) + 
  geom_tippoint(alpha=0.45, size=0.4)
ggsave("../results/rapidnj_tree_noneg_rooted_superfamily_bw.pdf", 
       plot = rapidnj_tree_noneg_rooted_plot_bw, 
       width = 20, height = 20)

# Legend of RapidNJ tree
library(grid)
library(gridExtra)
library(cowplot)
rapidnj_tree_noneg_rooted_legend <- get_legend(rapidnj_tree_noneg_rooted_plot)
grid.newpage()
grid.draw(rapidnj_tree_noneg_rooted_legend)
ggsave("../results/rapidnj_tree_noneg_rooted_superfamily_legend.pdf", 
       plot = rapidnj_tree_noneg_rooted_legend, 
       width = 3, height = 4)

# RapidNJ tree, no negative branch lengths, Malaysia and Panama local trees
rapidnj_tree_noneg_rooted_Malaysia <- keep.tip(rapidnj_tree_noneg_rooted, Malaysia_metadata$db_id)
rapidnj_tree_noneg_rooted_Panama <- keep.tip(rapidnj_tree_noneg_rooted, Panama_metadata$db_id)

rapidnj_tree_noneg_rooted_Malaysia_plot <- ggtree(rapidnj_tree_noneg_rooted_Malaysia, 
                                                  branch.length = "none",
                                                  layout = "circular") %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily)) + 
  geom_rootedge(rootedge = 20) + geom_tippoint(aes(color=superfamily), 
                                               alpha=0.45, size=0.4) + 
  theme(legend.position="none")
ggsave("../results/rapidnj_tree_noneg_rooted_Malaysia_superfamily.pdf", 
       plot = rapidnj_tree_noneg_rooted_Malaysia_plot, 
       width = 10, height = 10)
rapidnj_tree_noneg_rooted_Malaysia_plot_bw <- ggtree(rapidnj_tree_noneg_rooted_Malaysia, 
                                            branch.length = 'none',
                                            layout="circular") %<+% 
  AllMitoMetadata + geom_rootedge(rootedge = 20) + 
  geom_tippoint(alpha=0.45, size=0.4)
ggsave("../results/rapidnj_tree_noneg_rooted_Malaysia_bw.pdf", 
       plot = rapidnj_tree_noneg_rooted_Malaysia_plot_bw, 
       width = 10, height = 10)

rapidnj_tree_noneg_rooted_Panama_plot <- ggtree(rapidnj_tree_noneg_rooted_Panama, 
                                                branch.length = "none",
                                                layout = "circular") %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily)) + 
  geom_rootedge(rootedge = 20) + geom_tippoint(aes(color=superfamily), 
                                               alpha=0.45, size=0.4) + 
  theme(legend.position="none")
ggsave("../results/rapidnj_tree_noneg_rooted_Panama_superfamily.pdf", 
       plot = rapidnj_tree_noneg_rooted_Panama_plot, 
       width = 10, height = 10)
rapidnj_tree_noneg_rooted_Panama_plot_bw <- ggtree(rapidnj_tree_noneg_rooted_Panama, 
                                                     branch.length = 'none',
                                                     layout="circular") %<+% 
  AllMitoMetadata + geom_rootedge(rootedge = 20) + 
  geom_tippoint(alpha=0.45, size=0.4)
ggsave("../results/rapidnj_tree_noneg_rooted_Panama_bw.pdf", 
       plot = rapidnj_tree_noneg_rooted_Panama_plot_bw, 
       width = 10, height = 10)

# Mrbayes intermediate output tree
mb_sampletrees <- ape::read.nexus("../data/1000subset_1000mtaa_test.run1.nex")
mb_sampletrees_plot <- ggtree(mb_sampletrees[[length(mb_sampletrees)]], branch.length='none', 
                              layout="circular") %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily)) + geom_rootedge(rootedge = 5) + 
  geom_tippoint(aes(color=superfamily), alpha=0.8, size=1) + 
  theme(legend.position="right")
ggsave("../results/1000subset_1000mtaa_test.run1_superfamilytree.pdf", 
       plot = mb_sampletrees_plot, 
       width = 20, height = 20)

ape::write.tree(mb_sampletrees[[length(mb_sampletrees)]], 
                file = paste("../data/1000subset_1000mtaa_test.run1_", 
                             names(mb_sampletrees[length(mb_sampletrees)]), 
                             ".tre", sep = ""))

# Mrbayes5d intermediate output tree
mb5d_sampletrees <- ape::read.nexus("../data/1000subset_1000mtaa_5dtest.run1.nex")
mb5d_sampletrees_plot <- ggtree(mb5d_sampletrees[[length(mb5d_sampletrees)]], branch.length='none', 
                                layout="circular") %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily)) + geom_rootedge(rootedge = 5) + 
  geom_tippoint(aes(color=superfamily), alpha=0.8, size=1) + 
  theme(legend.position="right")
ggsave("../results/1000subset_1000mtaa_5dtest.run1_superfamilytree.pdf", plot = mb5d_sampletrees_plot, 
       width = 20, height = 20)

ape::write.tree(mb5d_sampletrees[[length(mb5d_sampletrees)]], 
                file = paste("../data/1000subset_1000mtaa_5dtest.run1_", 
                             names(mb5d_sampletrees[length(mb5d_sampletrees)]), 
                             ".tre", sep = ""))

## Tree subsetting
subset.tree <- function(filename, number, name){
  tree <- ape::read.tree(filename)
  tree$node.label <- NULL
  subset <- keep.tip(tree, sample(tree$tip.label, number))
  graph <- ggtree(subset) +
    geom_point2(aes(subset=!isTip), color="red", size=1)
  filename = substr(filename, 1, nchar(filename)-4)
  write.nexus.tree(subset, file = paste(filename, "_", number, "subset.nex", sep=""),
                   name = name)
  return(graph)
}# Subset trees in Newick format

subset.tree("../data/5_backbone_mtaa.tre", 200, name = "mtaa_200")
subset.tree("../data/5_backbone_mtaa.tre", 100, "mtaa_100")
subset.tree("../data/5_backbone_mtaa.tre", 50, "mtaa_50")

subset.tree.nex <- function(filename, number, name){
  tree <- ape::read.nexus(filename)
  tree$node.label <- NULL
  subset <- keep.tip(tree, sample(tree$tip.label, number))
  graph <- ggtree(subset) +
    geom_point2(aes(subset=!isTip), color="red", size=1)
  filename = substr(filename, 1, nchar(filename)-4)
  write.nexus.tree(subset, file = paste(filename, "_", number, "subset.nex", sep=""),
                   name = name)
  return(graph)
}# Subset trees in Nexus format

subset.tree.nex("../data/cynmix.nex.con.tre", 15, "cynmix_15subset")

## Rooting a tree
is.rooted(backbone_nuclear)
quick.tree.plot(backbone_nuclear)

backbone_nuclear_unrooted <- unroot(backbone_nuclear)
is.rooted(backbone_nuclear_unrooted)
quick.tree.plot(backbone_nuclear_unrooted)

backbone_nuclear_rerooted <- root(backbone_nuclear_unrooted, outgroup = "SRAA00045", 
                                  resolve.root = T)
is.rooted(backbone_nuclear_rerooted)
quick.tree.plot(backbone_nuclear_rerooted)
