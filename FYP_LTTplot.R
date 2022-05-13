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

##### Set-up #####
dev.off()
rm(list=ls())
setwd("")

##### Dating #####
rapidnj_tree_noneg <- read.newick("../data/rapidnj_test_5000taxa_noneg.txt")
rapidnj_tree_noneg$tip.label <- gsub("['\"]", "", rapidnj_tree_noneg$tip.label)

# Calibrating tree by assigning dates from Zhang et al to internal nodes 
calib <- rbind(cal1=c("GBDL00116","BIOD01546", "fixage", 216.1),
             cal2=c("GBDL00308","BIOD00861", "fixage", 96.5),
             cal3=c("SPSO00092","SPSO00075","fixage", 136.64),
             cal4=c("GBDL01373","BIOD02005", "fixage", 145.62),
             cal5=c("BIOD00338", "BIOD02122", "fixage", 111.91),
             cal6=c("BIOD02378","GBDL01078","fixage",80.01),
             cal7=c("BIOD00196","GBDL00422","fixage",127.55),
             cal8=c("CCCP00308","BIOD00777","fixage",123.67),
             cal9=c("GBDL01279", "GBDL01427", "fixage", 143.3),
             cal10=c("GBDL01121","GBDL01121","fixage",148.04))
colnames(calib) <- c("tax1", "tax2", "age_type", "age")

# Ultrametric tree with calibrated dates and pathD8
ultra_rapidnj_tree_noneg <- pathd8(rapidnj_tree_noneg, 
                      exec ="../../../pathd8/PATHd8.exe", 
                      calibration=calib, seql=6372)

# LTT plot
ltt.plot(ultra_rapidnj_tree_noneg$path_tree, backward = T, 
         xlab="Time before present (Ma)", 
         log="y", ylab="Number of lineages")

##### Split the tree #####
# Malaysia tree
Malaysia_tips <- Malaysia$db_id
ultra_Malaysia_tree <- keep.tip(ultra_rapidnj_tree_noneg$path_tree, Malaysia_tips)
ltt.plot(ultra_Malaysia_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Malaysia", col="forestgreen")

# Panama tree
Panama_tips <- Panama$db_id
ultra_Panama_tree <- keep.tip(ultra_rapidnj_tree_noneg$path_tree, Panama_tips)
ltt.plot(ultra_Panama_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Panama", col="red")

par(mfrow = c(1, 2))
ltt.plot(ultra_Malaysia_tree, backward=TRUE,log="y",ylab="Lineages",
                    xlab="Time before present (Millions of years)",col="forestgreen",lwd=1)
ltt.plot(ultra_Panama_tree, backward=TRUE, log="y", ylab="Lineages",
                     xlab="Time before present (Millions of years)", col="red")
par(mfrow = c(1, 1))
ltt.plot(ultra_Malaysia_tree, backward=TRUE,log="y", ylab="Lineages",
         xlab="Time before present (Millions of years)",col="forestgreen", main="Both")
ltt.lines(ultra_Panama_tree,backward=TRUE,ylab="Lineages",
          xlab="Time before present (Millions of years)",col="red")

?ltt.plot

# Generate random trees
ltt.plot(rtree_2_ultra$path_tree, backward = T,
                        ylab="Number of lineages",
                        xlab="Time before present (Millions of years)")
ltt.plot(rtree_1_ultra$path_tree, backward = T,
         ylab="Number of lineages",
         xlab="Time before present (Millions of years)")
dev.off()