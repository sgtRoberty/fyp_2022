install.packages("picante")

library(ape)
library(dplyr)
library(tidyr)
library(tibble)
library(picante)

# setwd("")

##### Jasper Ridge plant communities #####
# Available on: https://pedrohbraga.github.io/CommunityPhylogenetics-Workshop/CommunityPhylogenetics-Workshop.html#content
# And: https://github.com/pedrohbraga/PhyloCompMethods-in-R-workshop
dir.create("../data/Jasper")
download.file("https://raw.githubusercontent.com/pedrohbraga/PhyloCompMethods-in-R-workshop/master/data/Jasper/resources/data/Jasper/jasper_data.csv", 
              "../data/Jasper/jasper_data.csv")
download.file("https://raw.githubusercontent.com/pedrohbraga/PhyloCompMethods-in-R-workshop/master/data/Jasper/resources/data/Jasper/jasper_tree.phy", 
              "../data/Jasper/jasper_tree.phy")

JasperPlants.tree <- read.tree("../data/Jasper/jasper_tree.phy")
JasperPlants.comm <- read.csv("./data/Jasper/jasper_data.csv", 
                              row.names = 1)

JasperPlants.cleanTree <- drop.tip(phy = JasperPlants.tree, 
                                   tip = setdiff(JasperPlants.tree$tip.label,
                                                 colnames(JasperPlants.comm)))

JasperPlants.picCleanTree <- match.phylo.comm(phy = JasperPlants.tree, 
                                              comm = JasperPlants.comm)$phy
JasperPlants.picCleanComm <- match.phylo.comm(phy = JasperPlants.tree, 
                                              comm = JasperPlants.comm)$comm
plot(JasperPlants.picCleanTree,
     cex = 0.4)
JasperPlants.picCleanComm[1:4, 1:4]

Jasper.PD <- pd(samp = JasperPlants.picCleanComm, 
                tree = JasperPlants.picCleanTree,
                include.root = FALSE)

head(Jasper.PD)

cor.test(Jasper.PD$PD, Jasper.PD$SR)

##### Making picante table from metadata #####

metadata <- read.csv("SITE-100_Mitogenome_Metadata_2022-05-17.csv")


picante_metadata <-  metadata %>% select(family, country) %>% count(family, country)

############ this is table of family abundance by country
picante_table <- pivot_wider(picante_metadata, names_from = family, values_from=n)


############ get rid of NA
picante_table$country <- as.character(picante_table$country)
picante_table$country[is.na(picante_table$country) == TRUE] <- "N/A"
picante_table[is.na(picante_table)] <- 0
picante_table <- column_to_rownames(picante_table, var = "country")


############ change to presence/absence

#FOR LOOP?
#for all values in table
#IF i >= 1 change value to 1






