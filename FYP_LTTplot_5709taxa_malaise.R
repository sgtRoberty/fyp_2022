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
brlen_5709taxa_malaise_lbremoved_tree <- read.newick("../data/brlen_5709taxa_malaise_lbremoved.raxml.bestTree")
Ntip(brlen_5709taxa_malaise_lbremoved_tree) # 18622 taxa

non_coleoptera_mito_5709_metadata <- dplyr::filter(mito_5709_metadata, 
                                                   order != "Coleoptera")

is.rooted(brlen_5709taxa_malaise_lbremoved_tree)
brlen_5709taxa_malaise_lbremoved_tree_rooted <- root(brlen_5709taxa_malaise_lbremoved_tree, 
                                                     outgroup = c("GBDL01734",	
                                                                  "GBDL01737",
                                                                  "GBDL01738",
                                                                  "SRAA00097",
                                                                  "SRAA00098",
                                                                  "SRAA00101",
                                                                  "SRAA00102"), 
                                                     resolve.root = F)

# Calibrating tree by assigning dates from Zhang et al 2018 to internal nodes 
calib <- rbind(cal1=c("SRAA00052","SRAA00047", "minage", 155.7), # Derodontidae + Clambidae (fossil)
             cal2=c("GBDL00308","BIOD00861", "fixage", 96.5), # Lampyridae + Rhagophthalmidae
             cal3=c("SPSO00092","SPSO00075","fixage", 136.64), # Silvanidae + Cryptophagidae
             cal4=c("GBDL01373","BIOD02005", "fixage", 145.62), # Omethidae + Artematopodidae
             cal5=c("BIOD00338", "BIOD02122", "fixage", 111.91), # Prionoceridae + Melyridae
             cal6=c("BIOD02378","GBDL01078","fixage", 80.01), # Biphyllidae + Byturidae
             cal7=c("BIOD00196","GBDL00422","fixage", 127.55), # Nitidulidae + Kateretidae
             cal8=c("BIOD01722","BIOD00777","fixage", 123.69), # Brentidae + Curculionidae
             cal9=c("GBDL01279", "GBDL01427", "fixage", 143.43), # Geotrupidae + Passalidae
             cal10=c("BIOD00193","GBDL01121","fixage", 148.04), # Ptiliidae + Hydraenidae
             cal11=c("CERN00072","GBDL00619","minage", 122.5), # Cerambycidae + Chrysomelidae (fossil)
             cal12=c("BIOD00570", "GBDL00060", "minage", 155.7), # Ptiliidae + Staphylinidae (fossil)
             cal13=c("GBDL01338", "GBDL01326", "minage", 155.7), # Trogidae + Scarabaeidae (fossil)
             cal14=c("GBDL01172", "BIOD02122", "minage", 155.7), # Trogossitidae + Melyridae (fossil)
             cal15=c("GBDL01737", "BGLP001", "minage", 295), # Neuropterida/Myrmeleontidae + Coleoptera/Carabidae (fossil)
             cal16=c("GBDL01737", "BGLP001", "maxage", 323.2)
)

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
                         cal13=c("GBDL01737", "BGLP001", "minage", 295), # Neuropterida/Myrmeleontidae + Coleoptera/Carabidae (fossil)
                         cal14=c("GBDL01737", "BGLP001", "maxage", 323.2),
                         cal15=c("SRAA00052","SRAA00047", "minage", 155.7), # Derodontidae + Clambidae (fossil)
                         cal16=c("GBDL00308","BIOD00861", "fixage", 96.5), # Lampyridae + Rhagophthalmidae
                         cal17=c("SPSO00092","SPSO00075","fixage", 136.64), # Silvanidae + Cryptophagidae
                         cal18=c("GBDL01373","BIOD02005", "fixage", 145.62), # Omethidae + Artematopodidae
                         cal19=c("BIOD00338", "BIOD02122", "fixage", 111.91), # Prionoceridae + Melyridae
                         cal20=c("BIOD02378","GBDL01078","fixage", 80.01), # Biphyllidae + Byturidae
                         cal21=c("BIOD00196","GBDL00422","fixage", 127.55), # Nitidulidae + Kateretidae
                         cal22=c("BIOD01722","BIOD00777","fixage", 123.69), # Brentidae + Curculionidae
                         cal23=c("GBDL01279", "GBDL01427", "fixage", 143.43), # Geotrupidae + Passalidae
                         cal24=c("BIOD00193","GBDL01121","fixage", 148.04) # Ptiliidae + Hydraenidae
)

colnames(calib) <- c("tax1", "tax2", "age_type", "age")
colnames(calib_minmaxage) <- c("tax1", "tax2", "age_type", "age")

# Ultrametric tree with calibrated dates and pathD8
ultra_5709taxa_malaise_tree <- pathd8(brlen_5709taxa_malaise_lbremoved_tree, 
                                      exec ="../../../pathd8/PATHd8.exe", 
                                      calibration=calib_minmaxage, seql = 18603)

ultra_5709taxa_malaise_tree_rooted <- pathd8(brlen_5709taxa_malaise_lbremoved_tree_rooted, 
                                             exec ="../../../pathd8/PATHd8.exe", 
                                             calibration=calib_minmaxage, seql = 18603)

# LTT plot
png("../results/global_malaise_trap_ltt.png", 
    width = 650, height = 650)
ltt.plot(ultra_5709taxa_malaise_tree$path_tree, 
         backward = T, log="y",
         xlab="Time before present (Ma)", 
         ylab="log(Number of lineages)",
         main=paste("Global malaise trap tree \n(",
                    Ntip(ultra_5709taxa_malaise_tree$path_tree),
                    " taxa)", 
                    sep = ""))
dev.off()

png("../results/global_malaise_trap_rooted_ltt.png", 
    width = 650, height = 650)
ltt.plot(ultra_5709taxa_malaise_tree_rooted$path_tree, 
         backward = T, log="y",
         xlab="Time before present (Ma)", 
         ylab="log(Number of lineages)",
         main=paste("Global malaise trap tree rooted\n(",
                    Ntip(ultra_5709taxa_malaise_tree$path_tree),
                    " taxa)", 
                    sep = ""))
dev.off()

##### Split the tree #####
Malaysia_tips <- c(filter(mito_5709_metadata, country == "Malaysia")$db_id,
                   filter(MalaiseTrapMetadata, country == "Malaysia")$db_id)

Panama_tips <- c(filter(mito_5709_metadata, country == "Panama")$db_id,
                 filter(MalaiseTrapMetadata, country == "Panama")$db_id)

China_tips <- c(filter(mito_5709_metadata, country == "China")$db_id,
                filter(MalaiseTrapMetadata, country == "China")$db_id)

Canada_tips <- c(filter(mito_5709_metadata, country == "Canada")$db_id,
                 filter(MalaiseTrapMetadata, country == "Canada")$db_id)

CostaRica_tips <- c(filter(mito_5709_metadata, country == "Costa Rica")$db_id,
                 filter(MalaiseTrapMetadata, country == "Costa Rica")$db_id)

Honduras_tips <- c(filter(mito_5709_metadata, country == "Honduras")$db_id,
                    filter(MalaiseTrapMetadata, country == "Honduras")$db_id)

Australia_tips <- c(filter(mito_5709_metadata, country == "Australia")$db_id,
                    filter(MalaiseTrapMetadata, country == "Australia")$db_id)

USA_tips <- c(filter(mito_5709_metadata, country == "USA")$db_id,
                    filter(MalaiseTrapMetadata, country == "United States")$db_id)

Madagascar_tips <- c(filter(mito_5709_metadata, country == "Madagascar")$db_id,
              filter(MalaiseTrapMetadata, country == "Madagascar")$db_id)

Malaysia_tips <- Malaysia_tips[! Malaysia_tips %in% removed_longbr_names]
CostaRica_tips <- CostaRica_tips[! CostaRica_tips %in% removed_longbr_names]
Honduras_tips <- Honduras_tips[! Honduras_tips %in% removed_longbr_names]

# Malaysia tree
ultra_Malaysia_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, Malaysia_tips)
ltt.plot(ultra_Malaysia_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Malaysia", col="forestgreen")
dev.off()

# Panama tree
ultra_Panama_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, Panama_tips)
ltt.plot(ultra_Panama_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Panama", col="red")
dev.off()

# China tree
ultra_China_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, China_tips)
ltt.plot(ultra_China_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="China", col="blue")
dev.off()

# Canada tree
ultra_Canada_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, Canada_tips)
ltt.plot(ultra_Canada_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Canada", col="purple")
dev.off()

# Costa Rica tree
ultra_CostaRica_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, CostaRica_tips)
ltt.plot(ultra_CostaRica_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Costa Rica", col="orange")
dev.off()

# Honduras tree
ultra_Honduras_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, Honduras_tips)
ltt.plot(ultra_Honduras_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Honduras", col="cyan")
dev.off()

# Australia tree
ultra_Australia_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, Australia_tips)
ltt.plot(ultra_Australia_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Australia", col="brown")
dev.off()

# USA tree
ultra_USA_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, USA_tips)
ltt.plot(ultra_USA_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="USA", col="darkorchid")
dev.off()

# Madagascar tree
ultra_Madagascar_tree <- keep.tip(ultra_5709taxa_malaise_tree$path_tree, Madagascar_tips)
ltt.plot(ultra_Madagascar_tree, backward = T, log="y", 
         ylab="Number of lineages", 
         xlab="Time before present (Millions of years)",
         main="Madagascar", col="darkgreen")
dev.off()

###### LTT plots #####
# Multiple LTT plots - Malaysia and Panama
par(mfrow = c(1, 2))
ltt.plot(ultra_Malaysia_tree, backward=TRUE, log="y", ylab="Lineages",
         xlab="Time before present (Millions of years)", main="Malaysia",
         col="forestgreen", lwd=1)
ltt.plot(ultra_Panama_tree, backward=TRUE, log="y", ylab="Lineages",
         xlab="Time before present (Millions of years)", main="Panama",
         col="red")
dev.off()

# Multiple LTT plots - Four localities
par(mfrow = c(1, 1))
ltt.plot(ultra_5709taxa_malaise_tree$path_tree, backward=TRUE, 
         log="y", ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", col="black")
ltt.lines(ultra_Malaysia_tree, backward=TRUE, col="forestgreen")
ltt.lines(ultra_Panama_tree, backward=TRUE, col="red")
ltt.lines(ultra_China_tree, backward=TRUE, col="blue")
ltt.lines(ultra_Canada_tree, backward=TRUE, col="purple")
dev.off()

# Randomly sampled global tree and Malaysia tree
par(mfrow = c(1, 1))
png("../results/randsampled_global_malaysia_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_Malaysia_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Malaysia tree \n(",
                                 Ntip(ultra_Malaysia_tree)," taxa)", sep = ""))
ltt.lines(ultra_Malaysia_tree, backward=TRUE, col="forestgreen")
dev.off()

# Randomly sampled global tree and Panama tree
png("../results/randsampled_global_panama_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_Panama_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Panama tree \n(",
                                 Ntip(ultra_Panama_tree)," taxa)", sep = ""))
ltt.lines(ultra_Panama_tree, backward=TRUE, col="red")
dev.off()

# Randomly sampled global tree and Canada tree
png("../results/randsampled_global_canada_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_Canada_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Canada tree \n(",
                                 Ntip(ultra_Canada_tree)," taxa)", sep = ""))
ltt.lines(ultra_Canada_tree, backward=TRUE, col="purple")
dev.off()

# Randomly sampled global tree and Costa Rica tree
png("../results/randsampled_global_costarica_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_CostaRica_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Costa Rica tree \n(",
                                 Ntip(ultra_CostaRica_tree)," taxa)", sep = ""))
ltt.lines(ultra_CostaRica_tree, backward=TRUE, col="orange")
dev.off()

# Randomly sampled global tree and China tree
png("../results/randsampled_global_china_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_China_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and China tree \n(",
                                 Ntip(ultra_China_tree)," taxa)", sep = ""))
ltt.lines(ultra_China_tree, backward=TRUE, col="blue")
dev.off()

# Randomly sampled global tree and Honduras tree
png("../results/randsampled_global_honduras_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_Honduras_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Honduras tree \n(",
                                 Ntip(ultra_Honduras_tree)," taxa)", sep = ""))
ltt.lines(ultra_Honduras_tree, backward=TRUE, col="cyan")
dev.off()

# Randomly sampled global tree and Australia tree
png("../results/randsampled_global_australia_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_Australia_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Australia tree \n(",
                                 Ntip(ultra_Australia_tree)," taxa)", sep = ""))
ltt.lines(ultra_Australia_tree, backward=TRUE, col="brown")
dev.off()

# Randomly sampled global tree and USA tree
png("../results/randsampled_global_usa_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_USA_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and USA tree \n(",
                                 Ntip(ultra_USA_tree)," taxa)", sep = ""))
ltt.lines(ultra_USA_tree, backward=TRUE, col="darkorchid")
dev.off()

# Randomly sampled global tree and Madagascar tree
png("../results/randsampled_global_madagascar_ltt.png", 
    width = 650, height = 650)
ltt.plot(keep.tip(ultra_5709taxa_malaise_tree$path_tree, 
                  tip = sample(x = 1:Ntip(ultra_5709taxa_malaise_tree$path_tree),
                               size = Ntip(ultra_Madagascar_tree), replace = F)), 
         backward=TRUE, log="y", 
         ylab="log(Number of lineages)",
         xlab="Time before present (Millions of years)", 
         col="black", main=paste("Randomly sampled global tree and Madagascar tree \n(",
                                 Ntip(ultra_Madagascar_tree)," taxa)", sep = ""))
ltt.lines(ultra_Madagascar_tree, backward=TRUE, col="darkgreen")
dev.off()
