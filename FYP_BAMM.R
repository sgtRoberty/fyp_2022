##### BAMMtools #####
install.packages("BAMMtools")
library(BAMMtools)
source("../code/FYP/bammCheck.R")

##### Whales example analysis #####
# https://lukejharmon.github.io/ilhabela/instruction/2015/07/02/diversification-analysis-bamm-rpanda/
# http://bamm-project.org/bammgraph.html
data(whales, events.whales)
plot.phylo(whales, cex = 0.5, font = 4)
edata_whales <- BAMMtools::getEventData(whales, 
                                        events.whales, 
                                        burnin=0.1)
plot.bammdata(edata_whales, lwd=3, 
              pal="temperature",
              labels = T, cex = 0.5)

?subtreeBAMM
whales$tip.label[substr(whales$tip.label, 1, 10) == "Mesoplodon"]
whales$tip.label[substr(whales$tip.label, 1, 9) == "Delphinus"]
subset_edata_whales <- subtreeBAMM(edata_whales,
                                   tips = c(whales$tip.label[substr(whales$tip.label, 
                                                                  1, 10) == "Mesoplodon"],
                                            whales$tip.label[substr(whales$tip.label, 
                                                                    1, 9) == "Delphinus"]))
plot.bammdata(subset_edata_whales, lwd=3, 
              method="polar", pal="temperature",
              labels = T, cex = 0.6)

# Rate through time plots
plot.new()

pdf("../results/FYP_figs/whale_ratethroughtime.pdf",
    width = 17, height = 7)
par(mfrow=c(1,3))
st <- max(branching.times(whales))
plotRateThroughTime(edata_whales, intervalCol="red", avgCol="red", 
                    start.time=st, ylim=c(0,1), cex.axis=2)
text(x=30, y= 0.2, label="All whales", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, intervalCol="blue", avgCol="blue", 
                    start.time=st, node=140, ylim=c(0,1), xlim=c(35,0), cex.axis=1.5)
text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, intervalCol="darkgreen", avgCol="darkgreen", 
                    start.time=st, node=140, nodetype = "exclude", ylim=c(0,1), cex.axis=1.5)
text(x=30, y= 0.8, label="Non-dolphins", font=4, cex=2.0, pos=4)
dev.off()

# Evolutionary rate variation through time: grayscale
plot.new()
par(mfrow=c(1,3))
st <- max(branching.times(whales))

plotRateThroughTime(edata_whales, avgCol="black", start.time=st, 
                    ylim=c(0,1), cex.axis=1, intervalCol='gray80', 
                    intervals=c(0.05, 0.95), opacity=1, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.1, padj = 1)
mtext(side=2, line=3, "Speciation rate", font=2, cex=1.1, padj = -0.5)

text(x=30, y= 0.8, label="All whales", font=4, cex=2.0, pos=4)
plotRateThroughTime(edata_whales, avgCol="black", start.time=st, 
                    node=140, ylim=c(0,1), xlim=c(35,0), cex.axis=1.5,intervalCol='gray80', 
                    intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Dolphins only", font=4, cex=2.0, pos=4)

plotRateThroughTime(edata_whales, avgCol="black", start.time=st, node=140, 
                    nodetype = "exclude", ylim=c(0,1), cex.axis=1.5, intervalCol='gray80', 
                    intervals=c(0.05, 0.95), opacity=1)
text(x=30, y= 0.8, label="Non-dolphins", font=4, cex=2.0, pos=4)

# Extracting tip rates
tip.rates <- getTipRates(edata_whales)
str(tip.rates)
hist(tip.rates$lambda.avg,xlab="average lambda")
hist(tip.rates$mu.avg,xlab="average mu")

boxplot(tip.rates$lambda.avg[53:87], tip.rates$lambda.avg[1:52], 
        col = c("red", "blue"), 
        names = c("dolphins", "other cetaceans"))

# Marginal shift probabilities
marg_probs <- marginalShiftProbsTree(edata_whales)
plot.phylo(marg_probs, use.edge.length = F,
           no.margin = T, show.tip.label = F)

# Distinct shift configurations and their frequencies
cset <- credibleShiftSet(edata_whales, expectedNumberOfShifts=1, threshold=3)
plot.credibleshiftset(cset, lwd=2.5)

# Bayesian credible sets of shift configurations
css <- credibleShiftSet(edata_whales, 
                        expectedNumberOfShifts=1, 
                        threshold=5, set.limit = 0.95)
css$number.distinct # Number of distinct shift configurations in the data
summary(css)
plot.credibleshiftset(css)


##### ultra_allmito_5682taxa_tree$path_tree #####
##### allmito family taxon count #####
allmito_5682taxa_family_table <- as.data.frame(table(aa_supermatrix_metadata$family))
colnames(allmito_5682taxa_family_table) <- c("family", "5682taxa_count")
write_csv(allmito_5682taxa_family_table, 
          file = "../results/FYP_family_table/allmito_5682taxa_family_table.csv")

# Manually add total global count to allmito_5682taxa_family_table_globalcount.csv"
allmito_5682taxa_global_family_table <- 
  read_csv("../results/FYP_family_table/allmito_5682taxa_family_table_globalcount.csv")
for (i in 1:length(allmito_5682taxa_global_family_table$alt_family)){
  if (is.na(allmito_5682taxa_global_family_table$alt_family[i])){
    allmito_5682taxa_global_family_table$alt_family[i] <- 
      allmito_5682taxa_global_family_table$family[i]
  }
}

library(tidyverse)
allmito_5682taxa_global_alt_family_table <- 
  allmito_5682taxa_global_family_table[!names(allmito_5682taxa_global_family_table) 
                                       == "family"]
allmito_5682taxa_global_alt_family_table <-
  allmito_5682taxa_global_alt_family_table %>% 
  group_by(alt_family) %>% 
  summarize_all(sum, na.rm = T)

sum(allmito_5682taxa_global_alt_family_table$global_count)
sum(allmito_5682taxa_global_alt_family_table$`5682taxa_count`)

allmito_5682taxa_global_alt_family_table$Sample_Probs <- 
  allmito_5682taxa_global_alt_family_table$`5682taxa_count` / allmito_5682taxa_global_alt_family_table$global_count

allmito_5682taxa_global_alt_family_table <- 
allmito_5682taxa_global_alt_family_table %>%
  filter(!alt_family %in% c("Megaloptera", "Neuroptera",
                            "Raphidioptera", "Strepsiptera"))
write_csv(allmito_5682taxa_global_alt_family_table, 
          "../results/FYP_family_table/allmito_5682taxa_global_alt_family_table.csv")

##### allmito sample probabilities ######
# Read in allmito_5682taxa_global_alt_family_table.csv
allmito_5682taxa_global_alt_family_table <- 
  read_csv("../results/FYP_family_table/allmito_5682taxa_global_alt_family_table.csv")
# # Malaysia sample probabilities
allmito_sampleProbs <- data.frame(speciesName=aa_supermatrix_metadata$db_id, 
                                   cladeName=aa_supermatrix_metadata$family)
allmito_sampleProbs <- drop_na(allmito_sampleProbs)
allmito_sampleProbs <- 
  filter(allmito_sampleProbs, !cladeName %in% c("Chrysopidae",
                                                "Coniopterygidae",
                                                "Corydalidae",
                                                "Inocelliidae",
                                                "Myrmeleontidae",
                                                "Osmylidae",
                                                "Raphidiidae",
                                                "Stylopidae",
                                                "Xenidae"))
allmito_sampleProbs$cladeName <-
  str_replace_all(allmito_sampleProbs$cladeName, 
                  c("Anobiidae" = "Ptinidae", "Apionidae" = "Brentidae",
                    "Bolboceratidae" = "Geotrupidae", "Brachyceridae" = "Curculionidae",
                    "Cicindelidae" = "Carabidae", 
                    "Cybocephalidae" = "Nitidulidae", "Dryophthoridae" = "Curculionidae",
                    "Helophoridae" = "Hydrophilidae", "Hydrochidae" = "Hydrophilidae",
                    "Platypodidae" = "Curculionidae", "Scydmaenidae" = "Staphylinidae",
                    "Spercheidae" = "Hydrophilidae"))

length(unique(allmito_sampleProbs$cladeName))

for (i in 1:length(allmito_sampleProbs$cladeName)){
  Name <- allmito_sampleProbs$cladeName[i]
  allmito_sampleProbs$samplingFraction[i] <- 
    filter(allmito_5682taxa_global_alt_family_table, alt_family == Name)$Sample_Probs
}

allmito_sampleProbs <- rbind(
  c(as.character(length(allmito_sampleProbs$speciesName)/sum(allmito_5682taxa_global_alt_family_table$global_count)), 
    "", ""), 
  allmito_sampleProbs)
# Write out sample probability files
write_tsv(allmito_sampleProbs,
          "../results/bamm_aa_5682mito_all/allmito_sampleProbs.txt",
          col_names = F)

##### allmito subset tree for BAMM and set priors #####
# read in allmito_sampleProbs
allmito_sampleProbs <- read_tsv("../results/bamm_aa_5682mito_all/allmito_sampleProbs.txt",
                                col_names = F)
colnames(allmito_sampleProbs) <- c("speciesName", "cladeName", "samplingFraction")
# Force ultrametric by extend
allmito_4924taxa_tree_extd <- ultra_allmito_5682taxa_tree$path_tree
allmito_4924taxa_tree_extd <- keep.tip(allmito_4924taxa_tree_extd, 
                                       tip = allmito_sampleProbs$speciesName[-1])
Ntip(allmito_4924taxa_tree_extd) # 4924 taxa, excluding non-Coleoptera, all family level taxonomised
allmito_4924taxa_tree_extd$edge.length[allmito_4924taxa_tree_extd$edge.length == 0] <- 1e-8
head(table(allmito_4924taxa_tree_extd$edge.length), 20)
allmito_4924taxa_tree_extd <- 
  phytools::force.ultrametric(allmito_4924taxa_tree_extd,
                              method = "extend")
# Check the tree for binary, ultrametric, zero branch length
is.binary(allmito_4924taxa_tree_extd)
is.ultrametric(allmito_4924taxa_tree_extd)
max(allmito_4924taxa_tree_extd$edge.length)
min(allmito_4924taxa_tree_extd$edge.length)
head(table(allmito_4924taxa_tree_extd$edge.length), 20)

# Write out 4924-taxon, extended, no-zero-branch-length tree
write.tree(allmito_4924taxa_tree_extd, 
           "../results/bamm_aa_5682mito_all/ultra_aa_mitogenome_all_tree_extd_nozero.tre")

# Set priors
allmito_5682taxa_tree_extd_priors <- setBAMMpriors(allmito_4924taxa_tree_extd, 
                                                   outfile = NULL)
allmito_5682taxa_tree_extd_priors

generateControlFile('../results/bamm_aa_5682mito_all/aa_5682mito_all_diversification.txt', 
                    type = 'diversification', 
                    params = list(
                      treefile = 'ultra_aa_mitogenome_all_tree_extd_nozero.tre',
                      #loadEventData = 1,
                      #eventDataInfile = "bamm_aa_5682mito_all_event_data_in.txt",
                      seed = 12345,
                      useGlobalSamplingProbability = 0,
                      sampleProbsFilename = "allmito_sampleProbs.txt",
                      overwrite = 1,
                      expectedNumberOfShifts = "65",
                      lambdaInitPrior = as.numeric(allmito_5682taxa_tree_extd_priors['lambdaInitPrior']),
                      lambdaShiftPrior = as.numeric(allmito_5682taxa_tree_extd_priors['lambdaShiftPrior']),
                      muInitPrior = as.numeric(allmito_5682taxa_tree_extd_priors['muInitPrior']),
                      numberOfGenerations = "4500000",
                      mcmcWriteFreq = "1000",
                      eventDataWriteFreq = "1000",
                      printFreq = "1000",
                      acceptanceResetFreq = "1000",
                      deltaT = 0.1,
                      outName = "bamm_aa_5682mito_all",
                      minCladeSizeForShift = 10
                    )
)

##### allmito BAMM output #####
allmito_4924taxa_tree_extd <- 
  read.tree("../results/bamm_aa_5682mito_all/ultra_aa_mitogenome_all_tree_extd_nozero.tre")
allmito_edata <- getEventData(allmito_4924taxa_tree_extd, 
                              eventdata = "../results/bamm_aa_5682mito_all/bamm_aa_5682mito_all_event_data.txt", 
                              burnin=0.1,
                              nsamples = 2500,
                              verbose = T)
malaysia_tree_extd <- 
  read.tree("../results/bamm_aa_5682mito_Malaysia/ultra_aa_mitogenome_Malaysia_tree_extd_nozero.tre")
panama_tree_extd <- 
  read.tree("../results/bamm_aa_5682mito_Panama/ultra_aa_mitogenome_Panama_tree_extd_nozero.tre")
subset_malaysia_edata <- subtreeBAMM(allmito_edata,
                                     tips = malaysia_tree_extd$tip.label)
subset_panama_edata <- subtreeBAMM(allmito_edata,
                                     tips = panama_tree_extd$tip.label)
# Mean phylorate plot
plot.bammdata(allmito_edata, lwd=2,
              pal="temperature",
              legend=T)
addBAMMshifts(allmito_edata, cex=2)

# Assessing MCMC convergence
mcmcout <- read.csv("../results/bamm_aa_5682mito_all/bamm_aa_5682mito_all_mcmc_out.txt", 
                    header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
geweke.diag(postburn)

# bammCheck.R
bammCheck(expectedNumberOfShifts = 50, burnin = 0.1,
          mcmcFile = "../results/bamm_aa_5682mito_all/bamm_aa_5682mito_all_mcmc_out.txt",
          chainFile = "../results/bamm_aa_5682mito_all/bamm_aa_5682mito_all_chain_swap.txt") 

# Analysis of rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
# The posterior probabilities of each rate shift count observed during 
# simulation of the posterior
shift_probs <- summary(allmito_edata)
# Posterior probabilities are dependent on the prior on the number of shifts, unlike Bayes factors

# Bayes factors for model comparison - Prior distribution in BAMM
bfmat <- computeBayesFactors(mcmcout, 
                             expectedNumberOfShifts=50, 
                             burnin=0.1)
bfmat
plotPrior(mcmcout, expectedNumberOfShifts=50)

# Bayesian credible sets of shift configurations
css <- credibleShiftSet(allmito_edata, 
                        expectedNumberOfShifts=50, 
                        threshold=5, set.limit = 0.95)
css$number.distinct # Number of distinct shift configurations in the data
summary(css)
plot.credibleshiftset(css)

# Finding the single best shift configuration
best <- getBestShiftConfiguration(allmito_edata, 
                                  expectedNumberOfShifts=35)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex = 1)

# Marginal shift probabilities
marg_probs <- marginalShiftProbsTree(allmito_edata)
plot.phylo(marg_probs, show.tip.label = F)

# Rate-through-time analysis
plotRateThroughTime(allmito_edata, ratetype="speciation")
plotRateThroughTime(allmito_edata, ratetype="extinction")
plotRateThroughTime(allmito_edata, ratetype="netdiv")

rtt <- getRateThroughTimeMatrix(allmito_edata)
plot(colMeans(rtt$lambda) ~ rtt$times)
plot(colMeans(rtt$mu) ~ rtt$times)
plot((colMeans(rtt$lambda) - colMeans(rtt$mu)) ~ rtt$times)

mean(rtt$lambda)
mean(rtt$mu)
mean(rtt$lambda) - mean(rtt$mu)

rtt_subset_malaysia <- getRateThroughTimeMatrix(subset_malaysia_edata)
plot(colMeans(rtt_subset_malaysia$lambda) ~ rtt_subset_malaysia$times)
plot(colMeans(rtt_subset_malaysia$mu) ~ rtt_subset_malaysia$times)
plot((colMeans(rtt_subset_malaysia$lambda) - colMeans(rtt_subset_malaysia$mu)) ~ rtt_subset_malaysia$times)
plotRateThroughTime(subset_malaysia_edata, ratetype="netdiv")

rtt_subset_panama <- getRateThroughTimeMatrix(subset_panama_edata)
plot(colMeans(rtt_subset_panama$lambda) ~ rtt_subset_panama$times)
plot(colMeans(rtt_subset_panama$mu) ~ rtt_subset_panama$times)
plot((colMeans(rtt_subset_panama$lambda) - colMeans(rtt_subset_panama$mu)) ~ rtt_subset_panama$times)
plotRateThroughTime(subset_panama_edata, ratetype="netdiv")

##### Incomplete, non-random sampling #####
##### Malaysia and Panama family taxon count# ####
Malaysia_metadata <- filter(aa_supermatrix_metadata, country == "Malaysia")
Panama_metadata <- filter(aa_supermatrix_metadata, country == "Panama")

Malaysia_family_table <- as.data.frame(table(Malaysia_metadata$family))
rownames(Malaysia_family_table) <- Malaysia_family_table$Var1
Malaysia_family_table$Var1 <- NULL
Panama_family_table <- as.data.frame(table(Panama_metadata$family))
rownames(Panama_family_table) <- Panama_family_table$Var1
Panama_family_table$Var1 <- NULL

Malaysia_Panama_family_table <- merge(Malaysia_family_table, Panama_family_table,
                                      by = 'row.names', all = T)
colnames(Malaysia_Panama_family_table) <- c("family", "Malaysia_count", "Panama_count")
write_csv(Malaysia_Panama_family_table, 
          file = "../results/FYP_family_table/Malaysia_Panama_family_table.csv")

# Manually add total global count to Malaysia_Panama_family_table_globalcount.csv
global_family_table <- read_csv("../results/FYP_family_table/Malaysia_Panama_family_table_globalcount.csv")
for (i in 1:length(global_family_table$alt_family)){
  if (is.na(global_family_table$alt_family[i])){
    global_family_table$alt_family[i] <- global_family_table$family[i]
  }
}

library(tidyverse)
global_alt_family_table <- global_family_table[!names(global_family_table) 
                                               == "family"]
global_alt_family_table <-
  global_alt_family_table %>% 
  group_by(alt_family) %>% 
  summarize_all(sum, na.rm = T)

global_alt_family_table$Malaysia_sampleProbs <- 
  global_alt_family_table$Malaysia_count / global_alt_family_table$global_count
global_alt_family_table$Panama_sampleProbs <- 
  global_alt_family_table$Panama_count / global_alt_family_table$global_count

View(global_alt_family_table)
write_csv(global_alt_family_table, 
          "../results/FYP_family_table/global_alt_family_table.csv")

# Read in global_alt_family_table.csv
global_alt_family_table <- read_csv("../results/FYP_family_table/global_alt_family_table.csv")
# Malaysia sample probabilities
Malaysia_sampleProbs <- data.frame(speciesName=Malaysia_metadata$db_id, 
                                   cladeName=Malaysia_metadata$family)
Malaysia_sampleProbs <- drop_na(Malaysia_sampleProbs)
Malaysia_sampleProbs$cladeName <-
  str_replace_all(Malaysia_sampleProbs$cladeName, 
                  c("Anobiidae" = "Ptinidae", "Apionidae" = "Brentidae",
                    "Bolboceratidae" = "Geotrupidae", "Brachyceridae" = "Curculionidae",
                    "Cicindelidae" = "Carabidae", "Cybocephalidae" = "Nitidulidae", 
                    "Dryophthoridae" = "Curculionidae","Platypodidae" = "Curculionidae"))
for (i in 1:length(Malaysia_sampleProbs$cladeName)){
  Name <- Malaysia_sampleProbs$cladeName[i]
  Malaysia_sampleProbs$samplingFraction[i] <- 
    filter(global_alt_family_table, alt_family == Name)$Malaysia_sampleProbs
}

Malaysia_sampleProbs <- rbind(c(as.character(length(Malaysia_sampleProbs$speciesName)/386500), 
                                "", ""), 
                              Malaysia_sampleProbs)

# Panama sample probabilities
Panama_sampleProbs <- data.frame(speciesName=Panama_metadata$db_id, 
                                 cladeName=Panama_metadata$family)
Panama_sampleProbs <- drop_na(Panama_sampleProbs)
Panama_sampleProbs$cladeName <-
  str_replace_all(Panama_sampleProbs$cladeName, 
                  c("Anobiidae" = "Ptinidae", "Apionidae" = "Brentidae",
                    "Bolboceratidae" = "Geotrupidae", "Brachyceridae" = "Curculionidae",
                    "Cicindelidae" = "Carabidae", "Cybocephalidae" = "Nitidulidae",
                    "Dryophthoridae" = "Curculionidae", "Platypodidae" = "Curculionidae"))
for (i in 1:length(Panama_sampleProbs$cladeName)){
  Name <- Panama_sampleProbs$cladeName[i]
  Panama_sampleProbs$samplingFraction[i] <- 
    filter(global_alt_family_table, alt_family == Name)$Panama_sampleProbs
}

Panama_sampleProbs <- rbind(c(as.character(length(Panama_sampleProbs$speciesName)/386500), 
                              "", ""), 
                            Panama_sampleProbs)

# Write out sample probability files
write_tsv(Malaysia_sampleProbs, 
          "../results/bamm_aa_5682mito_Malaysia/Malaysia_sampleProbs.txt",
          col_names = F)
write_tsv(Panama_sampleProbs, 
          "../results/bamm_aa_5682mito_Panama/Panama_sampleProbs.txt",
          col_names = F)

##### Set BAMM priors #####
# http://bamm-project.org/configuration.html
# Force ultrametric by extend
# Malaysia tree
Malaysia_sampleProbs <- read_tsv("../results/bamm_aa_5682mito_Malaysia/Malaysia_sampleProbs.txt",
                                 col_names = F)
colnames(Malaysia_sampleProbs) <- c("speciesName", "cladeName", "samplingFraction")

malaysia_tree_extd <- ultra_aa_mitogenome_Malaysia_tree
malaysia_tree_extd <- keep.tip(malaysia_tree_extd, 
                               tip = Malaysia_sampleProbs$speciesName[-1])
Ntip(malaysia_tree_extd)
malaysia_tree_extd$edge.length[malaysia_tree_extd$edge.length == 0] <- 1e-8
head(table(malaysia_tree_extd$edge.length), 20)
malaysia_tree_extd <- phytools::force.ultrametric(malaysia_tree_extd,
                                                  method = "extend")
is.binary(malaysia_tree_extd)
is.ultrametric(malaysia_tree_extd)
max(malaysia_tree_extd$edge.length)
min(malaysia_tree_extd$edge.length)
head(table(malaysia_tree_extd$edge.length), 20)

write.tree(malaysia_tree_extd, 
           "../results/bamm_aa_5682mito_Malaysia/ultra_aa_mitogenome_Malaysia_tree_extd_nozero.tre")
malaysia_tree_extd_priors <- setBAMMpriors(malaysia_tree_extd, outfile = NULL)
malaysia_tree_extd_priors

generateControlFile('../results/bamm_aa_5682mito_Malaysia/aa_5682mito_Malaysia_diversification.txt', 
                    type = 'diversification', 
                    params = list(
                      treefile = 'ultra_aa_mitogenome_Malaysia_tree_extd_nozero.tre',
                      #loadEventData = 1,
                      #eventDataInfile = "bamm_aa_5682mito_Malaysia_event_data_in.txt",
                      seed = 12345,
                      useGlobalSamplingProbability = 0,
                      sampleProbsFilename = "Malaysia_sampleProbs.txt",
                      overwrite = 1,
                      expectedNumberOfShifts = "20",
                      lambdaInitPrior = as.numeric(malaysia_tree_extd_priors['lambdaInitPrior']),
                      lambdaShiftPrior = as.numeric(malaysia_tree_extd_priors['lambdaShiftPrior']),
                      muInitPrior = as.numeric(malaysia_tree_extd_priors['muInitPrior']),
                      numberOfGenerations = "40 000 000",
                      mcmcWriteFreq = "10000",
                      eventDataWriteFreq = "10000",
                      printFreq = "5000",
                      acceptanceResetFreq = "10000",
                      outName = "bamm_aa_5682mito_Malaysia",
                      numberOfChains = 6,
                      minCladeSizeForShift = 10
                    )
)

# Panama tree
Panama_sampleProbs <- read_tsv("../results/bamm_aa_5682mito_Panama/Panama_sampleProbs.txt",
                               col_names = F)
colnames(Panama_sampleProbs) <- c("speciesName", "cladeName", "samplingFraction")

panama_tree_extd <- ultra_aa_mitogenome_Panama_tree
panama_tree_extd <- keep.tip(panama_tree_extd, 
                             tip = Panama_sampleProbs$speciesName[-1])
Ntip(panama_tree_extd)
panama_tree_extd$edge.length[panama_tree_extd$edge.length == 0] <- 1e-8
head(table(panama_tree_extd$edge.length), 20)
panama_tree_extd <- phytools::force.ultrametric(panama_tree_extd,
                                                method = "extend")
is.binary(panama_tree_extd)
is.ultrametric(panama_tree_extd)
max(panama_tree_extd$edge.length)
min(panama_tree_extd$edge.length)
head(table(panama_tree_extd$edge.length), 20)

write.tree(panama_tree_extd, "../results/bamm_aa_5682mito_Panama/ultra_aa_mitogenome_Panama_tree_extd_nozero.tre")

panama_tree_extd_priors <- setBAMMpriors(panama_tree_extd, outfile = NULL)
panama_tree_extd_priors

generateControlFile('../results/bamm_aa_5682mito_Panama/aa_5682mito_Panama_diversification.txt', 
                    type = 'diversification', 
                    params = list(
                      treefile = 'ultra_aa_mitogenome_Panama_tree_extd_nozero.tre',
                      #loadEventData = 1,
                      #eventDataInfile = "bamm_aa_5682mito_Panama_event_data_in.txt",
                      seed = 12345,
                      useGlobalSamplingProbability = 0,
                      sampleProbsFilename = "Panama_sampleProbs.txt",
                      overwrite = 1,
                      expectedNumberOfShifts = "10",
                      lambdaInitPrior = as.numeric(panama_tree_extd_priors['lambdaInitPrior']),
                      lambdaShiftPrior = as.numeric(panama_tree_extd_priors['lambdaShiftPrior']),
                      muInitPrior = as.numeric(panama_tree_extd_priors['muInitPrior']),
                      numberOfGenerations = "40000000",
                      mcmcWriteFreq = "10000",
                      eventDataWriteFreq = "10000",
                      printFreq = "5000",
                      acceptanceResetFreq = "10000",
                      outName = "bamm_aa_5682mito_Panama",
                      numberOfChains = 6,
                      minCladeSizeForShift = 10
                    )
)

##### ultra_aa_mitogenome_Malaysia_tree_extd_nozero.tre #####
malaysia_tree_extd <- 
  read.tree("../results/bamm_aa_5682mito_Malaysia/ultra_aa_mitogenome_Malaysia_tree_extd_nozero.tre")
malaysia_edata <- getEventData(malaysia_tree_extd, 
                               eventdata = "../results/bamm_aa_5682mito_Malaysia/bamm_aa_5682mito_Malaysia_event_data.txt", 
                               burnin=0.1)
# Mean phylorate plot
plot.bammdata(malaysia_edata, lwd=3, 
              pal="temperature",
              legend=T)
addBAMMshifts(malaysia_edata, cex=2)

# Assessing MCMC convergence
mcmcout <- read.csv("../results/bamm_aa_5682mito_Malaysia/bamm_aa_5682mito_Malaysia_mcmc_out.txt", 
                    header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
geweke.diag(postburn)

# bammCheck.R
bammCheck(expectedNumberOfShifts = 20, burnin = 0.1,
          mcmcFile = "../results/bamm_aa_5682mito_Malaysia/bamm_aa_5682mito_Malaysia_mcmc_out.txt",
          chainFile = "../results/bamm_aa_5682mito_Malaysia/bamm_aa_5682mito_Malaysia_chain_swap.txt") 

# Analysis of rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
# The posterior probabilities of each rate shift count observed during 
# simulation of the posterior
shift_probs <- summary(malaysia_edata)
# Posterior probabilities are dependent on the prior on the number of shifts, unlike Bayes factors

# Bayes factors for model comparison - Prior distribution in BAMM
bfmat <- computeBayesFactors(mcmcout, 
                             expectedNumberOfShifts=20, 
                             burnin=0.1)
bfmat
plotPrior(mcmcout, expectedNumberOfShifts=20)

# Bayesian credible sets of shift configurations
css <- credibleShiftSet(malaysia_edata, 
                        expectedNumberOfShifts=20, 
                        threshold=5, set.limit = 0.95)
css$number.distinct # Number of distinct shift configurations in the data
summary(css)
plot.credibleshiftset(css)

# Finding the single best shift configuration
best <- getBestShiftConfiguration(malaysia_edata, 
                                  expectedNumberOfShifts=20)
plot.bammdata(best, lwd = 2)
addBAMMshifts(best, cex = 1)

# Marginal shift probabilities
marg_probs <- marginalShiftProbsTree(malaysia_edata)
plot.phylo(marg_probs, show.tip.label = F)

# Rate-through-time analysis
plotRateThroughTime(malaysia_edata, ratetype="speciation")
plotRateThroughTime(malaysia_edata, ratetype="extinction")
plotRateThroughTime(malaysia_edata, ratetype="netdiv")

##### ultra_aa_mitogenome_Panam_tree_extd_nozero.tre #####
panama_tree_extd <- 
  read.tree("../results/bamm_aa_5682mito_Panama/ultra_aa_mitogenome_Panama_tree_extd_nozero.tre")
panama_edata <- getEventData(panama_tree_extd, 
                             eventdata = "../results/bamm_aa_5682mito_Panama/bamm_aa_5682mito_Panama_event_data.txt", 
                             burnin=0.1)
# Mean phylorate plot
plot.bammdata(panama_edata, lwd=2, 
              pal="temperature",
              legend=T)
addBAMMshifts(panama_edata, cex=2)

# Assessing MCMC convergence
mcmcout <- read.csv("../results/bamm_aa_5682mito_Panama/bamm_aa_5682mito_Panama_mcmc_out.txt", 
                    header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

# bammCheck.R
bammCheck(expectedNumberOfShifts = 10, burnin = 0.1,
          mcmcFile = "../results/bamm_aa_5682mito_Panama/bamm_aa_5682mito_Panama_mcmc_out.txt",
          chainFile = "../results/bamm_aa_5682mito_Panama/bamm_aa_5682mito_Panama_chain_swap.txt") 

# Analysis of rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)
shift_probs <- summary(panama_edata)

# Bayes factors for model comparison - Prior distribution in BAMM
bfmat <- computeBayesFactors(mcmcout, 
                             expectedNumberOfShifts=10, 
                             burnin=0.1)
bfmat
plotPrior(mcmcout, expectedNumberOfShifts=10)

# Bayesian credible sets of shift configurations
css <- credibleShiftSet(panama_edata, 
                        expectedNumberOfShifts=10, 
                        threshold=5, set.limit = 0.95)
css$number.distinct # Number of distinct shift configurations in the data
summary(css)
plot.credibleshiftset(css)

# Finding the single best shift configuration
best <- getBestShiftConfiguration(panama_edata, 
                                  expectedNumberOfShifts=10)
plot.bammdata(best, lwd = 1.5)
addBAMMshifts(best, cex = 1.5)

# Marginal shift probabilities
marg_probs <- marginalShiftProbsTree(panama_edata)
plot.phylo(marg_probs, show.tip.label = F)

# Rate-through-time analysis
plotRateThroughTime(panama_edata, ratetype="speciation")
plotRateThroughTime(panama_edata, ratetype="extinction")
plotRateThroughTime(panama_edata, ratetype="netdiv")

rtt <- getRateThroughTimeMatrix(panama_edata)
plot(colMeans(rtt$lambda) ~ rtt$times)
