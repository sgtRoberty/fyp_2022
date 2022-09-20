library(ape)
library(BAMMtools)
allmito_4924taxa_tree_extd <- 
  read.tree("ultra_aa_mitogenome_all_tree_extd_nozero.tre")
allmito_edata <- getEventData(allmito_4924taxa_tree_extd, 
                              eventdata = "bamm_aa_5682mito_all_event_data.txt", 
                              burnin=0.1,
                              verbose = T)
pdf("4500sample_bamm_all.pdf")
plot.bammdata(allmito_edata, lwd=2,
              pal="temperature",
              legend=T)
addBAMMshifts(allmito_edata, cex=2)
dev.off()

pdf("4500sample_bamm_all_speciation.pdf")
plotRateThroughTime(allmito_edata, ratetype="speciation")
dev.off()

pdf("4500sample_bamm_all_extinction.pdf")
plotRateThroughTime(allmito_edata, ratetype="extinction")
dev.off()

pdf("4500sample_bamm_all_netdiv.pdf")
plotRateThroughTime(allmito_edata, ratetype="netdiv")
dev.off()

rtt <- getRateThroughTimeMatrix(allmito_edata)

rm(allmito_edata)
save.image("4500sample_bamm_all.RData")

# Load RData
load("../results/analysis_bamm_all/4500sample_bamm_all.RData")
plot(colMeans(rtt$lambda) ~ rtt$times)
plot(colMeans(rtt$mu) ~ rtt$times)
plot((colMeans(rtt$lambda) - colMeans(rtt$mu)) ~ rtt$times)
mean(rtt$lambda)
mean(rtt$mu)
mean(rtt$lambda) - mean(rtt$mu)
