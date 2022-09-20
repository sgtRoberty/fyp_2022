rm(list = ls())
library(ape)
library(BAMMtools)

malaysia_tree_extd <- 
  read.tree("../aa_5682mito_Malaysia/ultra_aa_mitogenome_Malaysia_tree_extd_nozero.tre")
malaysia_edata <- getEventData(malaysia_tree_extd, 
                               eventdata = "../aa_5682mito_Malaysia/bamm_aa_5682mito_Malaysia_event_data.txt", 
                               burnin=0.1)
st_malaysia <- max(branching.times(malaysia_tree_extd))

panama_tree_extd <- 
  read.tree("../aa_5682mito_Panama/ultra_aa_mitogenome_Panama_tree_extd_nozero.tre")
panama_edata <- getEventData(panama_tree_extd, 
                             eventdata = "../aa_5682mito_Panama/bamm_aa_5682mito_Panama_event_data.txt", 
                             burnin=0.1)
st_panama <- max(branching.times(panama_tree_extd))


pdf("bamm_tworegions_speciation.pdf",
    width = 7, height = 14)
par(mfrow=c(2,1))
plotRateThroughTime(malaysia_edata, ratetype="speciation",
                    intervalCol="darkgreen", avgCol="forestgreen", 
                    intervals=c(0.05, 0.95), opacity = 0.095,
                    start.time=st_malaysia, ylim=c(0,0.6), cex.axis=1.3, 
                    useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Speciation rate", font=2, cex=1.4, padj = -0.5)
text(x=200, y= 0.35, label="Malaysia", font=2, cex=1.5, pos=4, col="forestgreen")

plotRateThroughTime(panama_edata, ratetype="speciation",
                    intervalCol="orange", avgCol="orange", 
                    intervals=c(0.05, 0.95), opacity = 0.12,
                    start.time=st_panama, ylim=c(0,0.6), cex.axis=1.3, 
                    useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Speciation rate", font=2, cex=1.4, padj = -0.5)
text(x=200, y= 0.55, label="Panama", font=2, cex=1.5, pos=4, col="orange")
dev.off()


pdf("bamm_tworegions_extinction.pdf",
    width = 7, height = 14)
par(mfrow=c(2,1))
plotRateThroughTime(malaysia_edata, ratetype="extinction",
                    intervalCol="darkgreen", avgCol="forestgreen", 
                    intervals=c(0.05, 0.95), opacity = 0.095,
                    start.time=st_malaysia, ylim=c(0,0.6), cex.axis=1.3, 
                    useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Extinction rate", font=2, cex=1.4, padj = -0.5)
text(x=200, y= 0.35, label="Malaysia", font=2, cex=1.5, pos=4, col="forestgreen")

plotRateThroughTime(panama_edata, ratetype="extinction",
                    intervalCol="orange", avgCol="orange", 
                    intervals=c(0.05, 0.95), opacity = 0.12,
                    start.time=st_panama, ylim=c(0,0.6), cex.axis=1.3, 
                    useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Extinction rate", font=2, cex=1.4, padj = -0.5)
text(x=200, y= 0.55, label="Panama", font=2, cex=1.5, pos=4, col="orange")
dev.off()


pdf("bamm_tworegions_netdiv.pdf",
    width = 7, height = 14)
par(mfrow=c(2,1))
plotRateThroughTime(malaysia_edata, ratetype="netdiv",
                    intervalCol="darkgreen", avgCol="forestgreen", 
                    intervals=c(0.05, 0.95), opacity = 0.095,
                    start.time=st_malaysia, ylim=c(-0.1,0.15), cex.axis=1.3, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Net diversification rate", font=2, cex=1.4, padj = -0.5)
text(x=200, y= 0.07, label="Malaysia", font=2, cex=1.5, pos=4, col="forestgreen")

plotRateThroughTime(panama_edata, ratetype="netdiv",
                    intervalCol="orange", avgCol="orange", 
                    intervals=c(0.05, 0.95), opacity = 0.12,
                    start.time=st_panama, ylim=c(-0.1,0.15), cex.axis=1.3, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Net diversification rate", font=2, cex=1.4, padj = -0.5)
text(x=200, y= 0.07, label="Panama", font=2, cex=1.5, pos=4, col="orange")
dev.off()
