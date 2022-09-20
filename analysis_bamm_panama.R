rm(list = ls())
library(ape)
library(BAMMtools)
panama_tree_extd <- 
  read.tree("ultra_aa_mitogenome_Panama_tree_extd_nozero.tre")
panama_edata <- getEventData(panama_tree_extd, 
                             eventdata = "bamm_aa_5682mito_Panama_event_data.txt", 
                             burnin=0.1)
st <- max(branching.times(panama_tree_extd))

pdf("4000sample_bamm_panama_speciation.pdf")
plotRateThroughTime(panama_edata, ratetype="speciation",
                    intervalCol="blue", avgCol="blue", 
                    intervals=c(0.05, 0.95), opacity = 0.075,
                    start.time=st, ylim=c(0,0.6), cex.axis=1, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.1, padj = 1)
mtext(side=2, line=3, "Speciation rate", font=2, cex=1.1)
dev.off()

pdf("4000sample_bamm_panama_extinction.pdf")
plotRateThroughTime(panama_edata, ratetype="extinction",
                    intervalCol="red", avgCol="red", 
                    intervals=c(0.05, 0.95), opacity = 0.075,
                    start.time=st, ylim=c(0,0.6), cex.axis=1, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.1, padj = 1)
mtext(side=2, line=3, "Extinction rate", font=2, cex=1.1)
dev.off()

pdf("4000sample_bamm_panama_netdiv.pdf")
plotRateThroughTime(panama_edata, ratetype="netdiv",
                    intervalCol="darkgreen", avgCol="darkgreen", 
                    intervals=c(0.05, 0.95), opacity = 0.08,
                    start.time=st, ylim=c(-0.1,0.15), cex.axis=1, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present (Millions of years)", 
      font=2, cex=1.1, padj = 1)
mtext(side=2, line=3, "Net diversification rate", font=2, cex=1.1)
dev.off()


pdf("4000sample_bamm_panama_ratethroughtime.pdf",
    width = 9, height = 16)
par(mfrow=c(3,1))
plotRateThroughTime(panama_edata, ratetype="speciation",
                    intervalCol="blue", avgCol="blue", 
                    intervals=c(0.05, 0.95), opacity = 0.075,
                    start.time=st, ylim=c(0,0.6), cex.axis=1.3, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Speciation rate", font=2, cex=1.4, padj = -0.5)
plotRateThroughTime(panama_edata, ratetype="extinction",
                    intervalCol="red", avgCol="red", 
                    intervals=c(0.05, 0.95), opacity = 0.075,
                    start.time=st, ylim=c(0,0.6), cex.axis=1.3, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Extinction rate", font=2, cex=1.4, padj = -0.5)
plotRateThroughTime(panama_edata, ratetype="netdiv",
                    intervalCol="darkgreen", avgCol="darkgreen", 
                    intervals=c(0.05, 0.95), opacity = 0.08,
                    start.time=st, ylim=c(-0.1,0.15), cex.axis=1.3, useMedian = F,
                    axis.labels = F)
mtext(side=1, line=2, "Time from the present\n(Millions of years)", 
      font=2, cex=1.4, padj = 1)
mtext(side=2, line=3, "Net diversification rate", font=2, cex=1.4, padj = -0.5)
dev.off()


rtt_panama <- getRateThroughTimeMatrix(panama_edata)
rm(panama_edata)
save.image("4000sample_bamm_panama.RData")