# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

## Load objects
load("outs/05.dds.se.rld.vsd.RData")
rm(ebg,ct)

## load packages 
#install.packages("BiocManager")
library(BiocManager)
#install(pkgs = "DESeq2")
library("DESeq2")

# ---------------------------------------------
# starlogs ----
source("https://raw.githubusercontent.com/JakeSaunders/STARlogFunctions/master/STARlogFunctions.R")

logs <- readSTARlogsFinalOut(logs = "aligned/")

#pdf("plots/05_2.STARreport.pdf",width = 8.5,height = 8.5)
makeSTARreport(final.logs = logs,samples.lab = list.files(path = "aligned/"))
#dev.off()
par(mfrow=c(1,1))
write.csv(x = logs,file = "outs/05_2.STARlogs.csv")


# ---------------------------------------------
# col sums ----

pdf("plots/05_2.SumCountedReads.pdf",width = 8,height = 8)
barplot( 
  height = colSums(assay(dds))/ 1e6,
  names.arg = dds$short.name,
  ylim=c(0,20),
  las=2,
  col=c(rep(x = "darkorchid", 5),rep(x = "DarkGreen", 5)),
  main = "Number of Counted Reads",
  ylab = "Millions of Counted Reads",
  xlab = "Sample"
)
legend(
  x = "topleft",
  legend = c("Fed", "Starved"),
  fill = c("darkorchid","DarkGreen"),
  border = c("darkorchid","DarkGreen"),bty = "n"
)
dev.off()

pdf(file = "plots/05_2.countsFedVsStarved.pdf")
par(mfrow=c(1,2))
library(scales)
boxplot((colSums(assay(dds))/ 1e6)[1:5], notch = T,
        col = "gold", main="Fed", ylim=c(2,18),frame=F,ylab = "Number of counts")
stripchart((colSums(assay(dds))/ 1e6)[1:5],col = scales::alpha("blue",0.5),pch=19,method = "jitter",
           vertical = T,add=T,cex=3)
boxplot((colSums(assay(dds))/ 1e6)[6:10], notch = T,
        col = "gold", main="Starved",ylim=c(2,18),frame = F,ylab = "Number of counts")
stripchart((colSums(assay(dds))/ 1e6)[6:10],col = scales::alpha("blue",0.5),pch=19,method = "jitter",
           vertical = T,add=T,cex=3)
dev.off()


        
(colSums(assay(dds))/ 1e6)[6:10]





# ---------------------------------------------
# notes ----
## this is an extra QC step that I added - CJS