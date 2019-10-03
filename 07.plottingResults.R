# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

## Load objects
load("outs/06.dds.rld.res.RData")

## load packages 
#install.packages("BiocManager")
library(BiocManager)
#install(pkgs = "DESeq2")
library("DESeq2")

# ---------------------------------------------
# 6.1 Counts plot ----
## Gene with lowest p-adjusted in base graphics
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dependent.var"),xlab = )


## Gene with ggplot
install.packages("ggbeeswarm")
library(ggbeeswarm)


padj <- sig.gene$padj[12]
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dependent.var"),
           returnData = TRUE)
ggplot(geneCounts, aes(x = dependent.var, y = count,color = rownames(geneCounts),size=4)) +
    scale_y_log10() +  geom_beeswarm(cex = 3) + 
    geom_text(aes(label=rownames(geneCounts)),hjust=0.5, vjust=0.5,col="black",position = "jitter") +
    ggtitle(topGene) + xlab(paste( "adjusted p-value =" ,round(padj,digits = 4)))

ggplot(geneCounts, aes(x = dependent.var, y = count,color = rownames(geneCounts),size=4)) +
    scale_y_log10() + geom_point(size = 3)

# ---------------------------------------------
# Make pdf report of sig genes counts plots ----

sig.gene <- read.csv(file = "outs/06.padj.less.0.01.csv")
sig.gene$X <- as.character(sig.gene$X)

pdf("plots/07.sig.exp.genes.pdf",width = 7,height = 10)
par(mfrow=c(4,2))
sig.gene <- sig.gene[order(sig.gene$padj),]
for (row in 1:nrow(sig.gene)) {
    print(
        sig.gene[row,]$padj
        )
    plotCounts(
        dds, gene = sig.gene[row,]$X, intgroup = c("dependent.var"),
        pch=18, cex = 2, col=alpha("darkgreen",alpha = 0.7),
        xlab = paste( "adjusted p-value =" ,round(sig.gene[row,]$padj,digits = 4))
    )
}
dev.off()

# ---------------------------------------------
# 6.2 MA-plot  ----
#   On the y-axis, the "M" stands for "minus" - subtraction of log values is equivalent to the log of the ratio 
#   - and on the x-axis, the "A" stands for "average".

##  load packages to use paticular shrinking method 
##  (I'm not sure why this is need, I think it just makes the data "look" better on the plot ~Jake)

resultsNames(dds)

pdf(file = "plots/07.MAplots.shrinks.pdf",width = 8.5,height = 11)
par(mfrow=c(3,1))
## example of no Shrink
res.noshr <- results(dds, name=resultsNames(dds)[2])
plotMA(res.noshr, 
       ylim = c(-20, 20),
       main= paste(resultsNames(dds)[2],": no shrink"), 
       colNonSig = "azure3",
       colSig = "black",
       colLine = ""
)

## run to apeglm Shrink, suggestted by pipeline write up
####install("apeglm")
library("apeglm")

res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
    ##  using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    ##  Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    ##  sequence count data: removing the noise and preserving large differences.
    ##  Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
plotMA(res, ylim = c(-6, 6),
       main=paste(resultsNames(dds)[2],": apeglm shrink (suggested)"), 
       colNonSig = "azure3",
       colSig = "black",
       colLine = ""
)

## run to ashr Shrink
#install.packages("ashr")
library("apeglm")
res.ashr <- lfcShrink(dds, coef=resultsNames(dds)[2], type="ashr")
plotMA(res.ashr, ylim = c(-10, 10),
       main=paste(resultsNames(dds)[2],": ashr shrink"), 
       colNonSig = "azure3",
       colSig = "black",
       colLine = ""
       )
    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##    Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ## https://doi.org/10.1093/biostatistics/kxw041
legend("bottomleft",legend = c("p > 0.1", "p < 0.1"),col = c("azure3","black"),pch=19,bty = "n",ncol=2,cex = 1)
dev.off()

#apeglm looked the best so make a single plot of this one and highlight sig genes
padj.sig <- read.csv("outs/06.padj.less.0.01.csv")

pdf(file = "plots/07.MAplots.apeglm.pdf",width = 11,height = 8.5)
par(mfrow=c(1,1))
plotMA(res, ylim = c(-6, 6),
       main="MA plot: AHK cells Fed vs Starved", 
       colNonSig = "azure3",
       colSig = "black",
       colLine = ""
)
with(res[padj.sig$X, ], {
    points(baseMean, log2FoldChange, col="dodgerblue", cex=1.5, lwd=2)
    text(baseMean, log2FoldChange, padj.sig$X, pos=2, col="black", cex=0.6)
})
with(res[rownames(res) == "CG3906", ], {
    points(baseMean, 6, col="dodgerblue", cex=1.5, lwd=2,pch = 2)
    text(baseMean, 6, padj.sig[12,]$X, pos=2, col="black", cex=0.6)
})

legend(x = "bottomleft",legend = "p < 0.1",pch = 1,col="dodgerblue",cex=1)

dev.off()

# ---------------------------------------------
#  6.3 Gene clustering ----

### let us select the 30 genes with the highest variance across samples

library(genefilter)
library(pheatmap)

pdf(file = "plots/07.heatmap.high.variance.genes.pdf", width = 8.5,height = 11)
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 40)

mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, "dependent.var"])
colnames(anno) <- "Feeding Status"
pheatmap(mat, annotation_col = anno,main = "Genes with highest variance across samples (n = 40)")
dev.off()


### heatmap of most signifcant

# load sig.gene list
sig.gene <- read.csv(file = "outs/06.padj.less.0.01.csv")

### select sig gene names in count table
# dfA[which(rownames(dfA) %in% rownames(dfB)),]
mat  <-assay(rld)[which(rownames(assay(rld)) %in% as.character(sig.gene$X) ),]

sig.gene

mat  <- assay(rld)[ topSigGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, "dependent.var"])
rownames(anno) <- c(1:10)
colnames(anno) <- "Feeding Status"
pdf(file = "plots/07.heatmap.sig.genes.pdf", width = 8.5,height = 11)
pheatmap(mat, annotation_col = anno,scale = "row",
         main = "Genes that are Signifcantly (p<0.1) Differentially Expressed")
dev.off()

# ---------------------------------------------
#  6.4 Independent filtering ----

# Consider this for the 2nd pass but skip now
# remove low rowMean gene before testing to lower FDR 
# lowly expressed gene aren't likely to be signicant anyway and can these genes have an influence 
# on the multiple testing adjustment, whose performance improves if such genes are removed.
# DESeq2() does some of this automaticly but we can tune automatic independent filtering 
# by, the results function.

# ---------------------------------------------
#  XY plot ----

#make count table object
ct <- assay(dds)

# cal mean and sd for fed vs starved
xy <- cbind.data.frame(
    names = rownames(ct),
    fed.mean = rowMeans(ct[,1:5]), 
    fed.sd = rowSds(ct[,1:5]),
    str.mean = rowMeans(ct[,6:10]), 
    str.sd = rowSds(ct[,6:10]),
    padj = res$padj,
    log2fold = res$log2FoldChange
    
)


xy[is.na(xy$padj),]$padj <- 1

### Tried in base graphs
# factorize pvals for color of dot in plot
# xy$color <- cut(x = xy$padj,
#     breaks = c(1.00,0.10,0.05,0.01,0.001),
#     labels = c("red", "orange","gold","black"),
#     right = FALSE,
#     include.lowest = TRUE
#     )
# grades$Letter = cut(grades$Final, 
#                     breaks = c(0,60,70,80,90, 100), 
#                     labels = c("F", "D", "E", "B", "A"),
#                     right = FALSE, 
#                     include.lowest = TRUE)

# make xy plot
#plot(x = xy$fed.mean, y = xy$str.mean, col=xy$color,pch=19,log="xy")

## make plot in ggplot2
library(ggplot2)

xy <- xy[!is.na(xy$padj),]
xy.sig <- xy[xy$padj <= 0.1 & !is.na(xy$padj),]


g <- ggplot(data = xy) + geom_point(size=2,alpha=0.5,aes(x=fed.mean,y=str.mean, col=padj) ) + 
    scale_x_log10() + scale_y_log10() + theme_dark() + 
    labs( x = "Fed", y = "Starved", title ="Gene Expression in AHK Cells (Mean Raw Counts)") + 
    scale_color_gradientn(
        name="adjusted p-value",trans="log10",
        colours = c("red","red","orange", "gold","purple", "darkgrey","darkgrey"),
        breaks = c(1.00,0.10,0.05,0.01,0.001))


g + geom_point(data = xy.sig, size=2,aes(x=fed.mean,y=str.mean, col=padj) ) +
    scale_color_gradientn(
        name="adjusted p-value",trans="log10",
        colours = c("red","red","orange", "gold","mediumorchid1", "darkgrey","darkgrey"),
        breaks = c(1.00,0.10,0.05,0.01,0.001)) +
    geom_text(data = xy.sig,  size = 3,
              # check_overlap = TRUE, 
              # hjust = +1, nudge_x = 0.05,
              position = position_nudge(x = 0.2, y = 0.05),
              # position = position_jitter( width=0.2, height=0.2),
              aes(x=fed.mean, y=str.mean, label=names,size=1,
                  # angle = 45
                  ))
ggsave(filename = "plots/07.xyplot.pdf",width = 8.5,height = 8.5 )

# ---------------------------------------------
# make 3d graph when rayshader package for ggplot comes out ----
#install.packages("rayshader")


# ---------------------------------------------
# chromesome package  ----

# ---------------------------------------------
## save objects for next step ----
#save(dds,rld,res,file = "outs/06.dds.rld.res.RData")
load("outs/06.dds.rld.res.RData")

# ---------------------------------------------
# Notes ----
## See 7 in: 
### https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
### start on 6 in next script
### for single cell:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis