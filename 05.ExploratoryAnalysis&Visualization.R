# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

## Load objects
load("outs/04.ebg.se.RData")
rm(ebg,ct)

## load packages 
#install.packages("BiocManager")
library(BiocManager)
#install(pkgs = "DESeq2")
library("DESeq2")

# ---------------------------------------------
# 3.1 Starting from SummarizedExperiment ----

SummarizedExperiment::colData(se)
dds <- DESeq2::DESeqDataSet(se, design = ~dependent.var )
dds
#save(dds,file = "outs/05.dds.raw.RData")

# ---------------------------------------------
# 4.1 Pre-filtering the dataset ----

## total number of rows/genes
nrow(dds)
## Number of FALSE = number of rows that have no data
table(rowSums(counts(dds)) > 1)
## remove those rows
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# ---------------------------------------------
# 4.2 Deciding between vsn, log(x+1) and rlog transformation ----
## NOTE: There are all kinds of typos and mistakes in this section 
## I rewroted the code around to get it to work and produce the same effect as the examples
#install(pkgs = "vsn")
library(vsn)

## create objects transformed both ways
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)

## make a dataframes for each transformation so let me do plot to see differences in transformations
sample.col <- c(2,8) # I selected sample 2 and 8 because they were different conditions but had similar total reads
log2.df <- as_data_frame(log2(counts(dds, normalized=TRUE)[, sample.col]+1)) %>%
  mutate(transformation = "log2(x + 1)")
colnames(log2.df) <- c("x","y","transformation")
vst.df <- as_data_frame(assay(vsd)[, sample.col]) %>% mutate(transformation = "vst")
colnames(vst.df) <- c("x","y","transformation")
rlog.df <- as_data_frame(assay(rld)[, sample.col]) %>% mutate(transformation = "rlog")
colnames(rlog.df) <- c("x","y","transformation")
df <- bind_rows( log2.df, vst.df, rlog.df)

## make plot showing effect of different transformations
ggplot(data = df, aes(x = df$x, y = df$y)) + geom_point(alpha = 0.10, color="darkgreen") + 
  coord_fixed() + facet_grid( . ~ transformation) + 
  labs(title ="Three Variance Stabilizing Transformation Strategies", 
         x = "Counts of 02.Fed", 
         y = "Counts of 08.Str")

## command to save plot
# ggsave("plots/05.TransformationStrategies.pdf", width = 11, height = 8.5)

## From pipeline write up:
###Scatterplot of transformed counts from two samples. Shown are scatterplots using the log2 transform of normalized counts (left), using the VST (middle), and using the rlog (right). While the rlog is on roughly the same scale as the log2 counts, the VST has a upward shift for the smaller values. It is the differences between samples (deviation from y=x in these scatterplots) which will contribute to the distance calculations and the PCA plot.
###We can see how genes with low counts (bottom left-hand corner) seem to be excessively variable on the ordinary logarithmic scale, while the VST and rlog compress differences for the low count genes for which the data provide little information about differential expression.

# rlog looks to be the best of the stragies because it shrinks the distance between 
# some of the lower expressed genes better, and those can drive differences that are not real
# also rlog tends to work well on small datasets (n < 30), 
# potentially outperforming the VST when there is a wide range of sequencing depth across samples 
# (an order of magnitude difference). Therefore use VST for medium-to-large datasets (n > 30). 

# I'm going to be use rld going forward

## clear works space of unneeded objects for next step
rm(sample.col,log2.df,rlog.df,vst.df,df,vsd)

# ---------------------------------------------
# 4.3 Sample distances - assess overall similarity between samples - Heatmaps ----
## We use the R function dist to calculate the Euclidean distance between samples.
### To ensure we have a roughly equal contribution from all genes, we use it on the RLD data. 
### We need to transpose the matrix of values using t, because the dist function expects the 
### different samples to be rows of its argument, and different dimensions (here, genes) to be columns.
sampleDists <- dist(t(assay(rld)))
sampleDists 

## We visualize the distances in a heatmap ---
#install.packages("pheatmap")
library("pheatmap")
library("RColorBrewer")

## To plot the sample distance matrix with the rows/columns arranged by the distances in our distance matrix, 
## we manually provide sampleDists to the clustering_distance argument of the pheatmap function.
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- dds$short.name
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#pdf(file = "plots/05.heatmap.simularity.pdf",width = 8,height = 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distances (Measure of Simularity)")
#dev.off()

## Poisson Distance---
###  Poisson Distance (Witten 2011) is a measure of dissimilarity between counts, 
### also takes the inherent variance structure of counts into consideration when calculating the distances between samples. 

install.packages("PoiClaClu")
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$dex, dds$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- dds$short.name
colnames(samplePoisDistMatrix) <- NULL

#pdf(file = "plots/05.heatmap.dissimilarity.pdf",width = 8,height = 8)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         main = "Sample-to-Sample Poisson Distance (Measure of Dissimilarity)")
#dev.off()

# ---------------------------------------------
# 4.4 principal components analysis (PCA) plot ----
## take a quick look at PCA
plotPCA(rld, intgroup = c("dependent.var"))

## save plot as object modify pca in ggplot
pcaData <- plotPCA(rld, intgroup = c("dependent.var"), returnData = TRUE)
pcaData

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = dependent.var, shape = dependent.var)) +
  geom_point(size =6) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(label = "PCA Plot of AHK Cells from Fed and Starved Flies") +
  coord_fixed() + 
  scale_color_manual(values = c("darkorchid","DarkGreen")) +
  scale_shape_manual(values = c(19,13)) +
  geom_text(aes(label=name),hjust=0.5, vjust=0.5,col="black")
## command to save plot
# ggsave("plots/05.PCA.pdf", width = 11, height = 8.5)

# ---------------------------------------------
## save objects for next step ----
#save(dds,se,rld,vsd,file = "outs/05.dds.se.rld.vsd.RData")
load("outs/05.dds.se.rld.vsd.RData")

# ---------------------------------------------
# Notes ----
## See 3 to 4 in: 
### https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#preparing-count-matrices
### After this read 5 Differential expression analysis