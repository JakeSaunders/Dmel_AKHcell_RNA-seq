# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

## Load objects
load("outs/04.ebg.se.RData")
load("outs/05.dds.se.rld.vsd.RData")
rm(ebg,ct,se,vsd)

## load packages 
#install.packages("BiocManager")
library(BiocManager)
#install(pkgs = "DESeq2")
library("DESeq2")

# ---------------------------------------------
#   5.1. Running the differential expression pipeline ----

dds <- DESeq(dds)
    ## The DESeq() function does the following on the raw 
    # using pre-existing size factors
    # estimating dispersions
    # gene-wise dispersion estimates
    # mean-dispersion relationship
    # final dispersion estimates
    # fitting model and testing

# ---------------------------------------------

dds <- DESeq(dds)
## The DESeq() function does the following on the raw 
# using pre-existing size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

# ---------------------------------------------
#   5.2 Building the results table -----

res <- results(dds)
    ## If we hadn't already defined the experimental design earlier we would have to do it now
    ## res <- results(dds, contrast=c("dependent.var"))

mcols(res, use.names = TRUE)
    ##  results retruns a special dataframe object with metadata on the meaning of the column names
     ##  mcol() returns that metadata
            # baseMean       intermediate            mean of normalized counts for all samples
            # log2FoldChange      results log2 fold change (MLE): dependent.var starved vs fed
            # lfcSE               results         standard error: dependent.var starved vs fed
            # stat                results         Wald statistic: dependent.var starved vs fed
            # pvalue              results      Wald test p-value: dependent.var starved vs fed
            # padj                results                                 BH adjusted p-values

## view summary of results
summary(res)

## how many sig genes?
table(res$padj < 0.1)

## which genes are sig diff?
padj.list <- rows <- res$padj 
padj.list[ is.na(padj.list) ] <- 1
sig.genes <- as.data.frame(res[padj.list < 0.1,])
write.csv(sig.genes,file = "outs/06.padj.less.0.01.csv")

# ---------------------------------------------
## save objects for next step ----
#save(dds,rld,res,file = "outs/06.dds.rld.res.RData")
load("outs/06.dds.rld.res.RData")

# ---------------------------------------------
# Notes ----
## See 5 in: 
### https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#preparing-count-matrices
### start on 6 in next script

