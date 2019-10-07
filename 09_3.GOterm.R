# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")


## load packages 
#install.packages("BiocManager")
library(BiocManager)
library("DESeq2")
library("GenomicFeatures", lib.loc="~/R/win-library/3.6")
library(tidyverse)
library("gage")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("org.Dm.eg.db")
library("GO.db")

# library(pathview) 

# ---------------------------------------------
## load for format count table 
cnts.norm <- read.csv("outs/09.summary.table.csv")
cnts.norm <- cnts.norm[, c(8,9:18)]


# ---------------------------------------------
# Generate up-to-date GO gene sets for 

go.sets.fly <- go.gsets(species = "Fly")

# ---------------------------------------------
## convert gene ids to GO terms
key.type <- "SYMBOL"
symbol.keys <- keys(org.Dm.eg.db,keytype = key.type) # list all of a paticular list from db

symbol.df <- AnnotationDbi::select(    
    org.Dm.eg.db,                       # db to look up information/symbols
    keys=symbol.keys,                   # select records from db
    columns=c("SYMBOL","GO","GOALL"),    # columns of dataframe to return
    keytype="SYMBOL"                    # this is the argument use in
)
save(symbol.df,file = "outs/09_3.symbol.df.Rdata")



cnts.norm <- left_join(cnts.norm, symbol.df[c(1:2,5)],by="SYMBOL")


##### START HERE
# rm(symbol.df)
# cnts.norm[ is.na(cnts.norm$GO) , ]$FLYBASECG <- cnts.norm[ is.na(cnts.norm$GO) , ]$GOALL


# ---------------------------------------------
# GAGE FUNCTION


test <- gage(
    cnts.norm, 
    gsets = go.sets.hs[go.subs.hs$MF],
    ref = 1:5, 
    samp = 6:10, 
    compare ="as.group"
)



#Molecular Function analysis is quicker, hence run as demo
# cnts.mf.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$MF],
#                   + ref = ref.idx, samp = samp.idx, compare ="unpaired")
#Biological Process analysis takes a few minutes if you try it
#cnts.bp.p <- gage(cnts.norm, gsets = go.sets.hs[go.subs.hs$BP],
# ref = ref.idx, samp = samp.idx, compare ="unpaired")




# ---------------------------------------------
# GO analysis and other gene set analyses ----



# ---------------------------------------------
# ---------------------------------------------
# ---------------------------------------------
# ---------------------------------------------
# NOTES:

# RNA-SEQ TO KEGG PATHWAY AND GO TERM 
# RNA-Seq Data Pathway and Gene-set Analysis Work
# https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf