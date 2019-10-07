# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")


## load packages 
#install.packages("BiocManager")
library(BiocManager)
library("DESeq2")
library("GenomicFeatures", lib.loc="~/R/win-library/3.6")
library(tidyverse)

## Load objects
load("outs/06.dds.rld.res.RData")

### GAGE (what I used last time ) ----
###### https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

#install("gage")
library("gage")

# install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

# ---------------------------------------------
# Convert row names to CG numbers ----

# install bioconductor package that has fly gene names
#install("org.Dm.eg.db")
library("org.Dm.eg.db")

# Using org.Dm.eg.db to get a dataframe of different gene names
columns(org.Dm.eg.db)   # see info stored in package/database

key.type <- "SYMBOL"
symbol.keys <- keys(org.Dm.eg.db,keytype = key.type) # list all of a paticular list from db

symbol.df <- AnnotationDbi::select(    
    org.Dm.eg.db,                       # db to look up information/symbols
    keys=symbol.keys,                   # select records from db
    columns=c("SYMBOL","GENENAME", "FLYBASECG", "UNIPROT"),    # columns of dataframe to return
    keytype="SYMBOL"                    # this is the argument use in
)

symbol.df <- symbol.df[ order(x = symbol.df$SYMBOL), ]
# write.csv(symbol.df,file = "outs/09.geneSymbolKey.csv",col.names = T,row.names = F)

# get other data frames formatted for merger 
res.ct <- as.data.frame(res)
res.ct$SYMBOL <- row.names(res.ct) # add rownames as column

rld.ct <- as.data.frame(assay(rld))
rld.ct$SYMBOL <- row.names(rld.ct)
colnames(rld.ct) <- c( paste0("norm.",rld$short.name), "SYMBOL")


dds.ct <- as.data.frame(assay(dds))
dds.ct$SYMBOL <- row.names(dds.ct)
colnames(dds.ct) <- c(paste0("raw.",dds$short.name), "SYMBOL")

# make one big summary file with rld counts, dds counts

master.df <- left_join(res.ct, rld.ct, by="SYMBOL")
master.df <- left_join(master.df,dds.ct,by="SYMBOL")
master.df <- left_join(master.df, symbol.df, by="SYMBOL")

# combine UNIPROT, FLYBASECG, GENENAME entries if more than one
# note the method for combining duplicates
master.df <- master.df %>% 
    dplyr::group_by( 
            baseMean , log2FoldChange , lfcSE , stat , pvalue , padj , SYMBOL , 
            norm.01.fed , norm.02.fed ,  norm.03.fed , norm.04.fed , norm.05.fed ,
            norm.06.str, norm.07.str , norm.08.str , norm.09.str , norm.10.str , 
            raw.01.fed , raw.02.fed , raw.03.fed , raw.04.fed , raw.05.fed , 
            raw.06.str , raw.07.str , raw.08.str , raw.09.str , raw.10.str 
    ) %>% 
    summarise(
        UNIPROT = paste0( unique(UNIPROT),collapse = ";"), 
        FLYBASECG = paste0( unique(FLYBASECG),collapse = ";"),
        GENENAME =  paste0( unique(GENENAME),collapse = ";")
    )

#order based on adjusted pvalue
master.df <- master.df[ base::order(master.df$padj),]
# add row names
rownames(master.df) <- master.df$SYMBOL
#save dataframe
#write.csv(x = master.df,file = "outs/09.summary.table.csv",row.names = T)


# Notes ----
## See 7 in: 
### http://www.circos.ca/
### start on 6 in next script
### for single cell:https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis