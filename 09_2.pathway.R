# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")


## load packages 
#install.packages("BiocManager")
library(BiocManager)
library("DESeq2")
library("GenomicFeatures", lib.loc="~/R/win-library/3.6")
library(tidyverse)



### GAGE (what I used last time ) ----
###### https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

#install("gage")
library("gage")

# install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

## Load objects
load("outs/06.dds.rld.res.RData")

# ---------------------------------------------
# Get kegg lists for metabolism & signaling pathways AND Disease pathways----

#Generate up-to-date KEGG pathway gene sets for any specified KEGG species
kg.dme <- kegg.gsets(species = "dme")   

############################################## NONE OF THIS MAKES SENSE IN FLIES B.C THERE IS ONLY ONE DIEASE PATHWAY
#### this code selects only the metabolism pathways, change sigmet.ix to below to change lists
# kg.sets - KEGG gene sets, a named list. Each element is a character vector of member gene IDs for a single KEGG pathway. The number of elements of this list is the total number of KEGG pathways defined for the specified species.
# sigmet.idx - integer indice, which elements in kg.sets are signaling or metabolism pathways.
# sig.idx - integer indice, which elements in kg.sets are signaling pathways.
# met.idx - integer indice, which elements in kg.sets are metabolism pathways.
# dise.idx - integer indice, which elements in kg.sets are disease pathways.

# KEGG metabolism & signaling pathways
#sigMet.kegg.gs <- kg.dme$kg.sets[kg.dme$sigmet.idx]
#
# KEGG metabolism & Disease pathways
#dise.kegg.gs <- kg.dme$kg.sets[kg.dme$dise.idx]
############################################ THIS WOULD BE USEFUL FOR OTHER SPECIES

length(kg.dme$kg.sets) # all pathways
length(kg.dme$sig.idx) # signaling
length(kg.dme$met.idx) # metabolic
length(kg.dme$dise.idx) # diease 

# just use all pathways not a subset
kegg.gs <- kg.dme$kg.sets

# ---------------------------------------------
# Get info to Convert row names to CG numbers ----

# install bioconductor package that has fly gene names
#install("org.Dm.eg.db")
library("org.Dm.eg.db")

# Using org.Dm.eg.db to get a dataframe of flybase names and gene symbols
columns(org.Dm.eg.db)   # see info stored in package/database

key.type <- "SYMBOL"
symbol.keys <- keys(org.Dm.eg.db,keytype = key.type) # list all of a paticular list from db
symbol.df <- AnnotationDbi::select(    
    org.Dm.eg.db,                       # db to look up information/symbols
    keys=symbol.keys,                   # select records from db
    columns=c("SYMBOL","FLYBASECG"),    # columns of dataframe to return
    keytype="SYMBOL"                    # this is the argument use in
)

symbol.df <- symbol.df[ order(x = symbol.df$SYMBOL), ]

# make normalized count table dataframe with symbol as column
cnts.norm <- as.data.frame(assay(rld))
colnames(cnts.norm) <- rld$short.name
cnts.norm$SYMBOL <- rownames(cnts.norm)

# there are some rownames that aren't in the symbol list 221 of 11233
table(rownames(cnts.norm) %in% symbol.df$SYMBOL)

# add flybasecg column
cnts.norm <- left_join(cnts.norm,symbol.df,by="SYMBOL")

# if there is a missing flybasecg name for gene add old symbol name to list
cnts.norm[ is.na(cnts.norm$FLYBASECG) , ]$FLYBASECG <- cnts.norm[ is.na(cnts.norm$FLYBASECG) , ]$SYMBOL

# add Dmel_ prefixt to cg's
cnts.norm$FLYBASECG<- paste0("Dmel_", cnts.norm$FLYBASECG)

#make FLYBASECG column new row names
rownames(cnts.norm) <- cnts.norm$FLYBASECG

#drop symbol and flybase columns
cnts.norm <- cnts.norm[,1:10]


# ---------------------------------------------
# run the GAGE (Generally Applicable Gene-set Enrichment) function ----


## same.dir will often change depending on source of pathways, FALSE for KEGG, BioCarta, TRUE for GO term groups
### run once as same.dir = T

cnts.kegg.p.SDT <- gage(
    exprs = cnts.norm,       # count table
    gsets = kegg.gs,         # list of kegg pathways
    ref = 1:5,               # column numbers for control samples - fed here
    samp = 6:10,             # column numbers for experimental samples - starved or str here
    same.dir = T,            # in KEGG pathways, genes frequently are not coregulated (e.g. up and down regulated) but mulitiple comparisons are better if TRUE
    compare = "as.group" ,   # compares ref to samp columns
    FDR.adj = TRUE,          # this data might be cleaner than most RNA seq data sets, consider reruning with false later
    use.fold = T,            # use t-test stat rather than fold change
    rank.test = F            # T use non-parametric test, since lots of small values
)

### run once as same.dir = F
cnts.kegg.p.SDF <- gage(
    exprs = cnts.norm,       # count table
    gsets = kegg.gs,         # list of kegg pathways
    ref = 1:5,               # column numbers for control samples - fed here
    samp = 6:10,             # column numbers for experimental samples - starved or str here
    same.dir = F,            # in KEGG pathways, genes frequently are not coregulated (e.g. up and down regulated) but mulitiple comparisons are better if TRUE
    compare = "as.group" ,   # compares ref to samp columns
    FDR.adj = TRUE,          # this data might be cleaner than most RNA seq data sets, consider reruning with false later
    use.fold = T,            # use t-test stat rather than fold change
    rank.test = F            # T use non-parametric test, since lots of small values
)
    ## what the output columns of gage mean
        # p.geomean	- geometric mean of the individual p-values from multiple single array based gene set tests
        # stat.mean	- mean of the individual statistics from multiple single array based gene set tests. Normally, its absoluate value measures the magnitude of gene-set level changes, and its sign indicates direction of the changes. When saaTest=gs.KSTest, stat.mean is always positive.
        # p.val	- gloal p-value or summary of the individual p-values from multiple single array based gene set tests. This is the default p-value being used.
        # q.val	- FDR q-value adjustment of the global p-value using the Benjamini & Hochberg procedure implemented in multtest package. This is the default q-value being used.
        # set.size	- the effective gene set size, i.e. the number of genes included in the gene set test
        # other columns	-  columns of the individual p-values or statistics, each measures the gene set perturbation in a single experiment (vs its control or all controls, depends on the "compare argument value)

# make and save one summary table to send to Erik
sdt.g <- as.data.frame(cnts.kegg.p.SDT$greater)
sdt.g$keggPathway <- rownames(sdt.g)
sdt.g$direction.change <- "greater"
sdt.l <- as.data.frame(cnts.kegg.p.SDT$less)
sdt.l$keggPathway <- rownames(sdt.l)
sdt.l$direction.change <- "lesser"
sdt.s <- as.data.frame(cnts.kegg.p.SDF$greater)
sdt.s$keggPathway <- rownames(sdt.s)
sdt.s$direction.change <- "absolute.diff"
gage.summary <- rbind.data.frame(sdt.g,sdt.l,sdt.s)
rm(sdt.s,sdt.g,sdt.l)

gage.summary <- gage.summary %>% select(direction.change,q.val,p.val,everything())

# write.csv(gage.summary,file = "outs/09_2.kegg.pathway.summary.csv",row.names = T)

# ---------------------------------------------
# sigGeneSet extracts sig data and plots heat (heatmap does not work) ----

SDF.kegg.sig <- sigGeneSet(setp = cnts.kegg.p.SDF,cutoff = 0.1,qpval = "p.val",
           heatmap = F) # heat map pdf can't be openned
SDT.kegg.sig <- sigGeneSet(setp = cnts.kegg.p.SDT,cutoff = 0.1,qpval = "p.val",
           heatmap = F) # heat map pdf can't be openned

# ---------------------------------------------
# Pathview ----

#install("pathview")
library(pathview) 

# select sig diff pathways (p <0.1) for ploting
table( gage.summary$p.val <= 0.1 )
sel <- gage.summary[ gage.summary$p.val <= 0.1 , ]
sel <- sel[ !is.na(sel$keggPathway), ]

# select unique kegg pathway names to plot, since duplicates increase running time
path.ids <- unique(sel$keggPathway)

# plot pathways
?pathview()

# plot most sig pathway
pathview(
    gene.data = cnts.norm,     # count table, (vignette did a weird normalize i skiped: cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx]))
    species = "dme",            # KEGG species code
    pathway.id = path.ids[1],   # pathway to be ploted, this code picks pathway with lowest pval
    gene.idtype = "KEGG",
    kegg.native = FALSE,        # plot as pdf if FALSE
    pdf.size = c(8.5,11),
    kegg.dir = "pathway"
        )
# loop to make pdfs for all pathways using graphviz 
sapply(
    X = path.ids, 
    FUN = function(pid) {
        pathview(
            gene.data = cnts.norm,   # count table, (vignette did a weird normalize i skiped: cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx]))
            species = "dme",         # KEGG species code
            pathway.id = pid,        # pathway to be ploted, this code picks pathway with lowest pval
            gene.idtype = "KEGG",
            kegg.native = FALSE,     # plot as pdf if FALSE
            pdf.size = c(8.5,11),
            kegg.dir = "pathway"
        )
    }
)

# loop to make pngs for all pathways using kegg graphs 
sapply(
    X = path.ids, 
    FUN = function(pid) {
        pathview(
            gene.data = cnts.norm,   # count table, (vignette did a weird normalize i skiped: cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx]))
            species = "dme",         # KEGG species code
            pathway.id = pid,        # pathway to be ploted, this code picks pathway with lowest pval
            gene.idtype = "KEGG",
            kegg.native = TRUE,     # plot as pdf if FALSE
            pdf.size = c(8.5,11),
            kegg.dir = "pathway"
        )
    }
)


# move pdf files into pathway dir
pdfs <- list.files(path = ".",pattern = "*.pdf")
file.rename(from = pdfs,to = paste0("pathway/",pdfs))
xmls <- list.files(path = ".",pattern = "*.xml")
file.rename(from = xmls,to = paste0("pathway/",xmls))
pngs <- list.files(path = ".",pattern = "*.png")
file.rename(from = pngs,to = paste0("pathway/",pngs))
rm(pdfs,xmls)



# ---------------------------------------------
# save items for next step ----

save(cnts.norm,gage.summary,SDF.kegg.sig,SDT.kegg.sig,file = "outs/09_2.pathway.RData")


# ---------------------------------------------
# NOTES:
# KEGG PATHWAY ANALYSIS AND PATHVIEW but uses microarray:
# Generally Applicable Gene-set/Pathway Analysis
# https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/gage.pdf

# RNA-SEQ TO KEGG PATHWAY AND GO TERM 
# RNA-Seq Data Pathway and Gene-set Analysis Work
# https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf

# PATHVIEW EXPLAINED
# Pathview: pathway based data integration and visualization
# https://www.bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf

# ---------------------------------------------
# OTHER NOTES
# https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/dataPrep.pdf
### New 2019: https://www.bioconductor.org/packages/devel/bioc//vignettes/ToPASeq/inst/doc/ToPASeq.html
# https://bioconductor.org/packages/release/bioc/vignettes/DEsubs/inst/doc/DEsubs.pdf
# https://www.bioconductor.org/packages/devel/bioc//vignettes/ToPASeq/inst/doc/ToPASeq.html
### https://bioconductor.org/packages/devel/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html#pathway-analysis
### https://bioconductor.org/packages/devel/bioc/vignettes/ideal/inst/doc/ideal-usersguide.html
### https://bioconductor.org/packages/3.9/bioc/html/goseq.html
