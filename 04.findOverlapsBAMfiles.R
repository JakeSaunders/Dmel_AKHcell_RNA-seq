# ---------------------------------------------
# set working dir, load packages, and load objects ----
setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

## Load objects
load("outs/03.TxDb.cdsBy.RData")

## load Genomic Alignments

#install.packages("BiocManager")
#library(BiocManager)
#install(pkgs = "GenomicAlignments")
library(GenomicAlignments)

# ---------------------------------------------
# make list of paths to bam files (each bam file is one sample) ----

bam.files <- list.files(
  path = "aligned/",
  pattern = ".bam$",
  recursive = T,
  full.names = T
)

# ---------------------------------------------
# find "overlaps between reference transcriptome (cdsBy.TxDb.Hsap.GRCh37.75) and bam files ----

## pick which TxDb object you want to use
ebg <- Dmel.TxDb.transcriptsBy          # This could really change results check back here

# use summarizeOverlaps to find "overlaps between reference transcriptome (cdsBy.TxDb.Hsap.GRCh37.75) and bam files    
## ignore.strand=TRUE - counting function to ignore the strand, i.e., to allow minus strand reads to count in plus strand genes, and vice versa.
### https://genomicsclass.github.io/book/pages/read_counting.html

se <- summarizeOverlaps(features=ebg, 
                        reads=bam.files,
                        mode="Union",   # Determines how reads that overlap are counted "Union" conservative
                        singleEnd=TRUE,
                        ignore.strand=TRUE
)

# takes about 20 minutes

# ------------------------------------------------
#Check and change column names to be more readable and add sampleTable ----
colnames(se)
colnames(se) <- sub(pattern = ".Aligned.sortedByCoord.out.bam",replacement = "",x = colnames(se))
colnames(se)


## use this code to add sample metadata to se object from either csv file or type dataframe
SummarizedExperiment::colData(se) <- DataFrame(
  short.name = colnames(se),
  file.name = paste0(colnames(se),".Aligned.sortedByCoord.out.bam"),
  sample.num = as.factor(sprintf("%02d", 1:10)),
  dependent.var = as.factor(c(rep("fed",5),rep("starved",5)))
)


# ------------------------------------------------
# check features of count table

dim(se)
assayNames(se)
head(assay(se))
tail(assay(se))
str(assay(se))
class(assay(se))
dim(assay(se))
colSums(assay(se))

reads.mills <- round( colSums(assay(se)) / 1e6, 1 )
boxplot(reads.mills[1:5], notch=T,ylim=c(0,18), main="Fed", ylab="# of Reads")
boxplot(reads.mills[6:10],notch=T,ylim=c(0,18), main="Starved", ylab="# of Reads")

rowRanges(se)
str(metadata(rowRanges(se)))

# make and save raw count table

count.table <- assay(se)
write.csv(ct, "outs/04.countTable.raw.csv",row.names = T,col.names = T)

# ------------------------------------------------
# Save and Load objects produced by this code ----
## Save 

save(
  ebg,
  se,
  ct,
  file = "outs/04.ebg.se.RData"
)
save(se,file = "outs/04.ebg.se.RData")

## Load

load("outs/04.ebg.se.RData")
# ---------------------------------------------
# Notes ----

## See 2.5-2.8 Defining gene models in" 
### https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#preparing-count-matrices
### After this read 3.0 but start with 4.0
