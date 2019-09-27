# ---------------------------------------------
# set working dir ----

setwd("C:/Users/saundecj/Google Drive/fly/AKHcells.rnaSeq")

# ---------------------------------------------
# Install bioconductor and Genomic Features Package ----
## commented out install package commands

#install.packages("BiocManager")
#library(BiocManager)
#install(pkgs = "GenomicFeatures")

## load packages

library(GenomicFeatures)

# ---------------------------------------------
# Download GTF or GFF file which corrseponds to genome used in alignment step (see 02.starAlign.slurm) ----
## fastq files were aligned to release 96 Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
## find gff3 file on ensembl.org

file.url <- "ftp://ftp.ensembl.org/pub/release-96/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gff3.gz"

### Download gtf file from ensembl

#curl::curl_download(url = file.url, destfile = "refs/Drosophila_melanogaster.BDGP6.22.96.gff3.gz")

### unzip file 

R.utils::gunzip("refs/Drosophila_melanogaster.BDGP6.22.96.gff3.gz")

# ---------------------------------------------
# Make a TxDb object from annotations available as a GFF3 or GTF file ----
## The TxDb class is a container for storing transcript annotations.
## helpful: https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf

### Use this dataframe to pass metadata to TxDb object constructed with makeTxDb functions
meta.data <- data.frame(
  name=c("Genome"),
  value=c("Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa")
)

Dmel.TxDb <- makeTxDbFromGFF(
  file = "refs/Drosophila_melanogaster.BDGP6.22.96.gff3",
  organism = "Drosophila melanogaster",
  metadata = meta.data
)

### Check output to make sure metadata in TxDb object is correct
Dmel.TxDb

# ---------------------------------------------
# Pull coding sequences from TxDb object ----
## cdsBy to look at protein coding transcripts
## by="gene" Determines the grouping of objects 
## look at use.names=TRUE to see if I need to use it here

## cdsBy for all protein coding transcripts

?cdsBy

Dmel.TxDb.cdsBy  <- cdsBy(Dmel.TxDb, by="gene")
Dmel.TxDb.exonsBy  <- exonsBy(Dmel.TxDb, by="gene")
Dmel.TxDb.transcriptsBy <- transcriptsBy(Dmel.TxDb, by ="gene")

# ---------------------------------------------
# How many genes? ----
## Reasonable amount for 
### Human ~23K
### Fly ~14K

length(Dmel.TxDb.cdsBy)

# ---------------------------------------------
# Save and Load objects produced by this code ----
## Save 

save(
  Dmel.TxDb.cdsBy,
  Dmel.TxDb.exonsBy,
  Dmel.TxDb.transcriptsBy,
  file = "outs/03.TxDb.cdsBy.RData"
)

## Load

load("outs/03.TxDb.cdsBy.RData")

# ---------------------------------------------
# Notes ----
## see make make "count tables"

## See 2.5 Defining gene models in" 
### https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#preparing-count-matrices

### https://bioconductor.org/packages/release/bioc/vignettes/GenomicAlignments/inst/doc/summarizeOverlaps.pdf

### https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#use-with-downstream-bioconductor-dge-packages
