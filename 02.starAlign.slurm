#!/bin/bash
#SBATCH --job-name="02.starAlign"
#SBATCH --partition=medium
#SBATCH --account="generalGrp"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=7-00:00:00
#SBATCH --mem=60GB
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output="02.starAlign.o"
#SBATCH --error="02.starAlign.e"
#SBATCH --mail-user=saundecj@wfu.edu
#################################### Now do your stuff ###################################


## run from command line before running slurm file
#wget ftp://ftp.ensembl.org/pub/release-96/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz -O refs/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
#gunzip refs/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa.gz
#STAR --runMode genomeGenerate --runThreadN 12 --genomeDir refs/ --genomeFastaFiles refs/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa

## Align fastq's to genome 

cd /deac/generalGrp/johnsonGrp/saundecj/projects/akhSingleCell/

module load star/2.5

for fq in `ls trim/*`; do

    sample=`basename $fq .trim.fastq.gz`
	mkdir aligned/$sample/
	
    STAR \
		--runThreadN 16 \
		--genomeDir refs/ \
		--readFilesIn $fq \
		--readFilesCommand zcat \
		--outFileNamePrefix aligned/$sample/$sample. \
		--outSAMtype BAM SortedByCoordinate 
done
 
module unload star/2.5

exit 0
