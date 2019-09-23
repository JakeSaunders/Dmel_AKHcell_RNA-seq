#!/bin/bash
#SBATCH --job-name="01.QCtrimmingQC"
#SBATCH --partition=medium
#SBATCH --account="generalGrp"
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=7-00:00:00
#SBATCH --mem=30GB
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --output="01.QCtrimQC.o"
#SBATCH --error="01.QCtrimQC.e"
#SBATCH --mail-user=saundecj@wfu.edu
#################################### Now do your stuff ###################################

cd /deac/generalGrp/johnsonGrp/saundecj/projects/akhSingleCell/

## Ran from command line before runing this script to QC raw reads
# fastqc -o QC/pretrim raws/*.*.fastq.gz &

## trim reads using trimmoatics
mkdir trim

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/01.fed.fastq.gz \
	trim/01.fed.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/02.fed.fastq.gz \
	trim/02.fed.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/03.fed.fastq.gz \
	trim/03.fed.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/04.fed.fastq.gz \
	trim/04.fed.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/05.fed.fastq.gz \
	trim/05.fed.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36
	
java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/06.str.fastq.gz \
	trim/06.str.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/07.str.fastq.gz \
	trim/07.str.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/08.str.fastq.gz \
	trim/08.str.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/09.str.fastq.gz \
	trim/09.str.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

java -jar  /deac/generalGrp/johnsonGrp/saundecj/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 \
	raws/10.str.fastq.gz \
	trim/10.str.trim.fastq.gz \
	ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:36

####### run slurm file
## Run QC on trimmed reads to check effect

fastqc -o QC/postrim trim/*.*.trim.fastq.gz &

exit 0

##Added following entries to base trimmomatic TruSeq3-SE file
# >TruSeqUniversalAdapter
# AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
# >TruSeqAdapterIndex23
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex10
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex2
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex3
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex4
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex8
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex9
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex6
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
# >TruSeq AdapterIndex11
# GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG

