#!/bin/bash

# Script for aligning raw RNAseq data (fastq files) to a reference genome using STAR
# Author: James Boot, 22/08/2023
# This script is designed to run on QMUL HPC

# Set queue options
#$ -m e 					                # Email when completed
#$ -M j.boot@qmul.ac.uk 	        # Send email to this address
#$ -cwd						                # Set the working directory for the job to current directory
#$ -pe smp 4				              # Request 4 cores
#$ -l h_vmem=8G				            # Request 8GB RAM
#$ -j y						                # Join stdout and stderr - standard output and standard error stream
#$ -l h_rt=24:00:00			          # Request 10 hrs runtime
#$ -t 1-58					              # Threading option - this should match the number of samples 

# Specify the parameters file containing sample names 
SAMPLES=$(sed -n "${SGE_TASK_ID}p" Samples.txt)
FASTQDIR=/data/WHRI-GenomeCentre/shares/Projects/Bioinformatic/Lin_Yung-Yao/BulkRNAseq/Data/SRA_Download

# Load STAR module
module load star/2.7.0f

# Run STAR with --quantMode GeneCounts so that gene level counts are automatically output - then htseq does not need to be run
STAR --outSAMstrandField intronMotif \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--genomeDir /data/WHRI-GenomeCentre/Genome/Human/hg19.GRCh38.104/star \
--readFilesIn ${FASTQDIR}/${SAMPLES}_1.fastq ${FASTQDIR}/${SAMPLES}_2.fastq \
--outFileNamePrefix ${SAMPLES}_ \
--runThreadN ${NSLOTS}	
