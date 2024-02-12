#!/bin/bash

# Script for running fastqc on fastq files
# Author: James Boot, Date: 22/08/2023
# This script is designed to run on QMUL HPC

#####
#$ -M j.boot@qmul.ac.uk                # Change to your email address 
#$ -m bes
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=2G  
#$ -j y
#$ -l h_rt=240:00:00
#$ -N FastQC                           # Change the name of the job accordingly
#####

# SECTION 1: Load module
module load fastqc

# SECTION 2: Loop through files
for INPUT_FILE in *.fastq.gz; do
  fastqc ${INPUT_FILE}
done

# SECTION 3: Run multiqc on above output
# Load gcenv containing multiqc
. /data/WHRI-GenomeCentre/gcenv/bin/activate

# Run multiqc
multiqc .
