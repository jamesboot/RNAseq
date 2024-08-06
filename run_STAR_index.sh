#!/bin/bash

# Set queue options
#$ -m bes                                # Email when completed
#$ -M j.boot@qmul.ac.uk                   # Send email to this address
#$ -cwd                                   # Set the working directory for the job to current directory
#$ -pe smp 8                              # Request 8 cores
#$ -l h_vmem=32G                          # Request 32GB RAM
#$ -j y                                   # Join stdout and stderr: standard output and standard error stream
#$ -l h_rt=06:00:00                       # Request 6 hrs runtime
#$ -l highmem                             # Join high memory queue

# Load module
module load star/2.7.0f

# Run STAR in genomeGenerate mode
# When generating a new index the --genomeDir --genomeFastaFiles and --sjdbGTFfile arguements will need updating
STAR --runThreadN ${NSLOTS} \
--runMode genomeGenerate \
--genomeDir /path/to/output/folder/star \
--genomeFastaFiles /path/to/genome/fasta/file.fa \
--sjdbGTFfile /path/to/genome/gtf/file.gtf \
--sjdbOverhang 100 \
--limitGenomeGenerateRAM=168632718037