# Script for aligning trimmed reads to reference genome using STAR
# This script will only report unique alignments

#!/bin/bash

#$ -M j.boot@qmul.ac.uk
#$ -m bes
#$ -l h_rt=240:00:0
#$ -pe smp 2
#$ -l h_vmem=8G
#$ -N STAR_unique
#$ -cwd
#$ -j y

# User: Edit here only
# ANALYSISDIR should be the directory where the script will be saved
# TRIMDIR should be the directory where trimmed files were previously
# STARDIR should be the directory where the STAR output will be saved
# STARDIR does not have to exist - it will be made when the script runs
# GENOME should be the full directory to the STAR index 
# R1PATTERN and R2PATTERN should be a regular expression of the naming convention of read 1 and read 2
ANALYSISDIR=
TRIMDIR=${ANALYSISDIR}/trimgalore_out
STARDIR=${ANALYSISDIR}/STAR_unique_out
GENOME=
FASTA=
GTF=
R1PATTERN=".*_1.fq.gz"
R2PATTERN=".*_2.fq.gz"


# User: Do not edit beyond this point

# Log
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : ${JOB_ID}"
echo "=========================================================="

# Load the programs that will be required
module load star/2.7.0f
module load samtools/1.19.2

# Specify output files to be generated 
R1LISTFILE=${ANALYSISDIR}/trim_R1_files_uniq.txt
R2LISTFILE=${ANALYSISDIR}/trim_R2_files_uniq.txt
SAMPLELISTFILE=${ANALYSISDIR}/trim_samples_uniq.txt

# Locate files 
# Find read 1's
find ${TRIMDIR} -regex ${R1PATTERN} | sort >> ${R1LISTFILE}
# Find read 2's
find ${TRIMDIR} -regex ${R2PATTERN} | sort >> ${R2LISTFILE}
# Get sample names 
#sed 's:.*/::' trim_R1_files_uniq.txt | sed 's:_.*::' >> ${SAMPLELISTFILE}
sed 's:.*/::' trim_R1_files_uniq.txt | sed 's:_1_val_1.fq.gz::' >> ${SAMPLELISTFILE}

# Now load saved file names and locations
R1FILES=${ANALYSISDIR}/trim_R1_files_uniq.txt
R2FILES=${ANALYSISDIR}/trim_R2_files_uniq.txt
SAMPLE_NAMES=${ANALYSISDIR}/trim_samples_uniq.txt

# Creating a STARDIR directory (if it doesn't exist), 
# where the results of running STAR will be saved.
mkdir -p ${STARDIR}

# For loop to run through all samples and files
ITERATIONS=$(wc -l < ${SAMPLE_NAMES})

for i in $(seq ${ITERATIONS}); do

    R1=$(sed -n "${i}p" ${R1FILES})	
	R2=$(sed -n "${i}p" ${R2FILES})	
	SAMPLE=$(sed -n "${i}p" ${SAMPLE_NAMES})

    echo "=========================================================="
    echo "Starting on : $(date)"
    echo "Starting analysis for sample: ${SAMPLE}"
    echo "Read 1: ${R1}"
    echo "Read 2: ${R2}"
    echo "=========================================================="

    echo "=========================================================="
    echo "Running STAR for: ${SAMPLE}"
    echo "=========================================================="

    # STAR alignment
    STAR --runThreadN ${NSLOTS} \
    --genomeDir ${GENOME} \
    --readFilesIn ${R1} ${R2} \
    --outSAMtype BAM SortedByCoordinate \
    --runMode alignReads \
    --outFileNamePrefix ${STARDIR}/${SAMPLE}_ \
    --outFilterMultimapNmax 1  \
    --outFilterMismatchNmax 3 \
    --alignEndsType EndToEnd \
    --readFilesCommand zcat \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF}

    samtools index ${STARDIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam

    echo "=========================================================="
    echo "Finished analysis for sample: ${SAMPLE}"
    echo "Finished on : $(date)"
    echo "=========================================================="

done