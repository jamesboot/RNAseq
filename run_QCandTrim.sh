# Script for running fastqc, trimming and factqc on trimmed files

#!/bin/bash

#$ -M j.boot@qmul.ac.uk
#$ -m bes
#$ -l h_rt=240:00:0
#$ -pe smp 2
#$ -l h_vmem=8G
#$ -N QCandTrim
#$ -cwd
#$ -j y

# User: Edit here only
# ANALYSISDIR should be the directory where the script will be saved
# FASTQDIR should be the directory where the FASTQ files are saved
# QCDIR should be the directory where the QC of the raw FASTQ files should be saved
# TRIMDIR should be the directory where trimmed files should be saved
# QCDIR and TRIMDIR do not have to exist - they will be made when the script runs
# R1PATTERN and R2PATTERN should be a regular expression of the naming convention of read 1 and read 2
ANALYSISDIR=
FASTQDIR=
QCDIR=${ANALYSISDIR}/1.fastQC_out
TRIMDIR=${ANALYSISDIR}/2.trim_out
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
module load trimgalore/0.6.5
module load fastqc/0.11.9

# Specify output files to be generated 
R1LISTFILE=${ANALYSISDIR}/R1_file_list.txt
R2LISTFILE=${ANALYSISDIR}/R2_file_list.txt
SAMPLELISTFILE=${ANALYSISDIR}/sample_name_list.txt

# Locate files 
# Find read 1's
find ${FASTQDIR} -regex ${R1PATTERN} | sort >> ${R1LISTFILE}
# Find read 2's
find ${FASTQDIR} -regex ${R2PATTERN} | sort >> ${R2LISTFILE}
# Get sample names 
#sed 's:.*/::' R1_file_list.txt | sed 's:_.*::' >> ${SAMPLELISTFILE}
sed 's:.*/::' R1_file_list.txt | sed 's:_1.fq.gz::'>> ${SAMPLELISTFILE}

# Now load saved file names and locations
R1FILES=${ANALYSISDIR}/R1_file_list.txt
R2FILES=${ANALYSISDIR}/R2_file_list.txt
SAMPLE_NAMES=${ANALYSISDIR}/sample_name_list.txt

# Creating a QCDIR and a TRIMDIR directory (if they don't exist), where the results of running
# fastqc and trimgalore, respectively, will be saved.
mkdir -p ${QCDIR}
mkdir -p ${TRIMDIR}

# For loop to run through all samples and files
ITERATIONS=$(wc -l < ${SAMPLE_NAMES})

for i in $(seq ${ITERATIONS}); do

	R1=$(sed -n "${i}p" ${R1FILES})	
	R2=$(sed -n "${i}p" ${R2FILES})	
	SAMPLE=$(sed -n "${i}p" ${SAMPLE_NAMES})

    echo "*******************"
    echo "Starting on : $(date)"
    echo "Starting analysis for sample: ${SAMPLE}"
    echo "Read 1: ${R1}"
    echo "Read 2: ${R2}"
    echo "*******************"

    echo "*******************"
    echo "Running FASTQC for Read 1"
    echo "*******************"
    fastqc ${R1} -t 1 -o ${QCDIR}

    echo "*******************"
    echo "Running FASTQC for Read 2"
    echo "*******************"
    fastqc ${R1} -t 1 -o ${QCDIR}

    # Running trimgalore to remove adapters. The --fastqc option will cause fastqc to run after trimming.
    # This should allow you to check that adapters are no longer present.
    echo "*******************"
    echo "Running TrimGalore"
    echo "*******************"
    trim_galore --fastqc --paired -o ${TRIMDIR} ${R1} ${R2}

    echo "*******************"
    echo "Finished analysis for sample: ${SAMPLE}"
    echo "Finished on : $(date)"
    echo "*******************"

done