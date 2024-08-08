#!/bin/bash

# Note that this scripts requires a conda environment with squire installed

#$ -l h_rt=240:00:0
#$ -l h_vmem=16G
#$ -N sMap
#$ -pe smp 4
#$ -cwd
#$ -j y

# User: Edit here only
# ANALYSISDIR should be the directory where the script will be saved
# TRIMDIR should be the directory where trimmed files were previously - Not sure if this is what should be used!
# SQUIREDIR should be the directory where the SQUIRE output will be saved
# SQUIREDIR does not have to exist - it will be made when the script runs
# FETCHFOLDER is the folder location of outputs from SQuIRE Fetch
# READLENGTH is the read length of RNAseq reads
# GTF is GTF file of genome transcripts 
# R1PATTERN and R2PATTERN should be a regular expression of the naming convention of read 1 and read 2
ANALYSISDIR=
TRIMDIR=${ANALYSISDIR}/trimgalore_out
SQUIREDIR=${ANALYSISDIR}/squire_map_out
FETCHFOLDER=
READLENGTH=
GTF=
R1PATTERN=".*_1.fq.gz"
R2PATTERN=".*_2.fq.gz"

# Log
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "=========================================================="

# Load modules and activate conda env
# Recommend specifying the version of miniconda and always using this
module load miniconda
conda activate squire

# Specify output files to be generated 
R1LISTFILE=${ANALYSISDIR}/trim_R1_files.txt
R2LISTFILE=${ANALYSISDIR}/trim_R2_files.txt
SAMPLELISTFILE=${ANALYSISDIR}/trim_samples.txt

# Locate files 
# Find read 1's
find ${TRIMDIR} -regex ${R1PATTERN} | sort >> ${R1LISTFILE}
# Find read 2's
find ${TRIMDIR} -regex ${R2PATTERN} | sort >> ${R2LISTFILE}
# Get sample names 
#sed 's:.*/::' trim_R1_files.txt | sed 's:_.*::' >> ${SAMPLELISTFILE}
sed 's:.*/::' trim_R1_files.txt | sed 's:_1_val_1.fq.gz::' >> ${SAMPLELISTFILE}

# Now load saved file names and locations
R1FILES=${ANALYSISDIR}/trim_R1_files.txt
R2FILES=${ANALYSISDIR}/trim_R2_files.txt
SAMPLE_NAMES=${ANALYSISDIR}/trim_samples.txt

# Creating a SQUIREDIR directory (if it doesn't exist), 
# where the results of running STAR will be saved.
mkdir -p ${SQUIREDIR}

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
    echo "Running SQUIRE for: ${SAMPLE}"
    echo "=========================================================="

    # Mapping
    squire Map -1 ${R1} \
    -2 ${R2} \
    -r ${READLENGTH} \
    -b hg38 \
    -p ${NSLOTS} \
    -o ${SQUIREDIR} \
    -n ${SAMPLE} \
    -f ${FETCHFOLDER} \
    -g ${GTF}

    echo "=========================================================="
    echo "Finished analysis for sample: ${SAMPLE}"
    echo "Finished on : $(date)"
    echo "=========================================================="

done