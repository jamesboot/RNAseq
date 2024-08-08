#!/bin/bash

# Note that this scripts requires a conda environment with squire installed

#$ -l h_rt=240:00:0
#$ -l h_vmem=15G
#$ -N sCount
#$ -pe smp 1
#$ -cwd
#$ -j y

# User: Edit here only
# ANALYSISDIR should be the directory where the script will be saved
# MAPDIR should be the directory of outputs from SQuIRE Map
# SQUIREDIR should be the directory where the SQUIRE output will be saved
# SQUIREDIR does not have to exist - it will be made when the script runs
# FETCHFOLDER is the folder location of outputs from SQuIRE Fetch
# CLEANFOLDER is the folder location of outputs from SQuIRE Clean
# READLENGTH is the read length of RNAseq reads
# STRANDEDNESS '0' if unstranded eg Standard Illumina, 1 if first- strand eg Illumina Truseq, dUTP, NSR, NNSR, 2 if second-strand, eg Ligation, Standard SOLiD (optional,default=0)
ANALYSISDIR=
MAPDIR=${ANALYSISDIR}/squire_map_out
SQUIREDIR=${ANALYSISDIR}/squire_out
FETCHFOLDER=
CLEANFOLDER=
READLENGTH=
STANDEDNESS=1

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

# Load trim sample names used for squire map
SAMPLE_NAMES=${ANALYSISDIR}/trim_samples.txt

# For loop to run through all samples and files
ITERATIONS=$(wc -l < ${SAMPLE_NAMES})

for i in $(seq ${ITERATIONS}); do

	SAMPLE=$(sed -n "${i}p" ${SAMPLE_NAMES})

    echo "=========================================================="
    echo "Starting on : $(date)"
    echo "Starting analysis for sample: ${SAMPLE}"
    echo "=========================================================="

    echo "=========================================================="
    echo "Running SQUIRE for: ${SAMPLE}"
    echo "=========================================================="

    # Squire Count
    squire Count -m ${MAPDIR} \
    -r ${READLENGTH} \
    -c ${CLEANFOLDER} \
    -f ${FETCHFOLDER} \
    -n ${SAMPLE} \
    -b hg38 \
    -s ${STRANDEDNESS}

    echo "=========================================================="
    echo "Finished analysis for sample: ${SAMPLE}"
    echo "Finished on : $(date)"
    echo "=========================================================="

done