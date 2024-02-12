#!/bin/bash 
 
# Script to download fastq files using SRA toolkit from NCBI

##### Set queue options
#$ -M j.boot@qmul.ac.uk            # Change to your email address 
#$ -m bes
#$ -cwd
#$ -pe smp 2
#$ -l h_vmem=8G
#$ -l h_rt=240:00:00
#$ -N YungYao_SRA-FASTQ-DUMP
#$ -j y
#####

module load sratools/2.10.8

for SRR in $(cat GSE159_SRR_Acc_List.txt)
do
	fasterq-dump --split-files $SRR
done

for SRR in $(cat GSE189_SRR_Acc_List.txt)
do
	fasterq-dump --split-files $SRR
done