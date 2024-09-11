#!/bin/bash
#SBATCH --job-name=nfcore-test
#SBATCH --nodes=1
#SBATCH -c=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --partition=ncpu
#SBATCH --wrap "sh /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/smrnaseq.sh"