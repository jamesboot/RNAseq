#!/bin/bash
#SBATCH --job-name=nfcore-test
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --partition=ncpu
#SBATCH --wrap "sh /nemo/stp/babs/working/bootj/nfcore-test/nfcore-test.sh"