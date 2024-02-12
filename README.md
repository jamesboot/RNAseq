# RNAseq
README last updated 12/02/2024
## 📝 Description
- *Example* scripts for bulk RNAseq analysis
- Scripts in this project were run on published datasets with GEO accession numbers: GSE189053 and GSE159273
- Results from this analysis are being published and are currently under review
## 🔩 Getting started
### Dependencies
- SRA Toolkit
- FastQC & MultiQC
- STAR
- Post-alignment analysis in R:
```
library(edgeR)
library(tidyverse)
library(gplots)
library(biomaRt)
library(data.table)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(fgsea)
```
### Process
1. Download samples using SRA Toolkit using `sra-fastq-dump.sh`
   - This script requires lists of sample accession numbers for files to download. Accession lists can be downloaded from NCBI SRA Run Selector.
2. Perform initial QC using `fastqc.sh`
   - This script should be run in the directory containing downloaded `.fastq` files.
3. Align samples to reference genome using STAR aligner in `star-align.sh`
   - This script is run as an array - creating a job for each sample to parallelise the process
   - This script requires a list of sample names in a `.txt` file as input an the folder containing `.fastq` files
4. Perform exploratory analysis, differential expression analysis and pathway analysis, along with visualisations in R using `Analysis.R`

