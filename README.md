# RNAseq
README last updated 08/10/2024
## 📝 Description
- Functions and analysis scripts for bulk RNAseq analysis
## crick-pipeline:
- Contains functions and analysis scripts for  post crick RNAseq pipeline
### extract_count_mats.R
- Script for extracting count matrices from DESEq2 object output from pipeline
### extractPCA.R
- Function to extract PCA and meta data dataframe from RNAseq pipeline DESeq2 object
### GSEA_analysis.R
- Analysis script for GSEA
### runGSEA.R
- Function for running GSEA and creating plots
- Required inputs:
    - *DE_Table* = the differential expression results table for the comparison of interest
    - *pathways* = gmt file read in by gmtfile
    - *min* = minimum pathway size
    - *max* = maximum pathway size
    - *saveRes* = whether results should be saved as CSV (True or False)
    - *plotTop* = whether top results should be plotted (True or False)
    - *compName* = name of differential comparison name DE table is from
    - *saveDir* = where outputs should be saved 
## custom-scripts: 
- Contains any scripts that have been made for projects where some sort of custom analysis was required
## qmul-hpc-scripts: 
- Is an archive of old scripts previously run on QMUL HPC.
## nf-core-test: 
- Contains scripts for testing the nf-core pipeline