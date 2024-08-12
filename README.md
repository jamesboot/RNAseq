# RNAseq

## 📝 Description
- Scripts for bulk RNAseq analysis
- Scripts in this repository have not been tested

## 🔩 Getting started
### Dependencies
- conda environment with Squire installed is required
- Use STAR installation available on Apocrita
- Please see in all scripts for parameters that need to be added by the user - after the Apocrita queue options

### Scripts

1. run_QCandTrim.sh
    - Trimming using trimgalore with FASTQC before and after trimming.
    - ❗ NOTE: multiqc will need to be run manually after this - script needs to be added for this.
    - See fastqc.sh script for another example 

2. STAR
    - Folder contains all STAR associated scripts
    - run_STAR_index.sh
        - If index is not already prepared
    - run_STAR_unique.sh
        - Alignment keeping only unique alignments.
    - run_STAR_unique.sh
        - Alignment keeping one random alignment.

3. Squire
    - Folder contains all Squire associated scripts
    - Documentation: https://github.com/wyang17/SQuIRE
    - squire_Map.sh
        - Map reads using Squire
    - suire_Count.sh
        - Count reads using Squire