## LOAD REQUIRED MODULES
ml purge
ml Nextflow/22.10.3
ml Singularity/3.6.4
ml CAMP_proxy

export NXF_HOME=/nemo/stp/babs/working/mitterr/.nextflow/

## Set NXF_WORK as the absolute NEMO path to the folder containing the pipeline scripts (i.e. not CAMP on a softlink)
cd /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq
export NXF_WORK=$(pwd -P)/work

## Uncomment the code below to set the NXF_WORK in the scratch space.
#export NXF_WORK=`echo $NXF_WORK | sed 's:^/nemo/stp/babs/working/:/flask/scratch/babs/:'`
#if [ ! -d "$NXF_WORK" ]; then
#  ln -s $NXF_WORK .
#fi

# Singularity is the system we use with nf-core to make software available in a consistent manner.
# You can leave this line as is - we are using the BABS library of singularity images
export NXF_SINGULARITY_LIBRARYDIR=/flask/apps/containers/all-singularity-images/

# You should create a folder in your working area to store any additional images that are not found in the BABS library folder above
export NXF_SINGULARITY_CACHEDIR=/flask/apps/misc/stp/babs/nf-core/singularity/

## MOTE: For the contaminant files (--rrna, --trna, --cdna, --ncrna, pirna, other_contamination), each FASTA file is first compared to the available miRNA sequences and overlaps are removed.  Thus, no need to pre-filter.

## RUN PIPELINE
nextflow run nf-core/smrnaseq \
  --input /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/samplesheet.csv \
  --protocol custom \
  --three_prime_adapter AACTGTAGGCACCATCAAT \
  --genome GRCh38 \
  --fasta /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/GRCh38.fa \
  --mirtrace_species hsa \
  --mirna_gtf /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/mirbase/hsa.gff3 \
  --hairpin /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/mirbase/hairpin.fa \
  --mature /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/mirbase/mature.fa \
  --email richard.mitter@crick.ac.uk \
  --outdir results \
  -profile crick \
  -c /nemo/stp/babs/working/mitterr/projects/patanir/hamish.crerar/hc682_deseq2/smrnaseq/nf_config.txt \
  -r 2.2.1 \
  -resume
