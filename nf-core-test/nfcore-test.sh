# LOAD REQUIRED MODULES
ml purge
ml Nextflow/24.04.2
ml Singularity/3.11.3
ml Java/21.0.2

# Singularity is the system we use with nf-core to make software available in a consistent manner.
# You can leave this line as is - we are using the BABS library of singularity images
export NXF_SINGULARITY_LIBRARYDIR=/flask/apps/containers/all-singularity-images/

# You should create a folder in your working area to store any additional images that are not found in the BABS library folder above
export NXF_SINGULARITY_CACHEDIR=/flask/apps/misc/stp/babs/nf-core/singularity/

# TEST PIPELINE
nextflow run nf-core/rnaseq \
  --outdir outs \
  -profile test,crick \
  -r 3.15.0 \
  -resume
