# Function to extract PCA and meta data dataframe from RNAseq pipeline DESeq2 object
extractPCA <- function(PCA_rda, sample_set) {
  load(PCA_rda)
  PC_df <- as.data.frame(ddsList[[sample_set]]@colData)
  return(PC_df)
}