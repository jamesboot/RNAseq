# Script for custom plots for af787

# Load packages
library(DESeq2)
library(ggplot2)

# Source functions
source('/nemo/stp/babs/working/bootj/github/RNAseq/crick-pipeline/extractPCA.R')

# Create outdir
outdir <- 'custom_outs'
if (!dir.exists(outdir)){
  dir.create(outdir)
}

# Extract PCA
pcadf <- extractPCA(PCA_rda = 'data/pca_GRCh38_analyse.rda',
                    sample_set = 'HD_ET_Only')

# Plot
ggplot(pcadf, aes(x = .PCA.PC1, y = .PCA.PC2)) +
  geom_point(size = 6,
    aes(
    colour = Pool,
    shape = Culture
  )) +
  ggtitle('PCA visualising all variables') +
  theme_bw(base_size = 14) +
    facet_wrap('CellType', nrow = 1)
ggsave(
  paste0(outdir, '/pca_celtype_culture_pool.pdf'),
  plot = last_plot(),
  height = 5,
  width = 10,
  units = 'in'
)
