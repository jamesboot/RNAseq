# Script for custom plots and analysis for Bauer project b830

# Check BiocManager install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load packages
library(DESeq2)
library(ggplot2)
library(readxl)
library(dplyr)
library(RColorBrewer)

# Create directory for outputs to be saved to if doesn't exist
outdir <- 'R_custom_outs'
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# Load data
load('babs/differential/data/pca_GRCm38_analyse.rda')

# Load in meta data
meta <- read_xlsx('meta/meta.xlsx')

# Define colours for later
# Variant colours
variantCols <- brewer.pal(5, "Set1")
variants <- factor(
  c('PBS', 'BA.1', 'BA.5', 'BA.2.86', 'XBB'),
  levels = c('PBS', 'BA.1', 'BA.5', 'BA.2.86', 'XBB')
)
names(variantCols) <- levels(variants)

# Dose colours
doseCols <- brewer.pal(3, "Set2")
doses <- factor(c('Low', 'Med', 'High'), levels = c('Low', 'Med', 'High'))
names(doseCols) <- levels(doses)

# Age colours
ageCols <- brewer.pal(2, "Dark2")
ages <- factor(c('Young', 'Old'), levels = c('Young', 'Old'))
names(ageCols) <- levels(ages)

####
# Plot PCAs ----
####

# EXPERIMENT 1 - Variant: plot Ct values onto PCA
# Extract PCA data
pcaDF <- data.frame(
  ID = ddsList$experiment1_variant$ID,
  PC1 = ddsList$experiment1_variant$.PCA$PC1,
  PC2 = ddsList$experiment1_variant$.PCA$PC2,
  Variant = ddsList$experiment1_variant$Variant
)

# Check levels
pcaDF$Variant <- factor(pcaDF$Variant, levels = c('PBS', 'BA.1', 'BA.5'))

# Add Ct values to pcaDF
pcaDF$Ct_TaqPath_N <- NA
for (x in c(1:nrow(pcaDF))) {
  pcaDF$Ct_TaqPath_N[x] <- as.numeric(meta$Ct_TaqPath_N[meta$`LIMS RNAseq` == as.character(pcaDF$ID[x])])
}

# Plot PCA coloured by Ct
mid <- mean(pcaDF$Ct_TaqPath_N)
ggplot(pcaDF, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Ct_TaqPath_N),
    colour = 'black',
    size = 6,
    shape = 21
  ) +
  scale_fill_gradient2(
    midpoint = mid,
    low = 'red',
    high = 'blue',
    mid = 'white'
  ) +
  ggtitle('Experiment 1 - Variant, Med Dose, DPI3, PCA coloured by Ct_TaqPath_N') +
  theme_bw(base_size = 18)
ggsave(
  paste0(outdir, '/exp1_variant_pca_Ct.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot PCA coloured by Variant
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = Variant)) +
  geom_point(size = 6) +
  ggtitle('Experiment 1 - Variant, Medium Dose, DPI3, PCA coloured by Variant') +
  theme_bw(base_size = 18) +
  scale_colour_manual(name = "Variant", values = variantCols)
ggsave(
  paste0(outdir, '/exp1_variant_pca_Variant.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# EXPERIMENT 1 - Dose: plot Ct values onto PCA
# Extract PCA data
pcaDF <- data.frame(
  ID = ddsList$experiment1_dose$ID,
  PC1 = ddsList$experiment1_dose$.PCA$PC1,
  PC2 = ddsList$experiment1_dose$.PCA$PC2,
  Variant = ddsList$experiment1_dose$Variant,
  Dose = ddsList$experiment1_dose$Dose
)

# Check levels
pcaDF$Variant <- factor(pcaDF$Variant, levels = c('PBS', 'BA.1', 'BA.5'))
pcaDF$Dose <- factor(pcaDF$Dose, levels = c('Low', 'Med'))

# Add Ct values to pcaDF
pcaDF$Ct_TaqPath_N <- NA
for (x in c(1:nrow(pcaDF))) {
  pcaDF$Ct_TaqPath_N[x] <- as.numeric(meta$Ct_TaqPath_N[meta$`LIMS RNAseq` == as.character(pcaDF$ID[x])])
}

# Plot PCA coloured by Ct
mid <- mean(pcaDF$Ct_TaqPath_N)
ggplot(pcaDF, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Ct_TaqPath_N),
    colour = 'black',
    size = 6,
    shape = 21
  ) +
  scale_fill_gradient2(
    midpoint = mid,
    low = 'red',
    high = 'blue',
    mid = 'white'
  ) +
  ggtitle('Experiment 1 - Dose, DPI3, PCA coloured by Ct_TaqPath_N') +
  theme_bw(base_size = 18)
ggsave(
  paste0(outdir, '/exp1_dose_pca_Ct.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot PCA coloured by Variant
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = Variant)) +
  geom_point(size = 6) +
  ggtitle('Experiment 1 - Dose, Medium Dose, DPI3, PCA coloured by Variant') +
  theme_bw(base_size = 18) +
  scale_colour_manual(name = "Variant", values = variantCols)
ggsave(
  paste0(outdir, '/exp1_dose_pca_Variant.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot PCA coloured by Dose
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = Dose)) +
  geom_point(size = 6) +
  ggtitle('Experiment 1 - Dose, Medium Dose, DPI3, PCA coloured by Variant') +
  theme_bw(base_size = 18) +
  scale_colour_manual(name = "Dose", values = doseCols)
ggsave(
  paste0(outdir, '/exp1_dose_pca_Dose.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# EXPERIMENT 2 - Variant: plot Ct values onto PCA
# Extract PCA data
pcaDF <- data.frame(
  ID = ddsList$experiment2_variant$ID,
  PC1 = ddsList$experiment2_variant$.PCA$PC1,
  PC2 = ddsList$experiment2_variant$.PCA$PC2,
  Variant = ddsList$experiment2_variant$Variant
)

# Check levels
pcaDF$Variant <- factor(pcaDF$Variant, levels = c('PBS', 'BA.2.86', 'XBB'))

# Add Ct values to pcaDF
pcaDF$Ct_TaqPath_N <- NA
for (x in c(1:nrow(pcaDF))) {
  pcaDF$Ct_TaqPath_N[x] <- as.numeric(meta$Ct_TaqPath_N[meta$`LIMS RNAseq` == as.character(pcaDF$ID[x])])
}

# Plot PCA coloured by Ct
mid <- mean(pcaDF$Ct_TaqPath_N)
ggplot(pcaDF, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Ct_TaqPath_N),
    colour = 'black',
    size = 6,
    shape = 21
  ) +
  scale_fill_gradient2(
    midpoint = mid,
    low = 'red',
    high = 'blue',
    mid = 'white'
  ) +
  ggtitle('Experiment 2 - Variant, Med Dose, DPI3, PCA coloured by Ct_TaqPath_N') +
  theme_bw(base_size = 18)
ggsave(
  paste0(outdir, '/exp2_variant_pca_Ct.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot PCA coloured by Variant
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = Variant)) +
  geom_point(size = 6) +
  ggtitle('Experiment 2 - Variant, Medium Dose, DPI3, PCA coloured by Variant') +
  theme_bw(base_size = 18) +
  scale_colour_manual(name = "Variant", values = variantCols)
ggsave(
  paste0(outdir, '/exp2_variant_pca_Variant.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# EXPERIMENT 3 - Variant: plot Ct values onto PCA
# Extract PCA data
pcaDF <- data.frame(
  ID = ddsList$experiment3$ID,
  PC1 = ddsList$experiment3$.PCA$PC1,
  PC2 = ddsList$experiment3$.PCA$PC2,
  Dose = ddsList$experiment3$Dose,
  Age = ddsList$experiment3$Age
)

# Check levels
pcaDF$Dose <- factor(pcaDF$Dose, levels = c('Low', 'Med', 'High'))
pcaDF$Age <- factor(pcaDF$Age, levels = c('Young', 'Old'))

# Add Ct values to pcaDF
pcaDF$Ct_TaqPath_N <- NA
for (x in c(1:nrow(pcaDF))) {
  pcaDF$Ct_TaqPath_N[x] <- as.numeric(meta$Ct_TaqPath_N[meta$`LIMS RNAseq` == as.character(pcaDF$ID[x])])
}

# Plot PCA coloured by Ct
mid <- mean(pcaDF$Ct_TaqPath_N)
ggplot(pcaDF, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(fill = Ct_TaqPath_N),
    colour = 'black',
    size = 6,
    shape = 21
  ) +
  scale_fill_gradient2(
    midpoint = mid,
    low = 'red',
    high = 'blue',
    mid = 'white'
  ) +
  ggtitle('Experiment 3, PCA coloured by Ct_TaqPath_N') +
  theme_bw(base_size = 18)
ggsave(
  paste0(outdir, '/exp3_pca_Ct.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot PCA coloured by Dose
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = Dose)) +
  geom_point(size = 6) +
  ggtitle('Experiment 3, PCA coloured by Dose') +
  theme_bw(base_size = 18) +
  scale_colour_manual(name = "Dose", values = doseCols)
ggsave(
  paste0(outdir, '/exp3_pca_Dose.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot PCA coloured by Age
ggplot(pcaDF, aes(x = PC1, y = PC2, colour = Age)) +
  geom_point(size = 6) +
  ggtitle('Experiment 3, PCA coloured by Age') +
  theme_bw(base_size = 18) +
  scale_colour_manual(name = "Age", values = ageCols)
ggsave(
  paste0(outdir, '/exp3_pca_Age.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

####
# Correlate Viral counts and Ct ----
####

# Load the raw counts matrix
tpm <- read.delim('babs/nfcore/results/GRCm38/star_rsem/rsem.merged.gene_tpm.tsv',
                  sep = '\t')
row.names(tpm) <- tpm$gene_id
tpm$gene_id <- NULL
tpm$transcript_id.s. <- NULL

# Extract viral counts
row <- grep('NC', rownames(tpm))
viraltpm <- data.frame(ID = colnames(tpm), ViralTPM = as.numeric(tpm[row, ]))

# Load experiment table
experimentTable <- read.csv('babs/docs/experiment_table.csv')

# Create Ct Vals dataframe
CtVals <- data.frame(ID = meta$`LIMS RNAseq`,
                     Ct_TaqPath_N = as.numeric(meta$Ct_TaqPath_N))

# Any undertemined Ct values are converted to NA - replace with Ct of 36
CtVals$Ct_TaqPath_N[is.na(CtVals$Ct_TaqPath_N)] <- 36

# Merge CtVals and viraltpm
correlation <- merge(CtVals, viraltpm, by = 'ID')

# Merge experimentTable and correlation
correlation <- merge(correlation, experimentTable, by = 'ID')

# Log transform viral TPM 
correlation$logViralTPM <- log(correlation$ViralTPM + 1)

# Filter to Experiment 1 - Variant samples only
exp1varCorr <- correlation %>%
  filter(Experiment == 1, Dose == 'Med', DPI == 'DPI3')

# Check levels
exp1varCorr$Variant <- factor(exp1varCorr$Variant, levels = c('PBS', 'BA.1', 'BA.5'))
exp1varCorr$Dose <- factor(exp1varCorr$Dose, levels = c('Low', 'Med'))

# Plot correlation for Experiment 1 - Variant
ggplot(exp1varCorr, aes(x = Ct_TaqPath_N, y = logViralTPM)) +
  geom_point(aes(colour = Variant), size = 6) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 1 - Variant, correlation of logViralTPM and Ct_TaqPath_N') +
  scale_colour_manual(name = "Variant", values = variantCols)
ggsave(
  paste0(outdir, '/exp1_variant_corr_Variant.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot Ct value as paired box plot
ggplot(exp1varCorr, aes(x = Variant, y = Ct_TaqPath_N)) +
  geom_boxplot(aes(color = Variant), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Variant),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Variant", values = variantCols) +
  scale_fill_manual(name = "Variant", values = variantCols) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 1 - Variant, distribution of Ct_TaqPath_N values')
ggsave(
  paste0(outdir, '/exp1_variant_Cts.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot TPM value as paired box plot
ggplot(exp1varCorr, aes(x = Variant, y = logViralTPM)) +
  geom_boxplot(aes(color = Variant), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Variant),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Variant", values = variantCols) +
  scale_fill_manual(name = "Variant", values = variantCols) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 1 - Variant, distribution of logViralTPM values')
ggsave(
  paste0(outdir, '/exp1_variant_TPM.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Filter to Experiment 1 - Dose samples only
exp1doseCorr <- correlation %>%
  filter(Experiment == 1, DPI == 'DPI3')

# Check levels
exp1doseCorr$Variant <- factor(exp1doseCorr$Variant, levels = c('PBS', 'BA.1', 'BA.5'))
exp1doseCorr$Dose <- factor(exp1doseCorr$Dose, levels = c('Low', 'Med'))

# Plot correlation for Experiment 1 - Dose
ggplot(exp1doseCorr, aes(x = Ct_TaqPath_N, y = logViralTPM)) +
  geom_point(aes(colour = Variant), size = 6) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 1 - Dose, correlation of logViralTPM and Ct_TaqPath_N') +
  scale_colour_manual(name = "Variant", values = variantCols) +
  facet_wrap( ~ Dose)
ggsave(
  paste0(outdir, '/exp1_dose_corr_Variant.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot Ct value as paired box plot
ggplot(exp1doseCorr, aes(x = Variant, y = Ct_TaqPath_N)) +
  geom_boxplot(aes(color = Variant), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Variant),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Variant", values = variantCols) +
  scale_fill_manual(name = "Variant", values = variantCols) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 1 - Dose, distribution of Ct_TaqPath_N values') +
  facet_wrap( ~ Dose)
ggsave(
  paste0(outdir, '/exp1_dose_Cts.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot TPM value as paired box plot
ggplot(exp1doseCorr, aes(x = Variant, y = logViralTPM)) +
  geom_boxplot(aes(color = Variant), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Variant),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Variant", values = variantCols) +
  scale_fill_manual(name = "Variant", values = variantCols) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 1 - Dose, distribution of logViralTPM values') +
  facet_wrap( ~ Dose)
ggsave(
  paste0(outdir, '/exp1_dose_TPM.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Filter to Experiment 2 - Variant samples only
exp2varCorr <- correlation %>%
  filter(Experiment == 2, DPI == 'DPI3')

# Check levels
exp2varCorr$Variant <- factor(exp2varCorr$Variant, levels = c('PBS', 'BA.2.86', 'XBB'))

# Plot correlation for Experiment 2 - Variant
ggplot(exp2varCorr, aes(x = Ct_TaqPath_N, y = logViralTPM)) +
  geom_point(aes(colour = Variant), size = 6) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 2 - Variant, correlation of logViralTPM and Ct_TaqPath_N') +
  scale_colour_manual(name = "Variant", values = variantCols)
ggsave(
  paste0(outdir, '/exp2_variant_corr_Variant.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot Ct value as paired box plot
ggplot(exp2varCorr, aes(x = Variant, y = Ct_TaqPath_N)) +
  geom_boxplot(aes(color = Variant), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Variant),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Variant", values = variantCols) +
  scale_fill_manual(name = "Variant", values = variantCols) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 2 - Variant, distribution of Ct_TaqPath_N values')
ggsave(
  paste0(outdir, '/exp2_variant_Cts.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot TPM value as paired box plot
ggplot(exp2varCorr, aes(x = Variant, y = logViralTPM)) +
  geom_boxplot(aes(color = Variant), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Variant),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Variant", values = variantCols) +
  scale_fill_manual(name = "Variant", values = variantCols) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 2 - Variant, distribution of logViralTPM values')
ggsave(
  paste0(outdir, '/exp2_variant_TPM.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Filter to Experiment 3 samples only
exp3Corr <- correlation %>%
  filter(Experiment == 3)

# Check levels
exp3Corr$Dose <- factor(exp3Corr$Dose, levels = c('Low', 'Med', 'High'))
exp3Corr$Age <- factor(exp3Corr$Age, levels = c('Young', 'Old'))

# Plot correlation for Experiment 3
ggplot(exp3Corr, aes(x = Ct_TaqPath_N, y = logViralTPM)) +
  geom_point(aes(colour = Age), size = 6) +
  theme_bw(base_size = 18) +
  ggtitle('Experiment 3, correlation of logViralTPM and Ct_TaqPath_N') +
  scale_colour_manual(name = "Age", values = ageCols) +
  facet_wrap( ~ Dose)
ggsave(
  paste0(outdir, '/exp3_corr_Dose_Age.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot Ct value as paired box plot
ggplot(exp3Corr, aes(x = Dose, y = Ct_TaqPath_N)) +
  geom_boxplot(aes(color = Dose), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Dose),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Dose", values = doseCols) +
  scale_fill_manual(name = "Dose", values = doseCols) +
  theme_bw(base_size = 18) +
  facet_wrap( ~ Age) +
  ggtitle('Experiment 3, distribution of Ct_TaqPath_N values')
ggsave(
  paste0(outdir, '/exp3_Cts.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)

# Plot TPM value as paired box plot
ggplot(exp3Corr, aes(x = Dose, y = logViralTPM)) +
  geom_boxplot(aes(color = Dose), width = 0.5, size = 0.4) +
  geom_dotplot(
    aes(fill = Dose),
    trim = FALSE,
    binaxis = 'y',
    stackdir = 'center'
  ) +
  scale_colour_manual(name = "Dose", values = doseCols) +
  scale_fill_manual(name = "Dose", values = doseCols) +
  theme_bw(base_size = 18) +
  facet_wrap( ~ Age) +
  ggtitle('Experiment 3, distribution of logViralTPM values')
ggsave(
  paste0(outdir, '/exp3_TPM.pdf'),
  plot = last_plot(),
  height = 8,
  width = 12,
  units = 'in'
)
