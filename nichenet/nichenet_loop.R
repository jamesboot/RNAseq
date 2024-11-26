# Script for performing nichenet analysis on bulk RNAseq data
# Adapted from Pro's script
# 08/10/24

# Install packages if required
#devtools::install_github("saeyslab/nichenetr")
#devtools::install_github("jokergoo/ComplexHeatmap")

# Load packages
library(nichenetr)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Prep - load everything according to nichenet vignette
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md
organism <- "mouse"

if (organism == "human") {
  lr_network <-
    readRDS(url(
      "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
    ))
  ligand_target_matrix <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
      )
    )
  weighted_networks <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"
      )
    )
} else if (organism == "mouse") {
  lr_network <-
    readRDS(url(
      "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"
    ))
  ligand_target_matrix <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
      )
    )
  weighted_networks <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"
      )
    )
}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

# Target genes in rows, ligands in columns
ligand_target_matrix[1:5, 1:5]

# Interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$lr_sig) 

# Interactions and their weights in the gene regulatory network
head(weighted_networks$gr) 

# Load expression matrix
expression <- read.csv('log_vst_symbols.csv', row.names = NULL)

# Load DGE analysis results
#dge <- readxl::read_xlsx('results/latest/differential_GRCm39_analyse_treatment.xlsx', sheet = 'C1.coculture - mono|MOC1|DMSO')

# Remove NA gene symbols
expression <- expression[!is.na(expression$symbol), ]

# Remove unwanted cols, group by symbol and summarise using mean
expression <- expression[, c(4:28)] %>%
  group_by(symbol) %>%
  summarise_all(mean) %>%
  column_to_rownames(var = 'symbol')

# Transpose
expression <- t(expression)

# Load experiment table (meta)
meta <- read.csv('extdata/metadata.csv')

# Modify ID column to match expression matrix
meta$ID <- gsub('GIA7103', '', meta$ID)

# Update gene symbols
# If this is not done, there will be 35 genes fewer in lr_network_expressed!
colnames(expression) <-
  convert_alias_to_symbols(colnames(expression), "mouse", verbose = T)

# Source custom functions
source('/nemo/stp/babs/working/bootj/github/RNAseq/nichenet/geneSelection.R')
source('/nemo/stp/babs/working/bootj/github/RNAseq/nichenet/performBulkNicheNet.R')

# Run through all comparisons:
# Cancer cells co-culture to fibroblasts co-culture (all DMSO)
# Fibroblasts co-culture to cancer cells co-culture (all DMSO)
# Cancer cells co-culture to fibroblasts co-culture (all Nintedanib)
# Fibroblasts co-culture to cancer cells co-culture (all Nintedanib)
# Cancer cells mono-culture to fibroblasts mono-culture (all DMSO)
# Fibroblasts mono-culture to cancer cells mono-culture (all DMSO)
# Cancer cells mono-culture to fibroblasts mono-culture (all Nintedanib)
# Fibroblasts mono-culture to cancer cells mono-culture (all Nintedanib)

# Make a list of all comparisons
culture <- unique(meta$culture)
treatment <- unique(meta$treatment)
sender <- unique(meta$cell_type)
receiver <- unique(meta$cell_type)
compTable <- crossing(culture,treatment,sender,receiver)

# Save comparison table for reference 
write.csv(compTable, 'nichenet_compTable.csv')

# Remove comparisons between the same cell type
compTable <- compTable %>%
  mutate(Remove = case_when(sender == receiver ~ 'Y',
                            sender != receiver ~ 'N')) %>%
  filter(Remove == 'N')

# Load gene sets of interest
fibroReceiveSig <-
  unlist(as.vector(
    read.delim(
      'gene_sig_to_test_for_fibroblasts_as_receivers.txt',
      header = F
    )
  ))

cancerReceiveSig <-
  unlist(as.vector(
    read.delim(
      'gene_sig_to_test_for_cancercell_as_receivers.txt',
      header = F
    )
  ))

# Run through all comparisons in table
for (x in 1:nrow(compTable)) {
  # Get the Sender IDs and the Receiver IDs
  senderID <- meta %>%
    filter(
      culture == compTable$culture[x],
      treatment == compTable$treatment[x],
      cell_type == compTable$sender[x]
    ) %>%
    pull(ID)
  
  receiverID <- meta %>%
    filter(
      culture == compTable$culture[x],
      treatment == compTable$treatment[x],
      cell_type == compTable$receiver[x]
    ) %>%
    pull(ID)
  
  # Try Strat1 - top 2000 genes in each group (Sender or Receiver)
  Strat1GeneSets <- geneSelection(
    ExpressionMatrix = expression,
    Strategy = 'Strat1',
    Sender = senderID,
    Receiver = receiverID,
    SamplesPerGroup = 3
  )
  # Run NicheNetAnalysis - use different gene set of interest depending on receiver
  # If cell type is MOCAF1
  if (unique(meta$cell_type[meta$ID == receiverID]) == 'MOCAF1_M1') {
    performBulkNicheNet(
      GeneSets = Strat1GeneSets,
      Strategy = 'Strat1',
      gmt = fibroReceiveSig,
      SenderName = compTable$sender[x],
      ReceiverName = compTable$receiver[x],
      SenderSamples = senderID,
      ReceiverSamples = receiverID,
      expressionMat = expression,
      CompName = paste(
        compTable$culture[x],
        compTable$treatment[x],
        compTable$sender[x],
        'to',
        compTable$receiver[x],
        sep = '_'
      ),
      geneSetName = 'fibroReceiveSig'
    )
  } else if (unique(meta$cell_type[meta$ID == receiverID]) == 'MOC1') {
    performBulkNicheNet(
      GeneSets = Strat1GeneSets,
      Strategy = 'Strat1',
      gmt = cancerReceiveSig,
      SenderName = compTable$sender[x],
      ReceiverName = compTable$receiver[x],
      SenderSamples = senderID,
      ReceiverSamples = receiverID,
      expressionMat = expression,
      CompName = paste(
        compTable$culture[x],
        compTable$treatment[x],
        compTable$sender[x],
        'to',
        compTable$receiver[x],
        sep = '_'
      ),
      geneSetName = 'cancerReceiveSig'
    )
  }
}