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
source('geneSelection.R')
source('performBulkNicheNet.R')

# 1. DMSO COCULTURE CAFS > COCULTURE CANCER

# Get the CAF IDs and the Cancer IDs
cafID <- meta %>%
  filter(culture == 'coculture',
         treatment == 'DMSO',
         cell_type == 'MOCAF1_M1') %>%
  pull(ID)

cancerID <- meta %>%
  filter(culture == 'coculture',
         treatment == 'DMSO',
         cell_type == 'MOC1') %>%
  pull(ID)

# Try Strat1 
Strat1GeneSets <- geneSelection(
  ExpressionMatrix = expression,
  Strategy = 'Strat1',
  Sender = cafID,
  Receiver = cancerID,
  SamplesPerGroup = 3
)
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat1GeneSets,
  Strategy = 'Strat1',
  gmt = 'mh.all.v2024.1.Mm.symbols.gmt',
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'coculture_CAF_to_Cancer',
  geneSetName = 'GO_Hallmark'
)
# Try Strat1 with different gene set of interest - top 2000 genes in co-culture cancer cells
topCancerGenes <- 
  names(sort(apply(expression[cancerID,], 2, mean), decreasing = T))[1:1500]
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat1GeneSets,
  Strategy = 'Strat1',
  gmt = topCancerGenes,
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'coculture_CAF_to_Cancer',
  geneSetName = 'topCancerGenes'
)

# Try Strat2
Strat2GeneSets <- geneSelection(
  ExpressionMatrix = expression,
  Strategy = 'Strat2',
  Sender = cafID,
  Receiver = cancerID,
  SamplesPerGroup = 3,
  Strat2Threshold = 6
)
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat2GeneSets,
  Strategy = 'Strat2',
  gmt = 'mh.all.v2024.1.Mm.symbols.gmt',
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'coculture_CAF_to_Cancer',
  geneSetName = 'GO_Hallmark'
)

# Try Strat2 with different gene set of interest - top 2000 genes in co-culture cancer cells
topCancerGenes <- 
  names(sort(apply(expression[cancerID,], 2, mean), decreasing = T))[1:1500]
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat2GeneSets,
  Strategy = 'Strat2',
  gmt = topCancerGenes,
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'coculture_CAF_to_Cancer',
  geneSetName = 'topCancerGenes'
)

# 2. DMSO MONOCULTURE CAFS > MONOCULTURE CANCER

# Get the CAF IDs and the Cancer IDs
cafID<- meta %>%
  filter(culture == 'mono',
         treatment == 'DMSO',
         cell_type == 'MOCAF1_M1') %>%
  pull(ID)

cancerID <- meta %>%
  filter(culture == 'mono',
         treatment == 'DMSO',
         cell_type == 'MOC1') %>%
  pull(ID)

# Try strat1 
Strat1GeneSets <- geneSelection(
  ExpressionMatrix = expression,
  Strategy = 'Strat1',
  Sender = cafID,
  Receiver = cancerID,
  SamplesPerGroup = 3
)
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat1GeneSets,
  Strategy = 'Strat1',
  gmt = 'mh.all.v2024.1.Mm.symbols.gmt',
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'monoculture_CAF_to_Cancer',
  geneSetName = 'GO_Hallmark'
)
# Try Strat1 with different gene set of interest - top 2000 genes in co-culture cancer cells
topCancerGenes <- 
  names(sort(apply(expression[cancerID,], 2, mean), decreasing = T))[1:1500]
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat1GeneSets,
  Strategy = 'Strat1',
  gmt = topCancerGenes,
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'monoculture_CAF_to_Cancer',
  geneSetName = 'topCancerGenes'
)

# Try strat2
Strat2GeneSets <- geneSelection(
  ExpressionMatrix = expression,
  Strategy = 'Strat2',
  Sender = cafID,
  Receiver = cancerID,
  SamplesPerGroup = 3,
  Strat2Threshold = 6
)
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat2GeneSets,
  Strategy = 'Strat2',
  gmt = 'mh.all.v2024.1.Mm.symbols.gmt',
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'monoculture_CAF_to_Cancer',
  geneSetName = 'GO_Hallmark'
)

# Try Strat2 with different gene set of interest - top 2000 genes in mono-culture cancer cells
topCancerGenes <- 
  names(sort(apply(expression[cancerID,], 2, mean), decreasing = T))[1:1500]
# Run NicheNetAnalysis
performBulkNicheNet(
  GeneSets = Strat2GeneSets,
  Strategy = 'Strat2',
  gmt = topCancerGenes,
  SenderName = 'CAF',
  ReceiverName = 'CancerCell',
  SenderSamples = cafID,
  ReceiverSamples = cancerID,
  expressionMat = expression,
  CompName = 'monoculture_CAF_to_Cancer',
  geneSetName = 'topCancerGenes'
)

 