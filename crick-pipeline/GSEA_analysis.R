# Script for pathway analysis for projects

# Load packages
library(DESeq2)
library(ggplot2)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(fgsea)

# Locate excel files
infile <- list.files('results/latest', pattern = '*.xlsx', full.names = T)

# Load comparison results
comparisons <- list(
  Mono_ET_v_HD = read_xlsx(infile, sheet = 'C1.ET_MSC - HD_MSC|Mono-culture'),
  Coculture_ET_v_HD = read_xlsx(infile, sheet = 'A'),
  HD_Coculture_v_Mono = read_xlsx(infile, sheet = 'B'),
  ET_Coculture_v_Mono = read_xlsx(infile, sheet = 'C')
)

# Locate GMT files
gmtfile <- list.files('.', pattern = '*.gmt', full.names = T)
gmtNames <- gsub('./', '', gsub('.gmt', '', gmtfile))

# Load GMT
gmt <- lapply(gmtfile, gmtPathways)
names(gmt) <- gmtNames

# Load custom function
source('runGSEA.R')

# Make a list for results to go in
results <- list()

# Run through all comparisons and gmt files
for (x in names(comparisons)) {
  message(paste('Starting for comparison:', x))
  for (pathSet in names(gmt)) {
    message(paste('Starting for pathway:', pathSet))
    results[[x]] <- runGSEA(
      DE_Table = comparisons[[x]],
      pathways = gmt[[pathSet]],
      min = 10,
      max = 250,
      saveRes = T,
      plotTop = T,
      saveDir = 'custom_outs/',
      compName = paste0(x, '_', pathSet)
    )
  }
}