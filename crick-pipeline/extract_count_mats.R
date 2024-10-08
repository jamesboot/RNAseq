# Script for extracting bits and bobs for Giovanni

# Load packages
library(DESeq2)

# Create non-log transformed matrix from DESeq2 object (note gene IDs will be used) ----
# Read in DESeq2 object
obj <- readRDS('data/analyse_x_GRCm39.rds')

# View counts matrix 
obj$treatment$model1$dds@assays@data$vst

# Extract vst matrix and raise to power 2
nonlog_vsd <- 2^obj$treatment$model1$dds@assays@data$vst
View(nonlog_vsd)

# Save
write.csv(nonlog_vsd, 'nonlog_vst_mat.csv')

# Create log(vst) matrix and non-log vst matrix from txt pipeline output ----
# Read in text file from results
test <- read.delim('results/latest/vst_GRCm39_analyse_treatment.txt',
                   skip = 32)[,c(1,2,3,28:51)]

# Save log(vst) matrix
write.csv(test, 'log_vst_symbols.csv')

# Create a function to raise to power 2 
raise2 <- function(x){2^x}

# Raise to power 2 to transform to non-log
nonlog_vst_withSymbols <- test
nonlog_vst_withSymbols[c(4:27)] <- lapply(nonlog_vst_withSymbols[c(4:27)], raise2)

# Save 
write.csv(nonlog_vst_withSymbols, 'nonlog_vst_symbols.csv')
