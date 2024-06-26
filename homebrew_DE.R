# Develop home-made method for differential gene expression
# Compare results with EdgeR
# Use data from old project

# Load required packages
library(tidyverse)
library(ggplot2)
library(MASS)
library(lmtest)

# Set working directory
setwd('/Users/jamesboot/Documents/GitHubEnterprise/RNAseq/')

# Specify data and meta directory
dataDir <- '/Users/jamesboot/Documents/9.Genome Centre Files/YungYao_Analysis/STAR'
metaDir <- '/Users/jamesboot/Documents/9.Genome Centre Files/YungYao_Analysis/SRA'

# Find files
inFiles <- list.files(path = dataDir,
                      pattern = '*_ReadsPerGene.out.tab',
                      full.names = T)

# Get sample names from file names
samples <- gsub(
  '_ReadsPerGene.out.tab',
  '',
  list.files(path = dataDir,
             pattern = '*_ReadsPerGene.out.tab')
)

# Load files
rawDat <- lapply(inFiles, function(x){
  read.delim(x,
             sep = '\t',
             skip = 4,
             header = F)
})

# Rename elements of list
names(rawDat) <- samples

# Tidy up each of the dataframes in list
# Use col 4 of the dataframe, see link below for reason
# https://www.biostars.org/p/285757/
# https://www.biostars.org/p/407788/ 
for (x in samples) {
  row.names(rawDat[[x]]) <- rawDat[[x]]$V1
  rawDat[[x]] <- rawDat[[x]][, c(1,4)]
  colnames(rawDat[[x]]) <- c('GeneID', x)
}

# Merge list elements
countMat <- rawDat %>% 
  reduce(left_join, by = "GeneID")

# Make gene names row names
row.names(countMat) <- countMat$GeneID
countMat$GeneID <- NULL

# Need to update the sample names (colnames)
GSE159_SRA <- read.delim(file.path(metaDir, 'GSE159_SraRunTable.txt'),
                         sep = ',')[, c('Run', 'Sample.Name')]
GSE159_SampNames <- read.csv(file.path(metaDir, 'GSE159_SampNames.csv'),
                             header = F,
                             col.names = c(
                               'SRA.Sample',
                               'Real.Sample'
                             ))

GSE189_SRA <- read.delim(file.path(metaDir, 'GSE189_SraRunTable.txt'),
                         sep = ',')[, c('Run', 'Sample.Name')]
GSE189_SampNames <- read.csv(file.path(metaDir, 'GSE189_SampNames.csv'),
                             header = F,
                             col.names = c(
                               'SRA.Sample',
                               'Real.Sample'
                             ))

All_SRA <- rbind(GSE159_SRA, GSE189_SRA)
colnames(All_SRA) <- c('Run', 'SRA.Sample')
All_SampNames <- rbind(GSE159_SampNames, GSE189_SampNames)

SampleNameLookup <- list(All_SRA, All_SampNames) %>% 
  reduce(left_join, by = "SRA.Sample")

colnames(countMat) <- 
  plyr::mapvalues(colnames(countMat),
                  as.character(SampleNameLookup$Run), 
                  as.character(SampleNameLookup$Real.Sample))

# Subset counts matrix to samples of itnerest
soi <- c("A_ctrl_0", "A_mut_0", "B_ctrl_0", "B_mut_0", "C_ctrl_0", "C_mut_0")
CountMat <- countMat[, soi] 

# Create a meta data dataframe
meta <- data.frame('Sample' = colnames(CountMat),
                   'Group' = rep(c('CORR', 'DMD'), 3))

# Plot distribution of counts
CountMatLog <- log10(CountMat)
ggplot(CountMatLog, aes(x=A_ctrl_0)) + geom_histogram()
ggplot(CountMat, aes(x=A_ctrl_0)) + geom_histogram(binwidth = 10)

# Filter - remove genes with max counts of 10 across all samples
maximums <- apply(CountMat, MARGIN = 1, max)
filtCountMat <- CountMat[maximums >= 10, ]

# Normalisation (convert to counts per million)
# Calculate normalisation factors
normFactors <- apply(filtCountMat, MARGIN = 2, function(x){
  1e6/sum(x)
})

# Apply normalisation factors
normCountMat <- filtCountMat
for (x in 1:ncol(normCountMat)) {
  normCountMat[, x] <- normCountMat[, x]*normFactors[x]
}

# Perform differential expression with glm poisson

# Test on a single row

# Use group from meta data as grouping factor
Group <- meta$Group

# Fit models
fit0 <- glm.nb(unlist(normCountMat[1, ]) ~ 1)
fit1 <- glm.nb(unlist(normCountMat[1, ]) ~ Group)

# Perform likelihood ratio test for differences in models
lrtResult <- lrtest(fit1, fit0)
lrtResult$`Pr(>Chisq)`[2]

# Apply to whole data
prelimRes <- apply(normCountMat, 1,
      function(y) {
        # Fit null model
        fit0 <- glm.nb(y ~ 1)
        # Fit proposed model
        fit1 <- glm.nb(y ~ Group)
        # Compare fit of models
        lrtResult <- lrtest(fit1, fit0)
        # Return p-value
        return(lrtResult$`Pr(>Chisq)`[2])
      })

# Calculate adjusted p-values
results <- data.frame(Gene = names(prelimRes),
                      PVal = prelimRes,
                      FDR = p.adjust(prelimRes, method = 'fdr')
)

# Extract FDR <0.05
sigDEGs <- results[(results$FDR < 0.05), ]
sort(sigDEGs$FDR, decreasing = F)

# Look at expression of results
normCountMat[sigDEGs,]

# Now perform analysis in edegR to compare
library(edgeR)

y <- DGEList(counts=normCountMat,group=Group)
y <- calcNormFactors(y)
design <- model.matrix(~Group)
y <- estimateDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
edgeR_res <- topTags(lrt, p.value = 0.05, n = 1000)[['table']]

# Are sigDEGs all contains within EdgeR DEGs
sigDEGs %in% edgerDEGs
