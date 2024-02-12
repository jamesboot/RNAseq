# Script for the analysis of Yung-Yao's bulk RNAseq data 

# Setwd
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/RNAseq_Analysis/')

# Load packages
library(edgeR)
library(tidyverse)
library(gplots)
library(biomaRt)

# Find files
inFiles <- list.files(path = 'STAR',
                      pattern = '*_ReadsPerGene.out.tab',
                      full.names = T)

# Get sample names from file names
samples <- gsub(
  '_ReadsPerGene.out.tab',
  '',
  list.files(path = 'STAR',
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
GSE159_SRA <- read.delim('SRA/GSE159_SraRunTable.txt',
                         sep = ',')[, c('Run', 'Sample.Name')]
GSE159_SampNames <- read.csv('SRA/GSE159_SampNames.csv',
                             header = F,
                             col.names = c(
                               'SRA.Sample',
                               'Real.Sample'
                             ))

GSE189_SRA <- read.delim('SRA/GSE189_SraRunTable.txt',
                         sep = ',')[, c('Run', 'Sample.Name')]
GSE189_SampNames <- read.csv('SRA/GSE189_SampNames.csv',
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

# Pre-processing done - ready to go 

# Isolate just samples of interest
# Find sample names of D0 samples we need
R3381X_D0_Samps <- c("A_ctrl_0", "A_mut_0", "B_ctrl_0", "B_mut_0", "C_ctrl_0", "C_mut_0")
K2957fs_D0_Samps <- SampleNameLookup$Real.Sample[grep('Day 0', SampleNameLookup$Real.Sample)]

# Subset counts matrix to D0 samples
countMatD0 <- countMat[, c(R3381X_D0_Samps, K2957fs_D0_Samps)] 

# Create a meta data dataframe
metaD0 <- data.frame('Sample' = colnames(countMatD0),
                     'Cell.Line' = c(rep('R3381X', length(R3381X_D0_Samps)),
                                     rep('K2957fs', length(K2957fs_D0_Samps))),
                     'Group' = c('CORR', 'DMD',
                                 'CORR', 'DMD',
                                 'CORR', 'DMD',
                                 rep('DMD', 4),
                                 rep('CORR', 4)))

# Add Cell.Line+Group
metaD0 <- metaD0 %>%
  mutate(Cell.Group = paste0(Group, '_', Cell.Line))

# Make DGEList object for downstream
y <- DGEList(counts = countMatD0,
             group = metaD0$Cell.Group)

# Filter y - removes lowly expressed genes
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]

# Normalisation
y <- calcNormFactors(y)

# Extract CPM for PCA etc
NormCountsD0 <- cpm(y)
write.csv(NormCountsD0, 'NormCountsD0.csv')

# Subset down to most variable genes
# Calculate variability of every gene
geneVars <- apply(NormCountsD0, MARGIN = 1, var)
# Sort in decreasing order
geneVars <- sort(geneVars, decreasing = T)
topVarGenes <- names(geneVars[1:2000])

# PCA - rows need to be samples to transpose counts matrix
pcaMat <- t(NormCountsD0)[, topVarGenes]
samplePCA <- prcomp(pcaMat)

# Get PC Eigenvalues for later
pc_eigenvalues <- samplePCA$sdev^2

# Plot PCA
# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- samplePCA$x

# Convert to a tibble retaining the sample names as a new column
pc_scores <- pc_scores %>% 
  as_tibble(rownames = "sample")

# Add other meta data for plotting
pc_scores$`Cell Line` <- metaD0$Cell.Line
pc_scores$Group <- metaD0$Group

# Create the plot
pc_scores %>% 
  ggplot(aes(x = PC1, 
             y = PC2, 
             col = `Cell Line`, 
             shape = Group)) +
  geom_point(size = 4) +
  theme_classic()
ggsave(plot = last_plot(),
       filename = 'PCA.tiff',
       width = 6, height = 6, dpi = 300)

# Make heatmap dendrogram using most variable genes
library(ComplexHeatmap)

# Find most variable genes in normDat
vars <- apply(NormCountsD0, 1, var)

# Sort by variability and select top 'N''
topN <- 2000
varGenes <- names(sort(vars, decreasing = T))[1:topN]

# Subset normDat to varGenes only and convert to matrix
normDat_Subset <- as.matrix(NormCountsD0[varGenes, ])

# Calculate Z-Scores for genes in data set
# First calculate mean of each gene 
geneMeans <- apply(normDat_Subset, 1, mean)

# Second calculate stdev of each gene
geneStdev <- apply(normDat_Subset, 1, sd)

# Now in for loop convert to Z-Scores
zscores <- normDat_Subset
for (x in 1:topN) {
  zscores[x, ] <- (zscores[x, ] - geneMeans[x])/geneStdev[x]
}

# Annotation for heatmap
anno1 <- HeatmapAnnotation(
  Genotype = metaD0$Cell.Line,
  Group = metaD0$Group,
  col = list(
    Genotype = c('R3381X' = 'black',
                 'K2957fs' = 'grey'),
    Group = c("CORR" = 'orange',
              'DMD' = 'forestgreen')
  )
)

# Heatmap
tiff(filename = paste0('Complex_Heatmap_', topN, '.tiff'),
     height = 8,
     width = 8,
     units = 'in',
     res = 300)
Heatmap(zscores,
        clustering_distance_rows = 'pearson',
        show_row_names = F,
        show_column_names = F,
        show_row_dend = F,
        top_annotation = anno1)
dev.off()

# Make heatmap of selected genes of interest
# Normalised counts is contained within: NormCountsD0
goi <- c('ENSG00000112062', 'ENSG00000185386', 'ENSG00000188130','ENSG00000156711')
rowlabs <- c('MAPK14','MAPK11','MAPK12','MAPK13')

# Subset to GOI
GOI_NormCounts <- NormCountsD0[goi,]

# Calculate Z-Scores for genes in data set
# First calculate mean of each gene 
geneMeans <- apply(GOI_NormCounts, 1, mean)

# Second calculate stdev of each gene
geneStdev <- apply(GOI_NormCounts, 1, sd)

# Now in for loop convert to Z-Scores
zscores <- GOI_NormCounts
for (x in 1:nrow(GOI_NormCounts)) {
  zscores[x, ] <- (zscores[x, ] - geneMeans[x])/geneStdev[x]
}

# Annotation for heatmap
anno1 <- HeatmapAnnotation(
  Genotype = metaD0$Cell.Line,
  Group = metaD0$Group,
  col = list(
    Genotype = c('R3381X' = 'black',
                 'K2957fs' = 'grey'),
    Group = c("CORR" = 'orange',
              'DMD' = 'forestgreen')
  )
)

# Heatmap
tiff(filename = paste0('Complex_Heatmap_GOI.tiff'),
     height = 8,
     width = 8,
     units = 'in',
     res = 300)
Heatmap(zscores,
        clustering_distance_rows = 'pearson',
        row_labels = rowlabs,
        show_row_names = T,
        show_column_names = F,
        show_row_dend = F,
        top_annotation = anno1)
dev.off()

# Differential analysis performed from here:
# Count matrix of interest is countMatD0
# Meta data is metaD0
row.names(metaD0) <- metaD0$Sample
metaD0$Sample <- NULL

# Create DGEList object
dge <- DGEList(counts = countMatD0, samples = metaD0, group = metaD0$Group)
dge$samples

# Filter lowly expressed genes 
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalisation
dge <- calcNormFactors(dge)
dge$samples

# Exploratory plot
plotMDS(dge)

# Create design matrix
design_df <- data.frame(Sample = row.names(metaD0),
                        Group = metaD0$Group,
                        Cell_Line = metaD0$Cell.Line)

design <- model.matrix(~Cell_Line+Group,
                       data = design_df)

row.names(design) <- design_df$Sample

# Estimate dispersion
dge <- estimateDisp(dge, design, robust=TRUE)
dge$common.dispersion
plotBCV(dge)

# Fit glm
fit <- glmFit(dge, design)

# Calculate differentially expressed genes
lrt <- glmLRT(fit)
topTags(lrt)
sig.res <- topTags(lrt, 
                   adjust.method = 'fdr', 
                   p.value = 0.05,
                   sort.by = 'PValue',
                   n = 100000)
DEGs <- row.names(sig.res$table)

# Add gene names to the differential expression result table
# Get the gene names using biomaRt
ensembl <- useEnsembl('genes')
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
attributes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                    filters = 'ensembl_gene_id',
                    values = rownames(countMat),
                    mart = ensembl)
colnames(attributes) <- c('GeneID', 'GeneName')

# Extract DEG result table and add gene names using biomaRt object above
degResTable <- sig.res$table
degResTable$GeneID <- row.names(degResTable)
degResTable <- merge(degResTable, attributes, by = 'GeneID')

# Export DEGs
write.csv(degResTable, file = 'DMD_v_CORR_DEGs.csv')

# Have a look at the expression of the top DEGs
o <- order(lrt$table$PValue)
cpm(dge)[o[1:10],]

# Report DEGs at 5% FDR
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

# Extract CPM for significant DEGs for plotting
DEGs_counts <- cpm(dge)[DEGs,]

# Calculate Z-Scores for genes in data set
# First calculate mean of each gene 
geneMeans <- apply(DEGs_counts, 1, mean)

# Second calculate stdev of each gene
geneStdev <- apply(DEGs_counts, 1, sd)

# Now in for loop convert to Z-Scores
zscores <- DEGs_counts
for (x in 1:nrow(zscores)) {
  zscores[x, ] <- (zscores[x, ] - geneMeans[x])/geneStdev[x]
}

# Annotation for heatmap
anno1 <- HeatmapAnnotation(
  Genotype = metaD0$Cell.Line,
  Group = metaD0$Group,
  col = list(
    Genotype = c('R3381X' = 'black',
                 'K2957fs' = 'grey'),
    Group = c("CORR" = 'orange',
              'DMD' = 'forestgreen')
  )
)

# Heatmap
tiff(filename = paste0('Complex_Heatmap_DEGs.tiff'),
     height = 8,
     width = 8,
     units = 'in',
     res = 300)
pdf(file = paste0('Complex_Heatmap_DEGs.pdf'),
     height = 8,
     width = 8)
Heatmap(zscores,
        clustering_distance_rows = 'pearson',
        show_row_names = F,
        show_column_names = F,
        show_row_dend = F,
        top_annotation = anno1)
dev.off()

# Volcano plot of DEGs
library(EnhancedVolcano)
degResTable <- read.csv('R_outs/DMD_v_CORR_DEGs.csv')
EnhancedVolcano(degResTable,
                lab = degResTable$GeneName,
                x = 'logFC',
                y = 'FDR',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3,
                drawConnectors = T)
ggsave(plot = last_plot(),
       filename = 'DMD_v_CORR_Volcano.tiff',
       width = 6, height = 6, dpi = 300, units = 'in')

# Check DEGs for genes of interest
# Check in unfiltered results 
unfilt.res <- topTags(lrt, 
                   adjust.method = 'fdr', 
                   p.value = 1,
                   sort.by = 'PValue',
                   n = 100000)

# Add gene names to the differential expression result table
# Get the gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
attributes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                    filters = 'ensembl_gene_id',
                    values = rownames(countMat),
                    mart = ensembl)
colnames(attributes) <- c('GeneID', 'GeneName')

# Extract DEG result table and add gene names using biomaRt object above
unfilt.res <- unfilt.res$table
unfilt.res$GeneID <- row.names(unfilt.res)
unfilt.res <- merge(unfilt.res, attributes, by = 'GeneID')

# Now search for goi
goi <- c('PAX7', 'MYF5', 'MYOD1', 'MYOG', 'MYF6', 'MEF2C', 'ERBB3',
         'SPRY1', 'DLL1', 'JAG2', 'HES6', 'HEYL', 'MAPK11', 'MAPK13',
         'MAPK12', 'JAK2', 'CDKN2A')
degsGOI <- unfilt.res[unfilt.res$GeneName %in% goi, ]
write.csv(degsGOI, file = 'DMD_v_CORR_GOIs.csv')

EnhancedVolcano(degResTable,
                lab = degResTable$GeneName,
                selectLab = goi,
                x = 'logFC',
                y = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3,
                drawConnectors = T,
                max.overlaps = 25)
ggsave(plot = last_plot(),
       filename = 'DMD_v_CORR_GOI_Volcano.tiff',
       width = 6, height = 6, dpi = 300, units = 'in')
ggsave(plot = last_plot(),
       device = 'pdf',
       filename = 'DMD_v_CORR_GOI_Volcano.pdf',
       width = 6, height = 6, dpi = 300, units = 'in')

# Pathway analysis to start here
library(fgsea)

# Use unfiltered results!
# Sort result by fold change (desc)
# Get order
order <- order(unfilt.res$logFC,
               decreasing = T)

# Set order
unfilt.res <- unfilt.res[order,]

# Create list of ranks
ranks <- as.numeric(unfilt.res$logFC)
names(ranks) <- unfilt.res$GeneName
ranks <- na.omit(ranks)

# Load pathways
base.dir <- '/Users/jamesboot/Documents/9.Genome Centre Files/GC-MV-10494/GC-MV-10495_Analysis_Output3'
HallmarkPathways <- gmtPathways(file.path(base.dir, 'h.all.v2023.1.Hs.symbols.gmt.txt'))
CuratedPathways <- gmtPathways(file.path(base.dir, 'c2.cp.v2023.1.Hs.symbols.gmt.txt'))
GOPathways <- gmtPathways(file.path(base.dir, 'c5.go.bp.v2023.1.Hs.symbols.gmt.txt'))

pathways <- list(H = HallmarkPathways,
                 C2 = CuratedPathways,
                 C5 = GOPathways)

# Nested for loop to run GSEA on all comparisons and all pathway sets and save results
for (pathwaySet in names(pathways)) {
  print(paste('Running pathway analysis for set:', pathwaySet))
    # Run GSEA
    fgseaRes <- fgsea(
      pathways = pathways[[pathwaySet]],
      stats    = ranks,
      minSize  = 15,
      maxSize  = 500
    )
    
    # See top results: sort by pval
    ResOrder <- fgseaRes[order(pval), ]
    
    # Make new column to unlist gene names into 
    ResOrder$LeadingEdge <- NA
    
    # For loop to go through and unlist each element and add as new element in col
    for (x in 1:nrow(ResOrder)) {
      listElement <- ResOrder$leadingEdge[[x]]
      test <- paste(listElement, collapse = ', ')
      ResOrder$LeadingEdge[x] <- test
    }

    print(paste('Saving result...'))
    write.csv(ResOrder[, c(1:7, 9)],
              file = paste0('DMD_v_CORR_', pathwaySet, '.csv'))
    print(paste('Iteration end.'))
}

# Plot the top pathways in a table 
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
# topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
tiff('GSEA_Table.tiff',
     height = 12,
     width = 12,
     units = 'in',
     res = 300)
plotGseaTable(GOPathways[topPathwaysUp], ranks[[1]], fgseaRes, 
              gseaParam=1)
dev.off()

# Perform differential expression separately for each line ----

# First for R3 Line
# Differential analysis performed from here:
# Count matrix of interest is countMatD0

# Subset counts and meta down to cell line R3381X
R3_D0_meta <- metaD0[metaD0$Cell.Line == 'R3381X', ] 
R3_D0_count <- countMatD0[, row.names(R3_D0_meta)]

# Create DGEList object
dge <- DGEList(counts = R3_D0_count, samples = R3_D0_meta, group = R3_D0_meta$Group)
dge$samples

# Filter lowly expressed genes 
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalisation
dge <- calcNormFactors(dge)
dge$samples

# Exploratory plot
plotMDS(dge)

# Create design matrix
design_df <- data.frame(Sample = row.names(R3_D0_meta),
                        Group = R3_D0_meta$Group)

design <- model.matrix(~Group,
                       data = design_df)

row.names(design) <- design_df$Sample

# Estimate dispersion
dge <- estimateDisp(dge, design, robust=TRUE)
dge$common.dispersion
plotBCV(dge)

# Fit glm
fit <- glmFit(dge, design)

# Calculate differentially expressed genes
lrt <- glmLRT(fit)
topTags(lrt)
sig.res <- topTags(lrt, 
                   adjust.method = 'fdr', 
                   p.value = 0.05,
                   sort.by = 'PValue',
                   n = 100000)
DEGs <- row.names(sig.res$table)

# Add gene names to the differential expression result table
# Get the gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
attributes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                    filters = 'ensembl_gene_id',
                    values = rownames(countMat),
                    mart = ensembl)
colnames(attributes) <- c('GeneID', 'GeneName')

# Extract DEG result table and add gene names using biomaRt object above
degResTable <- sig.res$table
degResTable$GeneID <- row.names(degResTable)
degResTable <- merge(degResTable, attributes, by = 'GeneID')

# Export DEGs
write.csv(degResTable, file = 'R3_DMD_v_CORR_DEGs.csv')


# Second for K2 Line
# Differential analysis performed from here:
# Count matrix of interest is countMatD0
# Meta data is metaD0

# Subset counts and meta down to cell line K2957fs
K2_D0_meta <- metaD0[metaD0$Cell.Line == 'K2957fs', ] 
K2_D0_count <- countMatD0[, row.names(K2_D0_meta)]

# Create DGEList object
dge <- DGEList(counts = K2_D0_count, samples = K2_D0_meta, group = K2_D0_meta$Group)
dge$samples

# Filter lowly expressed genes 
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalisation
dge <- calcNormFactors(dge)
dge$samples

# Exploratory plot
plotMDS(dge)

# Create design matrix
design_df <- data.frame(Sample = row.names(K2_D0_meta),
                        Group = K2_D0_meta$Group)

design <- model.matrix(~Group,
                       data = design_df)

row.names(design) <- design_df$Sample

# Estimate dispersion
dge <- estimateDisp(dge, design, robust=TRUE)
dge$common.dispersion
plotBCV(dge)

# Fit glm
fit <- glmFit(dge, design)

# Calculate differentially expressed genes
lrt <- glmLRT(fit)
topTags(lrt)
sig.res <- topTags(lrt, 
                   adjust.method = 'fdr', 
                   p.value = 0.05,
                   sort.by = 'PValue',
                   n = 100000)
DEGs <- row.names(sig.res$table)

# Add gene names to the differential expression result table
# Get the gene names using biomaRt
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
attributes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                    filters = 'ensembl_gene_id',
                    values = rownames(countMat),
                    mart = ensembl)
colnames(attributes) <- c('GeneID', 'GeneName')

# Extract DEG result table and add gene names using biomaRt object above
degResTable <- sig.res$table
degResTable$GeneID <- row.names(degResTable)
degResTable <- merge(degResTable, attributes, by = 'GeneID')

# Export DEGs
write.csv(degResTable, file = 'K2_DMD_v_CORR_DEGs.csv')





