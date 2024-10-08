# Script for pathway analysis of bulk RNAseq data
# Differential expression results required
# Run downloaded txt file through FilePreProcess.ipynb first to generate CSV

# Setwd
setwd('/Users/jamesboot/Documents/9.Genome Centre Files/GC-MI-10822/Pathway_Analysis/')

# Load packages
library(tidyverse)
library(ggplot2)
library(fgsea)
library(multienrichjam)
library(enrichplot)

# Read in DE results generated from the python script
inFiles <- list.files(path = '.',
                      pattern = '.csv')

# Contains 3 sets of results
resList <- lapply(inFiles, function(x) {
  read.csv(x)
})

# Name the list elemments by comparison name
names(resList) <- c('CVID_v_CVIDAutoimmunity',
                    'CVID_v_NonCVID',
                    'CVID_v_Healthy')

# Change colnames of each list element
newColNames <- c(
  'RowNum',
  'GeneID',
  'GeneName',
  "Pvalue",
  "FDR",
  "Ratio",
  "FoldChange",
  "LSMeanGroup1",
  "LSMeanGroup2"
)
resList <- lapply(resList, setNames, newColNames)

# Load pathways
pathwayFiles <- list.files(path = '.', pattern = '.gmt.txt')
pathways <- lapply(pathwayFiles, function(x) {
  gmtPathways(x)
})

names(pathways) <- c('Curated', 'GO_BP', 'Hallmark')

# Create folder for results to go in
if (!dir.exists('Results')) {
  dir.create('Results')
} 

# Now loop through comparisons and perform analysis
for (x in names(resList)) {
  
  # Log
  message(paste('Starting analysis for comparison:', x))
  
  # Use unfiltered results
  # Sort result by fold change (desc)
  # Get order
  order <- order(resList[[x]]$FoldChange,
                 decreasing = T)
  
  # Set order
  resList[[x]] <- resList[[x]][order, ]
  
  # Create list of ranks
  ranks <- resList[[x]]$FoldChange
  names(ranks) <- resList[[x]]$GeneName
  ranks <- na.omit(ranks)
  
  # Nested for loop to run GSEA on all comparisons and all pathway sets and save results
  for (pathwaySet in names(pathways)) {
    
    # Log
    message(paste('Running analysis for pathway set:', pathwaySet))
    
    # Run GSEA
    fgseaRes <- fgsea(
      pathways = pathways[[pathwaySet]],
      stats    = ranks,
      minSize  = 15,
      maxSize  = 500
    )
    
    # See top results: sort by pval
    ResOrder <- fgseaRes[order(pval),]
    
    # Make new column to unlist gene names into
    ResOrder$LeadingEdge <- NA
    
    # For loop to go through and unlist each element and add as new element in col
    for (i in 1:nrow(ResOrder)) {
      listElement <- ResOrder$leadingEdge[[i]]
      test <- paste(listElement, collapse = ', ')
      ResOrder$LeadingEdge[i] <- test
    }
    
    # Log and save 
    message(paste('Saving result for:', pathwaySet))
    write.csv(ResOrder[, c(1:7, 9)],
              file = paste0('Results/', x, '_', pathwaySet, '.csv'))
    
    ResOrder$leadingEdge <- NULL
    
    ResOrder <- as.data.frame(ResOrder)
    
    # VISUALISATION
    
    # Ensure order of pathway is determined by padj
    ResOrder$pathway <- factor(ResOrder$pathway, 
                          levels = ResOrder$pathway[order(ResOrder$padj, decreasing = TRUE)])
    
    # Bar plot
    p <- ggplot(data = ResOrder[c(1:20),], aes(x = pathway, y = size, fill = padj)) +
      geom_bar(stat = "identity") +
      theme(axis.text = element_text(size = 6)) +
      scale_fill_continuous(low="blue", high="red") +
      coord_flip()
    ggsave(p,
           filename = paste0('Results/', x, '_', pathwaySet, '_bar.tiff'),
           height = 8,
           width = 8,
           dpi = 300,
           units = 'in')
    
    # Add list of pathway genes to res
    # Add true pathway size, fgsea size is a filtered size
    ResOrder$pathwayGenes <- NA
    ResOrder$TrueSize <- NA
    for (i in 1:length(ResOrder$pathway)){
      ResOrder$pathwayGenes[i] <- paste(pathways[[pathwaySet]][[as.character(ResOrder$pathway[i])]], collapse = ',')
      ResOrder$TrueSize[i] <- length(pathways[[pathwaySet]][[as.character(ResOrder$pathway[i])]])
    }
    
    # Convert dataframe to enrichment result
    enrichRes <- enrichDF2enrichResult(ResOrder,
                                       pvalueCutoff = 1,
                                       pAdjustMethod = 'none',
                                       keyColname = 'pathway',
                                       pathGenes = 'TrueSize',
                                       geneColname = 'pathwayGenes',
                                       geneDelim = ',',
                                       pvalueColname = 'padj')
    
    # Create similarity matrix 
    enrichRes2 <- pairwise_termsim(enrichRes)
    
    # Enrichment map
    map1 <- emapplot(enrichRes2,
                    cex.params = list(category_label = 0.7))
    ggsave(map1,
           filename = paste0('Results/', x, '_', pathwaySet, '_map.tiff'),
           height = 16,
           width = 16,
           dpi = 300,
           units = 'in')
    
  }
  # Log
  message(paste('Finished analysis for:', x))
}
