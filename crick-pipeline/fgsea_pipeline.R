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

# Load comparisons
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

# Define function for running GSEA
runGSEA <- function(DE_Table, pathways, min, max, saveRes, plotTop, compName, saveDir) {
  
  # Create ranks vector
  ranks <- c(DE_Table$log2FoldChange)
  names(ranks) <- DE_Table$symbol
  
  # Omit genes with NA values
  ranks <- na.omit(ranks)
  
  # Omit genes without entrez ID
  ranks <- ranks[!is.na(names(ranks))]
  
  # Run fgsea
  res <- fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = min,
    maxSize = max
  )
  
  # Add control step for if no results are generated
  if (nrow(res) > 0) {
    
    # Plot top results if asked for
    if (plotTop == T) {
      topPathwaysUp <- res[ES > 0][head(order(pval), n = 10), pathway]
      topPathwaysDown <-
        res[ES < 0][head(order(pval), n = 10), pathway]
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
      p <-
        plotGseaTable(pathways[topPathways], ranks, res, gseaParam = 0.5)
      pdf(
        file = paste0(saveDir, compName, '.pdf'),
        height = 5,
        width = 15
      )
      plot(p)
      dev.off()

      topPathwaysUpRes <- res[ES > 0][head(order(padj), n = 20), ]
      topPathwaysUpRes <- topPathwaysUpRes[order(topPathwaysUpRes$NES, decreasing = F), ]
      topPathwaysUpRes$pathway <- factor(topPathwaysUpRes$pathway, levels = topPathwaysUpRes$pathway)
      pdf(
        file = paste0(saveDir, compName, '_bubble.pdf'),
        height = 5,
        width = 15
      )
      p <- ggplot(topPathwaysUpRes, aes(
        x = NES,
        y = pathway,
        size = size,
        color = padj
      )) +
        geom_point(alpha = 1) +
        scale_colour_gradient(low = "red", high = "blue")
      plot(p)
      dev.off()

    }
    
    if (saveRes == T) {
      # Make new column to unlist gene names into
      res$LeadingEdge <- NA
      
      # For loop to go through and unlist each element and add as new element in col
      for (i in 1:nrow(res)) {
        listElement <- res$leadingEdge[[i]]
        test <- paste(listElement, collapse = ', ')
        res$LeadingEdge[i] <- test
      }
      
      write.csv(res[, -c('leadingEdge')], file = paste0(saveDir, compName, '.csv'))
      
    }
  }
  return(res)
}

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