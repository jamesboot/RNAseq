# Define function for running GSEA
runGSEA <- function(DE_Table, pathways, min, max, saveRes, plotTop, compName, saveDir) {
  
  # Load packages
  library(DESeq2)
  library(ggplot2)
  library(readxl)
  library(dplyr)
  library(RColorBrewer)
  library(fgsea)

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
      
      # Get top upreg pathways
      topPathwaysUpRes <- res[ES > 0][head(order(padj), n = 20), ]
      # Get top downreg pathways
      topPathwaysDownRes <- res[ES < 0][head(order(padj), n = 20), ]
      topPathwaysDownRes <- topPathwaysDownRes[order(topPathwaysDownRes$NES, decreasing = F), ]
      # Bind together
      topUpDownPaths <- rbind(topPathwaysUpRes, topPathwaysDownRes)
      # Set order
      topUpDownPaths <- topUpDownPaths[order(topUpDownPaths$NES, decreasing = F), ]
      topUpDownPaths$pathway <- factor(topUpDownPaths$pathway, levels = topUpDownPaths$pathway)
      
      # Bubble plot
      pdf(
        file = paste0(saveDir, compName, '_bubble.pdf'),
        height = 5,
        width = 15
      )
      p <- ggplot(topUpDownPaths, aes(
        x = NES,
        y = pathway,
        size = size,
        fill = padj
      )) +
        geom_point(alpha = 1, shape = 21) +
        scale_fill_viridis_c(direction = -1,
                             limits = c(min(topUpDownPaths$padj), 1))
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
