# Script for pathway analysis for Bauer project b830

# Load packages
library(DESeq2)
library(ggplot2)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(fgsea)

# Excel files
exp1_variant <- 'babs/differential/results/v0.1.1/differential_analyse_experiment1_variant.xlsx'
exp1_dose <- 'babs/differential/results/v0.1.1/differential_analyse_experiment1_dose.xlsx'
exp2_variant <- 'babs/differential/results/v0.1.1/differential_analyse_experiment2_variant.xlsx'

# Load comparisons
comparisons <- list(
  BA1_PBS = read_xlsx(exp1_variant, sheet = 'C1.BA.1 - PBS'),
  BA5_PBS = read_xlsx(exp1_variant, sheet = 'C1.BA.5 - PBS'),
  BA5_BA1 = read_xlsx(exp1_variant, sheet = 'C1.BA.5 - BA.1'),
  Med_Low_BA1 = read_xlsx(exp1_dose, sheet = 'A'),
  Med_Low_BA5 = read_xlsx(exp1_dose, sheet = 'B'),
  Med_Low = read_xlsx(exp1_dose, sheet = 'C'),
  BA286_PBS = read_xlsx(exp2_variant, sheet = 'C4.BA.2.86 - PBS'),
  XBB_PBS = read_xlsx(exp2_variant, sheet = 'C4.XBB - PBS'),
  XBB_BA286 = read_xlsx(exp2_variant, sheet = 'C4.XBB - BA.2.86')
)

# Load GMT
gmt <- gmtPathways('R_pathway_analysis/m5.go.bp.v2024.1.Mm.symbols.gmt')

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
  
  # Plot top results if asked for
  if (plotTop == T) {
    
    topPathwaysUp <- res[ES > 0][head(order(pval), n = 10), pathway]
    topPathwaysDown <- res[ES < 0][head(order(pval), n = 10), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    p <- plotGseaTable(pathways[topPathways], ranks, res, gseaParam = 0.5)
    pdf(file = paste0(saveDir, compName, '.pdf'),
        height = 5,
        width = 15)
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
    
    write.csv(res[,-c('leadingEdge')], file = paste0(saveDir, compName, '.csv'))
    
  }
  
  return(res)

}

# Make a list for results to go in
results <- list()

# Run through all comparisons
for (x in names(comparisons)) {
  message(paste('Starting for comparison:', x))
  results[[x]] <- runGSEA(
    DE_Table = comparisons[[x]],
    pathways = gmt,
    min = 10,
    max = 500,
    saveRes = T,
    plotTop = T,
    saveDir = 'R_pathway_analysis/',
    compName = x
  )
}




