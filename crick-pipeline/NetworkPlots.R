# Create network plots for enrichment results

# Install package
#devtools::install_github('ievaKer/aPEAR')

# Load packages
library(aPEAR)
library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(fgsea)

# Load DEGs
# Excel files
exp1_variant <- 'babs/differential/results/latest/differential_analyse_experiment1_variant.xlsx'
exp1_dose <- 'babs/differential/results/latest/differential_analyse_experiment1_dose.xlsx'
exp2_variant <- 'babs/differential/results/latest/differential_analyse_experiment2_variant.xlsx'
exp3 <- 'babs/differential/results/latest/differential_analyse_experiment3.xlsx'

# Load comparisons
comparisons <- list(
  BA1_PBS = read_xlsx(exp1_variant, sheet = 'A'),
  BA5_PBS = read_xlsx(exp1_variant, sheet = 'B'),
  BA5_BA1 = read_xlsx(exp1_variant, sheet = 'C'),
  Exp1_Virus_PBS = read_xlsx(exp1_variant, sheet = 'D'),
  
  Med_Low_BA1 = read_xlsx(exp1_dose, sheet = 'A'),
  Med_Low_BA5 = read_xlsx(exp1_dose, sheet = 'B'),
  Exp1_Med_Low = read_xlsx(exp1_dose, sheet = 'C'),
  
  BA286_PBS = read_xlsx(exp2_variant, sheet = 'A'),
  XBB_PBS = read_xlsx(exp2_variant, sheet = 'B'),
  XBB_BA286 = read_xlsx(exp2_variant, sheet = 'C'),
  Exp2_Virus_PBS = read_xlsx(exp2_variant, sheet = 'D'),
  
  Old_Young_Low = read_xlsx(exp3, sheet = 'A'),
  Old_Young_Hi = read_xlsx(exp3, sheet = 'B'),
  Old_Young_Med = read_xlsx(exp3, sheet = 'C'),
  Exp3_Hi_Low = read_xlsx(exp3, sheet = 'D'),
  Exp3_Med_Low = read_xlsx(exp3, sheet = 'E'),
  Exp3_Med_Hi = read_xlsx(exp3, sheet = 'F')
)

# Run through all comparisons in for loop
for (x in names(comparisons)) {
  message(paste('Starting for comparison:', x))
  
  # Get ordered gene list
  message(paste('Ordering list...'))
  genelistDF <- comparisons[[x]] %>%
    arrange(desc(shrunkLFC)) %>%
    dplyr::select(id, shrunkLFC)
  
  geneList <- genelistDF$shrunkLFC
  names(geneList) <- genelistDF$id
  
  # Run gseGO to get enrichment result
  message(paste('Performing enrichment...'))
  enrich <- gseGO(
    geneList,
    OrgDb = org.Mm.eg.db,
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 200,
    keyType = 'ENSEMBL',
    by = 'fgsea',
    eps = 0,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01
  )
  
  # If no sig results skip
  if (nrow(enrich@result) < 5) {
    next
  } else {
    # Plot
    message(paste('Plotting...'))
    p <- enrichmentNetwork(
      enrich@result,
      drawEllipses = TRUE,
      fontSize = 2.5,
      simMethod = 'jaccard',
      clustMethod = 'hier',
      repelLabels = T,
      minClusterSize = 3
    )
    ggsave(
      plot = p,
      filename = paste0('R_pathway_analysis/', x, '_network.pdf'),
      height = 15,
      width = 15,
      units = 'in'
    )
  }
  
}

# If you want an interactive plot - you canuse plotly
library(plotly)

ggplotly(p, tooltip=c('ID', 'Cluster', 'Cluster size'))
