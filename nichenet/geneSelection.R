# Functions for NicheNet analysis on bulk RNAseq

# Strategies for determining expressed genes in senders (CAFs) and receivers (Cancers):
# Strat1 = top 2000 genes
# Strat2 = all genes with vst norm counts >5
# Strat3 = ?
# Make a function for gene selection
# ExpressionMatrix is the expression matrix 
# Strategy must be either 'Strat1' or 'Strat2'
# Senders and Receivers should be a vector of sample names for each group
# SamplesPerGroup is how many replicates there are in the sender and receiver groups
# Strat2Threshold is how many vst norm counts must samples be greater than - only required for strat2
# Will return a list - element 1 is expressed genes sender, element 2 is receiver
geneSelection <- function(ExpressionMatrix,
                          Strategy,
                          Sender,
                          Receiver,
                          SamplesPerGroup,
                          Strat2Threshold) {
  if (!Strategy %in% c('Strat1', 'Strat2')) {
    return(message('ERROR: Invalid Strategy provided.'))
  } else if (sum(Sender %in% rownames(ExpressionMatrix)) != SamplesPerGroup) {
    return(message('ERROR: Invalid Sender names provided.'))
  } else if (sum(Receiver %in% rownames(ExpressionMatrix)) != SamplesPerGroup) {
    return(message('ERROR: Invalid Receiver names provided.'))
  } else if (Strategy == 'Strat1') {
    message('Applying Strategy 1...')
    # Determine which genes are expressed in senders (CAFs)
    expressed_genes_sender <-
      names(sort(apply(ExpressionMatrix[Sender,], 2, mean), decreasing = T))[1:2000]
    # Determine which genes are expressed in receivers (Cancers)
    expressed_genes_receiver <-
      names(sort(apply(ExpressionMatrix[Receiver,], 2, mean), decreasing = T))[1:2000]
  } else if (Strategy == 'Strat2') {
    message(paste0('Applying Strategy 2, with vst threshold of: ', Strat2Threshold))
    # Need a function to test whether values are greater than Strat2Threshold
    gt <- function(x) {
      x > Strat2Threshold
    }
    # Find columns for senders where vst counts are greater than 5 in all 3 samples
    colsKeep <-
      colSums(apply(ExpressionMatrix[Sender,], 2, gt)) == SamplesPerGroup
    expressed_genes_sender <-
      colnames(ExpressionMatrix[Sender,])[colsKeep]
    # Find columns for senders where vst counts are greater than 5 in all 3 samples
    colsKeep <-
      colSums(apply(ExpressionMatrix[Receiver,], 2, gt)) == SamplesPerGroup
    expressed_genes_receiver <-
      colnames(ExpressionMatrix[Receiver,])[colsKeep]
  }
  # Create list of results
  res <- list(senderGenes = expressed_genes_sender,
              receiverGenes = expressed_genes_receiver)
  # Print summary of overlaps
  nOverlap <-
    length(intersect(res$senderGenes, res$receiverGenes))
  nUnion <-
    length(union(res$senderGenes, res$receiverGenes))
  message(paste('Total unique Genes =', nUnion))
  message(paste('n Genes overlapping =', nOverlap))
  message(paste0('% Genes overlapping = ', round((nOverlap / nUnion) * 100, digits = 1), '%'))
  # Return results
  return(res)
}