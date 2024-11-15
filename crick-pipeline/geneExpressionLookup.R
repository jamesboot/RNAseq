# Function to plot expression of genes of interest from deseq2 objects output from rnaseq pipeline
# exp = Define from what experiment we want to look at
# model = Which model we want to look at
# variable = Variable you want to group samples by
# dat = loaded deseq2 object (analysed) from pipeline
# goi = gene symbol of gene of interest
# save.dir = location to save the output

geneExpressionLookup <- function(dat, exp, model, variable, goi, save.dir) {
  
  # Extract vst counts matrix
  vst <- as.data.frame(t(dat[[exp]][[model]]$dds@assays@data$vst))
  
  # Extract gene symbols and IDs
  genesLookup <- dat[[exp]][[model]]$dds@rowRanges@elementMetadata@listData[["symbol"]]
  
  # If gene not found return error
  if (!goi %in% genesLookup){
    return(message('ERROR: Gene not found.'))
  }
  
  # Extract meta
  colDat <- as.data.frame(dat[[exp]][[model]]$dds@colData)
  meta <- data.frame(row.names = row.names(colDat), Variable = colDat[, variable])
  
  # Merge meta and vst
  vst_meta <- merge(vst, meta, by = 'row.names')
  vst_meta <- column_to_rownames(vst_meta, var = 'Row.names')
  
  # Find the gene ID for a given symbol
  symbol <- names(genesLookup)[genesLookup %in% goi]
  
  # Filter counts down to goi
  selection <- vst_meta[, c(symbol, "Variable")]
  
  # Plot
  p <- ggplot(data = selection, aes(x = Variable, y = .data[[symbol]], colour = Variable)) +
    geom_boxplot(outliers = F) +
    labs(y = paste(goi, 'Expression log(vst)'), x = variable) +
    guides(colour = "none") +
    geom_point(size = 5, alpha = 1) +
    theme_bw(base_size = 18)
  
  # Save
  ggsave(
    filename = paste0(save.dir, '/', model, '_', goi, '.pdf'),
    units = 'in',
    height = 5,
    width = 5
  )
}
