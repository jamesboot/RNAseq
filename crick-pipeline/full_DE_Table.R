# Function to extract the full DEG results from DESeq2 obj
full_DE_Table <- function(DESeqRDS) {
  # Libraries
  library(DESeq2)
  library(openxlsx)
  
  # Load
  obj <- readRDS(DESeqRDS)
  
  # Find all the datasets
  sets <- names(obj)
  
  # Next for loop to loop through all datasets then all models and all comparisons
  # Write xlsx book for each dataset_model combination - with a worksheet for each comparison
  for (x in sets) {
    models <- names(obj[[x]])
    for (y in models) {
      comps <- names(obj[[x]][[y]]$comps)
      wb <- buildWorkbook(sheetName = 'Key',
                          data.frame(Comparison = comps,
                                     Sheet = LETTERS[1:length(comps)]))
      for (z in 1:length(comps)) {
        df <- as.data.frame(obj[[x]][[y]]$comps[[comps[z]]])
        df$id <- rownames(df)
        addWorksheet(wb, sheetName = LETTERS[z])
        writeData(wb, sheet = LETTERS[z], df)
      }
      saveWorkbook(
        wb,
        file = paste(x, y, 'DEG_FullTable.xlsx', sep = '_'),
        overwrite = T
      )
    }
  }
}
