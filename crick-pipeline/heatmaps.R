# Script to lookup expression of DE genes from previous study in all AF787 data/samples

# Libs
library(tidyverse)
library(readxl)
library(ComplexHeatmap)

# Load previous DE results
sheets <- excel_sheets('Previous_DE_Res.xlsx')
sheets<- sheets[1:3]
prevRes <- lapply(sheets, function(x){
  read_xlsx('Previous_DE_Res.xlsx', sheet = x)
})
names(prevRes) <- sheets

# Filter all sheets for only sign results
prevRes <- lapply(prevRes, function(x){
  x %>%
    filter(padj < 0.01)
})

# Load DESeq2 object
obj <- readRDS('data/analyse_x_GRCh38.rds')

# Edit meta data
obj$All$Dummy$dds$CellType_Culture <-
  paste(obj$All$Dummy$dds$CellType,
        obj$All$Dummy$dds$Culture,
        sep = '_')

# exp = Define from what experiment we want to look at
# model = Which model we want to look at
# variable = Variable you want to group samples by
# dat = loaded deseq2 object (analysed) from pipeline
# goi = gene symbol of gene of interest
# save.dir = location to save the output
dat = obj
exp = 'All'
model = 'Dummy'
variable = 'CellType_Culture'
save.dir = 'custom_outs'

# Create heatmap for ET samples ----

# Define gene set of interest (goi)
goi <- prevRes$DE_ET_HD %>%
  pull(gene)

# Extract vst counts matrix
vst <- as.data.frame(t(dat[[exp]][[model]]$dds@assays@data$vst))

# Extract gene symbols and IDs
genesLookup <-
  dat[[exp]][[model]]$dds@rowRanges@elementMetadata@listData[["symbol"]]

# Extract meta
colDat <- as.data.frame(dat[[exp]][[model]]$dds@colData)

# Add outlier sample annotation
outliers <- c('A12','A14','A16','A18')
colDat$Outlier[rownames(colDat) %in% outliers] <- 'True'
colDat$Outlier[!rownames(colDat) %in% outliers] <- 'False'

# Find the gene ID for a given symbol
symbol <- names(genesLookup)[genesLookup %in% goi]

# Filter counts down to goi
selection <- as.matrix(vst[, c(symbol)])

# Subset down to select samples
soi <- rownames(colDat %>%
                  filter(
                    CellType_Culture %in% c(
                      'ET_MSC_ET_HSPC_Co-culture',
                      'ET_MSC_Mono-culture',
                      'HD_MSC_ET_HSPC_Co-culture',
                      'HD_MSC_Mono-culture'
                    )
                  ))
selection <- selection[soi, ]

# Scale
selection <- scale(selection, center = T, scale = T)

# Prep annotation
# Get rows of selection and colDat in same order
colDat <- colDat[rownames(selection),]

# Get annotation
ra = HeatmapAnnotation(
  CellType_Culture = colDat$CellType_Culture,
  Outlier = colDat$Outlier,
  col = list(
    CellType_Culture = c(
      'ET_MSC_ET_HSPC_Co-culture' = '#F0E442',
      'ET_MSC_Mono-culture' = '#E69F00',
      'HD_MSC_ET_HSPC_Co-culture' = '#0072B2',
      'HD_MSC_Mono-culture' = '#56B4E9'
    ),
    Outlier = c("True" = "black", "False" = "beige")
  ),
  gp = gpar(col = "black"),
  annotation_label = NULL,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 8),
                                 title_gp = gpar(fontsize = 10, fontface = "bold"))
)

# Plot unsupervised clusters
hm <- Heatmap(
  t(selection),
  show_row_names = F,
  top_annotation = ra,
  heatmap_legend_param = list(title = "Z-Score",
                              labels_gp = gpar(fontsize = 8),
                              title_gp = gpar(fontsize = 10, fontface = "bold")),
  column_title = 'Unsupervised clustering of samples based on ET vs HD DEGs',
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8)
)

# Save
pdf(file="custom_outs/ET_HD_Unsupervised.pdf")
draw(hm)
dev.off()  

# Plot supervised clusters
hm <- Heatmap(
  t(selection),
  show_row_names = F,
  top_annotation = ra,
  heatmap_legend_param = list(title = "Z-Score",
                              labels_gp = gpar(fontsize = 8),
                              title_gp = gpar(fontsize = 10, fontface = "bold")),
  column_title = 'Supervised clustering of samples based on ET vs HD DEGs',
  column_split = colDat$CellType_Culture,
  cluster_columns = F,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8)
)

# Save
pdf(file="custom_outs/ET_HD_Supervised.pdf")
draw(hm)
dev.off() 

# Now repeat for PV samples ----

# Define gene set of interest (goi)
goi <- prevRes$DE_PV_HD %>%
  pull(gene)

# Extract vst counts matrix
vst <- as.data.frame(t(dat[[exp]][[model]]$dds@assays@data$vst))

# Extract gene symbols and IDs
genesLookup <-
  dat[[exp]][[model]]$dds@rowRanges@elementMetadata@listData[["symbol"]]

# Extract meta
colDat <- as.data.frame(dat[[exp]][[model]]$dds@colData)

# Add outlier sample annotation
outliers <- c('A12','A14','A16','A18')
colDat$Outlier[rownames(colDat) %in% outliers] <- 'True'
colDat$Outlier[!rownames(colDat) %in% outliers] <- 'False'

# Find the gene ID for a given symbol
symbol <- names(genesLookup)[genesLookup %in% goi]

# Filter counts down to goi
selection <- as.matrix(vst[, c(symbol)])

# Subset down to select samples
soi <- rownames(colDat %>%
                  filter(
                    CellType_Culture %in% c(
                      'PV_MSC_PV_HSPC_Co-culture',
                      'PV_MSC_Mono-culture',
                      'HD_MSC_PV_HSPC_Co-culture',
                      'HD_MSC_Mono-culture'
                    )
                  ))
selection <- selection[soi, ]

# Scale
selection <- scale(selection, center = T, scale = T)

# Remove NA
selection <- selection[, !colSums(is.na(selection))]

# Prep annotation
# Get rows of selection and colDat in same order
colDat <- colDat[rownames(selection),]

# Get annotation
ra = HeatmapAnnotation(
  CellType_Culture = colDat$CellType_Culture,
  Outlier = colDat$Outlier,
  col = list(
    CellType_Culture = c(
      'PV_MSC_PV_HSPC_Co-culture' = '#7AD151FF',
      'PV_MSC_Mono-culture' = '#22A884FF',
      'HD_MSC_PV_HSPC_Co-culture' = '#414487FF',
      'HD_MSC_Mono-culture' = '#440154FF'
    ),
    Outlier = c("True" = "black", "False" = "beige")
  ),
  gp = gpar(col = "black"),
  annotation_label = NULL,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 8),
                                 title_gp = gpar(fontsize = 10, fontface = "bold"))
)

# Plot unsupervised clusters
hm <- Heatmap(
  t(selection),
  show_row_names = F,
  top_annotation = ra,
  heatmap_legend_param = list(title = "Z-Score",
                              labels_gp = gpar(fontsize = 8),
                              title_gp = gpar(fontsize = 10, fontface = "bold")),
  column_title = 'Unsupervised clustering of samples based on PV vs HD DEGs',
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8)
)

# Save
pdf(file="custom_outs/PV_HD_Unsupervised.pdf")
draw(hm)
dev.off()  

# Plot supervised clusters
hm <- Heatmap(
  t(selection),
  show_row_names = F,
  top_annotation = ra,
  heatmap_legend_param = list(title = "Z-Score",
                              labels_gp = gpar(fontsize = 8),
                              title_gp = gpar(fontsize = 10, fontface = "bold")),
  column_title = 'Supervised clustering of samples based on PV vs HD DEGs',
  column_split = colDat$CellType_Culture,
  cluster_columns = F,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8)
)

# Save
pdf(file="custom_outs/PV_HD_Supervised.pdf")
draw(hm)
dev.off() 

# Create heatmap for ET samples using custom gene list ----

# Define gene set of interest (goi)
goi <- c(
  'ADIPOQ',
  'PPARG',
  'CEBPA',
  'LPL',
  'ZEB2',
  'RUNX2',
  'SP7',
  'CTNNB1',
  'ALPL',
  'OCN',
  'COL2A1',
  'COL9A2',
  'SOX9',
  'ACAN'
)

# Extract vst counts matrix
vst <- as.data.frame(t(dat[[exp]][[model]]$dds@assays@data$vst))

# Extract gene symbols and IDs
genesLookup <-
  dat[[exp]][[model]]$dds@rowRanges@elementMetadata@listData[["symbol"]]

# Extract meta
colDat <- as.data.frame(dat[[exp]][[model]]$dds@colData)

# Add outlier sample annotation
outliers <- c('A12','A14','A16','A18')
colDat$Outlier[rownames(colDat) %in% outliers] <- 'True'
colDat$Outlier[!rownames(colDat) %in% outliers] <- 'False'

# Find the gene ID for a given symbol
symbol <- names(genesLookup)[genesLookup %in% goi]
genesLookupSub <- genesLookup[genesLookup %in% goi]

# Filter counts down to goi
selection <- as.matrix(vst[, c(symbol)])

# Subset down to select samples
soi <- rownames(colDat %>%
                  filter(
                    CellType_Culture %in% c(
                      'ET_MSC_ET_HSPC_Co-culture',
                      'ET_MSC_Mono-culture',
                      'HD_MSC_ET_HSPC_Co-culture',
                      'HD_MSC_Mono-culture'
                    )
                  ))
selection <- selection[soi, ]

# Change gene IDs to symbols
selection <- selection[, names(genesLookupSub)]
colnames(selection) <- genesLookupSub

# Scale
selection <- scale(selection, center = T, scale = T)

# Prep annotation
# Get rows of selection and colDat in same order
colDat <- colDat[rownames(selection),]

# Get annotation
ra = HeatmapAnnotation(
  CellType_Culture = colDat$CellType_Culture,
  Outlier = colDat$Outlier,
  col = list(
    CellType_Culture = c(
      'ET_MSC_ET_HSPC_Co-culture' = '#F0E442',
      'ET_MSC_Mono-culture' = '#E69F00',
      'HD_MSC_ET_HSPC_Co-culture' = '#0072B2',
      'HD_MSC_Mono-culture' = '#56B4E9'
    ),
    Outlier = c("True" = "black", "False" = "beige")
  ),
  gp = gpar(col = "black"),
  annotation_label = NULL,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 8),
                                 title_gp = gpar(fontsize = 10, fontface = "bold"))
)

# Plot unsupervised clusters
hm <- Heatmap(
  t(selection),
  show_row_names = T,
  top_annotation = ra,
  heatmap_legend_param = list(title = "Z-Score",
                              labels_gp = gpar(fontsize = 8),
                              title_gp = gpar(fontsize = 10, fontface = "bold")),
  column_title = 'Unsupervised clustering of samples based on custom list',
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)

# Save
pdf(file="custom_outs/ET_HD_Unsupervised_CustomList.pdf")
draw(hm)
dev.off()  

# Plot supervised clusters
hm <- Heatmap(
  t(selection),
  show_row_names = T,
  top_annotation = ra,
  heatmap_legend_param = list(title = "Z-Score",
                              labels_gp = gpar(fontsize = 8),
                              title_gp = gpar(fontsize = 10, fontface = "bold")),
  column_title = 'Supervised clustering of samples based on custom list',
  column_split = colDat$CellType_Culture,
  cluster_columns = F,
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)

# Save
pdf(file="custom_outs/ET_HD_Supervised_CustomList.pdf")
draw(hm)
dev.off() 



