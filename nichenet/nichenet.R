# Script for performing nichenet analysis on bulk RNAseq data
# Adapted from Pro's script
# 08/10/24

# Install packages if required
#devtools::install_github("saeyslab/nichenetr")
#devtools::install_github("jokergoo/ComplexHeatmap")

# Load packages
library(nichenetr)
library(tidyverse)
library(ComplexHeatmap)

# Prep - load everything according to nichenet vignette
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md
organism <- "mouse"

if (organism == "human") {
  lr_network <-
    readRDS(url(
      "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
    ))
  ligand_target_matrix <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
      )
    )
  weighted_networks <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"
      )
    )
} else if (organism == "mouse") {
  lr_network <-
    readRDS(url(
      "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"
    ))
  ligand_target_matrix <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
      )
    )
  weighted_networks <-
    readRDS(
      url(
        "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"
      )
    )
}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5, 1:5] # target genes in rows, ligands in columns

head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# Load expression matrix
expression <- read.csv('log_vst_symbols.csv', row.names = NULL)

# Load DGE analysis results
#dge <- readxl::read_xlsx('results/latest/differential_GRCm39_analyse_treatment.xlsx', sheet = 'C1.coculture - mono|MOC1|DMSO')

# Remove NA gene symbols
expression <- expression[!is.na(expression$symbol), ]

# Remove unwanted cols, group by symbol and summarise using mean
expression <- expression[, c(4:28)] %>%
  group_by(symbol) %>%
  summarise_all(mean) %>%
  column_to_rownames(var = 'symbol')

# Transpose
expression <- t(expression)

# Load experiment table (meta)
meta <- read.csv('extdata/metadata.csv')

# Modify ID column to match expression matrix
meta$ID <- gsub('GIA7103', '', meta$ID)

# Update gene symbols
# If this is not done, there will be 35 genes fewer in lr_network_expressed!
colnames(expression) <-
  convert_alias_to_symbols(colnames(expression), "mouse", verbose = T)

# START WITH DMSO COCULTURE CAFS > CANCER

# Get the CAF IDs and the Cancer IDs
cafID <- meta %>%
  filter(culture == 'coculture',
         treatment == 'DMSO',
         cell_type == 'MOCAF1_M1') %>%
  pull(ID)

cancerID <- meta %>%
  filter(culture == 'coculture',
         treatment == 'DMSO',
         cell_type == 'MOC1') %>%
  pull(ID)

# Determine which genes are expressed in senders (CAFs)
expressed_genes_sender <-
  names(sort(apply(expression[cafID, ], 2, mean), decreasing = T))[1:2000]

# Determine which genes are expressed in receivers (Cancers)
expressed_genes_receiver <-
  names(sort(apply(expression[cancerID, ], 2, mean), decreasing = T))[1:2000]

# Check overlap between two gene lists
length(intersect(expressed_genes_sender, expressed_genes_receiver))

# Plot heatmap of the expressed gene lists
heatmap_mat <- expression[,c(expressed_genes_receiver, expressed_genes_sender)]
rownames(heatmap_mat) == meta$ID
ha <- HeatmapAnnotation(Cell_Type = meta$cell_type)
ComplexHeatmap::Heatmap(t(heatmap_mat),
                        show_row_names = FALSE,
                        use_raster = T,
                        bottom_annotation = ha)

# Filter the expressed ligands and receptors to only those that putatively bind together
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands, expressed_genes_sender)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <-
  lr_network %>% filter(from %in% expressed_ligands &
                          to %in% expressed_receptors) %>%
  pull(from) %>% unique()

# Define gene set of interest
# The gene set of interest consists of genes for which the expression is possibly affected due to communication with other cells
genes <- read.delim('mh.all.v2024.1.Mm.symbols.gmt', header = F)
genevector <- c()
for (x in 3:50) {
  genevector <- append(genevector, genes[, x])
}
genevector <- genevector[!duplicated(genevector)]

# Set background
background_expressed_genes <-
  expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)

# Perform NicheNet analysis
ligand_activities <- predict_ligand_activities(
  geneset = genevector,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Sort ligand activities
ligand_activities <-
  ligand_activities %>% arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

# See top ligands
best_upstream_ligands <-
  ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands

# Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = genevector,
    ligand_target_matrix = ligand_target_matrix,
    n = 200
  ) %>% bind_rows()

# Prepare ligand target visualisation
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25
)

# Set order of ligands and targets
order_ligands <-
  intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <-
  active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

# Transpose the activeligand target links df
vis_ligand_target <-
  t(active_ligand_target_links[order_targets, order_ligands])

# Plot
p_ligand_target_network <-
  make_heatmap_ggplot(
    vis_ligand_target,
    "Prioritized CAF-ligands",
    "Hallmark genes in malignant cells",
    color = "purple",
    legend_title = "Regulatory potential"
  ) +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
p_ligand_target_network

# Also plot receptors which can bind to top ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands,
                                                               expressed_receptors,
                                                               lr_network,
                                                               weighted_networks$lr_sig)

# Prepare visualisation
vis_ligand_receptor_network <-
  prepare_ligand_receptor_visualization(ligand_receptor_links_df,
                                        best_upstream_ligands,
                                        order_hclust = "both")

# Plot
p_ligand_target_network2 <-
  make_heatmap_ggplot(
    t(vis_ligand_receptor_network),
    y_name = "Prioritized CAF-ligands",
    x_name = "Receptors expressed by malignant cells",
    color = "mediumvioletred",
    legend_title = "Prior interaction potential"
  )
p_ligand_target_network2

# More visualisations

# Required packages
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Prepare ligand activity matrix
vis_ligand_aupr <-
  ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

# Plot ligand activity
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritized CAF-ligands",
  "Ligand activity",
  color = "darkorange",
  legend_title = "AUPR"
) +
  theme(axis.text.x.top = element_blank())
p_ligand_aupr

# Prepare expression of ligands in CAFs
expression_df_CAF <- expression[c(cafID), best_upstream_ligands]
expression_df_CAF <- t(expression_df_CAF)

# Set order
vis_ligand_tumor_expression <-
  expression_df_CAF[rev(best_upstream_ligands), c(cafID)]

# Set colour pallete
color <-
  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# Plot
p_ligand_tumor_expression <-
  make_heatmap_ggplot(
    vis_ligand_tumor_expression,
    "Prioritized CAF-ligands",
    "Tumor",
    color = color[100],
    legend_title = "Expression"
  )
p_ligand_tumor_expression


# Prepare expression of target genes
expression_df_cancer <-
  expression[c(cancerID), colnames(expression) %in% order_targets]

# Scale
vis_target_tumor_expression_scaled <-
  expression_df_cancer %>% scale_quantile()

# Plot
p_target_tumor_scaled_expression <-
  make_threecolor_heatmap_ggplot(
    vis_target_tumor_expression_scaled,
    "Tumor",
    "Target",
    low_color = color[1],
    mid_color = color[50],
    mid = 0.5,
    high_color = color[100],
    legend_title = "Scaled expression"
  )
p_target_tumor_scaled_expression

# Combined plot
figures_without_legend = plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_ligand_tumor_expression + theme(legend.position = "none",
                                    axis.title.y = element_blank()),
  p_ligand_target_network + theme(legend.position = "none",
                                  axis.ticks = element_blank(),
                                  axis.title.y = element_blank()), 
  NULL,
  NULL,
  p_target_tumor_scaled_expression + theme(legend.position = "none",
                                           axis.title.x = element_blank()), 
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_aupr)+20, ncol(vis_ligand_tumor_expression)+40, ncol(vis_ligand_target))-2,
  rel_heights = c(nrow(vis_ligand_aupr), nrow(vis_target_tumor_expression_scaled)+3)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_aupr)),
  as_ggplot(get_legend(p_ligand_tumor_expression)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_target_tumor_scaled_expression)),
  nrow = 2,
  align = "h")

pdf(file = 'coculture_dmso_CAFs_to_malignant.pdf',
    height = 10,
    width = 40)
plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,2), nrow = 2, align = "hv")
dev.off()

# START WITH DMSO MONOCULTURE CAFS > CANCER

# Get the CAF IDs and the Cancer IDs
cafID<- meta %>%
  filter(culture == 'mono',
         treatment == 'DMSO',
         cell_type == 'MOCAF1_M1') %>%
  pull(ID)

cancerID <- meta %>%
  filter(culture == 'mono',
         treatment == 'DMSO',
         cell_type == 'MOC1') %>%
  pull(ID)

# Determine which genes are expressed in senders (CAFs)
expressed_genes_sender <-
  names(sort(apply(expression[cafID, ], 2, mean), decreasing = T))[1:2000]

# Determine which genes are expressed in receivers (Cancers)
expressed_genes_receiver <-
  names(sort(apply(expression[cancerID, ], 2, mean), decreasing = T))[1:2000]

# Check overlap between two gene lists
length(intersect(expressed_genes_sender, expressed_genes_receiver))

# Plot heatmap of the expressed gene lists
heatmap_mat <- expression[,c(expressed_genes_receiver, expressed_genes_sender)]
rownames(heatmap_mat) == meta$ID
ha <- HeatmapAnnotation(Cell_Type = meta$cell_type)
ComplexHeatmap::Heatmap(t(heatmap_mat),
                        show_row_names = FALSE,
                        use_raster = T,
                        bottom_annotation = ha)

# Filter the expressed ligands and receptors to only those that putatively bind together
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands, expressed_genes_sender)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <-
  lr_network %>% filter(from %in% expressed_ligands &
                          to %in% expressed_receptors) %>%
  pull(from) %>% unique()

# Define gene set of interest
# The gene set of interest consists of genes for which the expression is possibly affected due to communication with other cells
genes <- read.delim('mh.all.v2024.1.Mm.symbols.gmt', header = F)
genevector <- c()
for (x in 3:50) {
  genevector <- append(genevector, genes[, x])
}
genevector <- genevector[!duplicated(genevector)]

# Set background
background_expressed_genes <-
  expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
length(background_expressed_genes)

# Perform NicheNet analysis
ligand_activities <- predict_ligand_activities(
  geneset = genevector,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Sort ligand activities
ligand_activities <-
  ligand_activities %>% arrange(-aupr_corrected) %>%
  mutate(rank = rank(desc(aupr_corrected)))

# See top ligands
best_upstream_ligands <-
  ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)
best_upstream_ligands

# Infer target genes and receptors of top-ranked ligands
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = genevector,
    ligand_target_matrix = ligand_target_matrix,
    n = 200
  ) %>% bind_rows()

# Prepare ligand target visualisation
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25
)

# Set order of ligands and targets
order_ligands <-
  intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <-
  active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

# Transpose the activeligand target links df
vis_ligand_target <-
  t(active_ligand_target_links[order_targets, order_ligands])

# Plot
p_ligand_target_network <-
  make_heatmap_ggplot(
    vis_ligand_target,
    "Prioritized CAF-ligands",
    "Hallmark genes in malignant cells",
    color = "purple",
    legend_title = "Regulatory potential"
  ) +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
p_ligand_target_network

# Also plot receptors which can bind to top ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands,
                                                               expressed_receptors,
                                                               lr_network,
                                                               weighted_networks$lr_sig)

# Prepare visualisation
vis_ligand_receptor_network <-
  prepare_ligand_receptor_visualization(ligand_receptor_links_df,
                                        best_upstream_ligands,
                                        order_hclust = "both")

# Plot
p_ligand_target_network2 <-
  make_heatmap_ggplot(
    t(vis_ligand_receptor_network),
    y_name = "Prioritized CAF-ligands",
    x_name = "Receptors expressed by malignant cells",
    color = "mediumvioletred",
    legend_title = "Prior interaction potential"
  )
p_ligand_target_network2

# More visualisations

# Required packages
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Prepare ligand activity matrix
vis_ligand_aupr <-
  ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

# Plot ligand activity
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritized CAF-ligands",
  "Ligand activity",
  color = "darkorange",
  legend_title = "AUPR"
) +
  theme(axis.text.x.top = element_blank())
p_ligand_aupr

# Prepare expression of ligands in CAFs
expression_df_CAF <- expression[c(cafID), best_upstream_ligands]
expression_df_CAF <- t(expression_df_CAF)

# Set order
vis_ligand_tumor_expression <-
  expression_df_CAF[rev(best_upstream_ligands), c(cafID)]

# Set colour pallete
color <-
  colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

# Plot
p_ligand_tumor_expression <-
  make_heatmap_ggplot(
    vis_ligand_tumor_expression,
    "Prioritized CAF-ligands",
    "Tumor",
    color = color[100],
    legend_title = "Expression"
  )
p_ligand_tumor_expression


# Prepare expression of target genes
expression_df_cancer <-
  expression[c(cancerID), colnames(expression) %in% order_targets]

# Scale
vis_target_tumor_expression_scaled <-
  expression_df_cancer %>% scale_quantile()

# Plot
p_target_tumor_scaled_expression <-
  make_threecolor_heatmap_ggplot(
    vis_target_tumor_expression_scaled,
    "Tumor",
    "Target",
    low_color = color[1],
    mid_color = color[50],
    mid = 0.5,
    high_color = color[100],
    legend_title = "Scaled expression"
  )
p_target_tumor_scaled_expression

# Combined plot
figures_without_legend = plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_ligand_tumor_expression + theme(legend.position = "none",
                                    axis.title.y = element_blank()),
  p_ligand_target_network + theme(legend.position = "none",
                                  axis.ticks = element_blank(),
                                  axis.title.y = element_blank()), 
  NULL,
  NULL,
  p_target_tumor_scaled_expression + theme(legend.position = "none",
                                           axis.title.x = element_blank()), 
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_aupr)+20, ncol(vis_ligand_tumor_expression)+40, ncol(vis_ligand_target))-2,
  rel_heights = c(nrow(vis_ligand_aupr), nrow(vis_target_tumor_expression_scaled)+3)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_aupr)),
  as_ggplot(get_legend(p_ligand_tumor_expression)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_target_tumor_scaled_expression)),
  nrow = 2,
  align = "h")

pdf(file = 'mono_dmso_CAFs_to_malignant.pdf',
    height = 10,
    width = 40)
plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,2), nrow = 2, align = "hv")
dev.off()

 