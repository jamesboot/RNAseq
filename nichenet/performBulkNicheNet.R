# Gene sets is output from geneSelection function - list of sender and reveiver genes
# Strategy is either Strat1 or Strat2
# CompName is the name of the overall comparison for saving directory
# gmt is pathway gmt file
# SenderName and ReceiverName are the cell type names for plotting
# SenderSamples and ReceiverSamples should be a vector of sample names for each group
# expressionMat is expression matrix
performBulkNicheNet <-
  function(GeneSets,
           Strategy,
           CompName,
           gmt,
           SenderName,
           ReceiverName,
           SenderSamples,
           ReceiverSamples,
           expressionMat,
           geneSetName) {
    strat <- Strategy
    expressed_genes_sender <- GeneSets$senderGenes
    expressed_genes_receiver <- GeneSets$receiverGenes
    
    # Create out directory if needed
    if (!dir.exists(file.path('NicheNet_Outs', strat, CompName, geneSetName))) {
      dir.create(file.path('NicheNet_Outs', strat, CompName, geneSetName), recursive = T)
    }
    
    # Save outdir
    outdir <- file.path('NicheNet_Outs', strat, CompName, geneSetName)
    
    # Plot heatmap of the expressed gene lists
    heatmap_mat <-
      expression[, c(expressed_genes_receiver, expressed_genes_sender)]
    rownames(heatmap_mat) == meta$ID
    ha <- HeatmapAnnotation(Cell_Type = meta$cell_type)
    p <- ComplexHeatmap::Heatmap(
      t(heatmap_mat),
      show_row_names = FALSE,
      use_raster = T,
      bottom_annotation = ha
    )
    pdf(
      file = paste0(outdir, '/expressed_genes_hmap.pdf'),
      height = 10,
      width = 10
    )
    plot(p)
    dev.off()
    
    # Filter the expressed ligands and receptors to only those that putatively bind together
    ligands <- lr_network %>% pull(from) %>% unique()
    expressed_ligands <- intersect(ligands, expressed_genes_sender)
    
    receptors <- lr_network %>% pull(to) %>% unique()
    expressed_receptors <-
      intersect(receptors, expressed_genes_receiver)
    
    potential_ligands <-
      lr_network %>% filter(from %in% expressed_ligands &
                              to %in% expressed_receptors) %>%
      pull(from) %>% unique()
    
    # Define gene set of interest
    # The gene set of interest consists of genes for which the expression is possibly affected due to communication with other cells
    if (is.character(gmt) & length(gmt) == 1) {
      genes <- read.delim(gmt, header = F)
      genevector <- c()
      for (x in 3:50) {
        genevector <- append(genevector, genes[, x])
      }
      genevector <- genevector[!duplicated(genevector)]
    } else if (is.vector(gmt)) {
      genevector <- gmt
    }
    
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
      intersect(best_upstream_ligands,
                colnames(active_ligand_target_links)) %>% rev()
    order_targets <-
      active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
    # Only want stuff actually in the expression matrix 
    order_targets <- order_targets[order_targets %in% colnames(expressionMat)]
    
    # Transpose the activeligand target links df
    vis_ligand_target <-
      t(active_ligand_target_links[order_targets, order_ligands])
    
    # Plot
    p_ligand_target_network <-
      make_heatmap_ggplot(
        vis_ligand_target,
        paste0("Prioritized ", SenderName, "-ligands"),
        paste0(geneSetName, " in ", ReceiverName, " cells"),
        color = "purple",
        legend_title = "Regulatory potential"
      ) +
      theme(legend.key.width= unit(2, 'cm')) +
      scale_fill_gradient2(low = "whitesmoke",  high = "purple")
    pdf(
      file = paste0(outdir, '/ligand_target_network1.pdf'),
      height = 5,
      width = 30
    )
    plot(p_ligand_target_network)
    dev.off()
    
    
    # Also plot receptors which can bind to top ligands
    ligand_receptor_links_df <-
      get_weighted_ligand_receptor_links(
        best_upstream_ligands,
        expressed_receptors,
        lr_network,
        weighted_networks$lr_sig
      )
    
    # Prepare visualisation
    vis_ligand_receptor_network <-
      prepare_ligand_receptor_visualization(ligand_receptor_links_df,
                                            best_upstream_ligands,
                                            order_hclust = "both")
    
    # Plot
    p_ligand_target_network2 <-
      make_heatmap_ggplot(
        t(vis_ligand_receptor_network),
        y_name = paste0("Prioritized ", SenderName, "-ligands"),
        x_name = paste0("Receptors expressed by ", ReceiverName, " cells"),
        color = "mediumvioletred",
        legend_title = "Prior interaction potential"
      ) +
      theme(legend.key.width= unit(1, 'cm'))
    pdf(
      file = paste0(outdir, '/ligand_target_network2.pdf'),
      height = 5,
      width = 5
    )
    plot(p_ligand_target_network2)
    dev.off()
    
    # More visualisations
    
    # Prepare ligand activity matrix
    vis_ligand_aupr <-
      ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
      column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
    
    # Plot ligand activity
    p_ligand_aupr <- make_heatmap_ggplot(
      vis_ligand_aupr,
      paste0("Prioritized ", SenderName, "-ligands"),
      "Ligand activity",
      color = "darkorange",
      legend_title = "AUPR"
    ) +
      theme(axis.text.x.top = element_blank()) +
      theme(legend.key.width= unit(1, 'cm'))
    pdf(
      file = paste0(outdir, '/ligand_aupr.pdf'),
      height = 5,
      width = 5
    )
    plot(p_ligand_aupr)
    dev.off()
    
    # Prepare expression of ligands in CAFs
    expression_df_sender <-
      expressionMat[c(SenderSamples), best_upstream_ligands]
    expression_df_sender <- t(expression_df_sender)
    
    # Set order
    vis_ligand_receiver_expression <-
      expression_df_sender[rev(best_upstream_ligands), c(SenderSamples)]
    
    # Set colour pallete
    color <-
      colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    
    # Plot
    p_ligand_receiver_expression <-
      make_heatmap_ggplot(
        vis_ligand_receiver_expression,
        paste0("Prioritized ", SenderName, "-ligands"),
        paste0(ReceiverName),
        color = color[100],
        legend_title = "Expression"
      ) +
      theme(legend.key.width= unit(1, 'cm'))
    pdf(
      file = paste0(outdir, '/ligand_receiver_expression.pdf'),
      height = 5,
      width = 5
    )
    plot(p_ligand_receiver_expression)
    dev.off()
    
    # Prepare expression of target genes
    expression_df_receiver <-
      expressionMat[c(ReceiverSamples), order_targets]
    
    # Scale
    vis_target_receiver_expression_scaled <-
      expression_df_receiver %>% scale_quantile()
    
    # Plot
    p_target_receiver_scaled_expression <-
      make_threecolor_heatmap_ggplot(
        vis_target_receiver_expression_scaled,
        paste0(ReceiverName),
        "Target",
        low_color = color[1],
        mid_color = color[50],
        mid = 0.5,
        high_color = color[100],
        legend_title = "Scaled expression"
      ) +
      theme(legend.key.width= unit(2, 'cm'))
    pdf(
      file = paste0(outdir, '/target_receiver_scaled_expression.pdf'),
      height = 5,
      width = 30
    )
    plot(p_target_receiver_scaled_expression)
    dev.off()
    
    # Combined plot
    figures_without_legend = plot_grid(
      p_ligand_aupr + theme(legend.position = "none"),
      p_ligand_receiver_expression + theme(legend.position = "none",
                                        axis.title.y = element_blank()),
      p_ligand_target_network + theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.title.y = element_blank()
      ),
      NULL,
      NULL,
      p_target_receiver_scaled_expression + theme(legend.position = "none",
                                                  axis.title.x = element_blank()),
      align = "hv",
      nrow = 2,
      rel_widths = c(
        ncol(vis_ligand_aupr) + 20,
        ncol(vis_ligand_receiver_expression) + 40,
        ncol(vis_ligand_target)
      ) - 2,
      rel_heights = c(
        nrow(vis_ligand_aupr),
        nrow(vis_target_receiver_expression_scaled) + 3
      )
    )
    
    legends = plot_grid(
      as_ggplot(get_legend(p_ligand_aupr)),
      as_ggplot(get_legend(p_ligand_receiver_expression)),
      as_ggplot(get_legend(p_ligand_target_network)),
      as_ggplot(get_legend(p_target_receiver_scaled_expression)),
      nrow = 2,
      align = "h"
    )
    
    p <- plot_grid(
      figures_without_legend,
      legends,
      rel_heights = c(10, 2),
      nrow = 2,
      align = "hv"
    )
    
    ggsave2(filename = paste0(outdir, '/panel_plot.pdf'),
            plot = p,
            height = 10,
            width = 40)
    
    # pdf(
    #   file = paste0(outdir, '/panel_plot.pdf'),
    #   height = 10,
    #   width = 40
    # )
    # plot_grid(
    #   figures_without_legend,
    #   legends,
    #   rel_heights = c(10, 2),
    #   nrow = 2,
    #   align = "hv"
    # )
    # dev.off()
    
  }
