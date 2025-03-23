# Custom MA plots

# Load packages
library(DESeq2)
library(ggplot2)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(ggpp)

# Load data
ddsList <- readRDS('babs_v1/differential/data/analyse_x_GRCm38.rds')

# Create directory for outputs to be saved to if doesn't exist
outdir <- 'R_custom_outs2'
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# Nested for loop to run through all comparisons
for (experiment in names(ddsList)) {
  models <- names(ddsList[[experiment]])
  for (model in models) {
    comps <- names(ddsList[[experiment]][[model]]$comps)
    for (comparison in comps) {
      # Pull out results
      res <- as.data.frame(ddsList[[experiment]][[model]]$comps[[comparison]])
      # Define levels of class
      res$class <- factor(res$class,
                          levels = c("Up*",
                                     "Up",
                                     "Down",
                                     "Down*", 
                                     "Low Count",
                                     "Zero Count",
                                     "Outlier"))
      # Select labels to plot
      # Get top Up* genes
      topUp <- res %>%
        filter(class == 'Up*') %>%
        arrange(desc(shrunkLFC)) %>%
        slice_head(n = 10) %>%
        na.omit() %>%
        pull(symbol)
      # Get top Down* genes
      topDown <- res %>%
        filter(class == 'Down*') %>%
        arrange(shrunkLFC) %>%
        slice_head(n = 10) %>%
        na.omit() %>%
        pull(symbol)
      # Merge
      labs <- c(topUp, topDown)
      # Plot
      pl <- ggplot(res, aes(
        x = baseMean,
        y = shrunkLFC,
        colour = class
      )) +
        geom_point(size = 0.5, alpha = 0.5) +
        scale_color_manual(values = c(
          "Down" = 'lightblue',
          "Up" = 'tomato',
          "Low Count" = 'grey',
          "Up*" = 'darkred',
          "Zero Count" = 'black',
          "Down*" = 'royalblue',
          "Outlier" = 'green3'
        )) +
        geom_text_repel(
          data = subset(res, symbol %in% labs),
          aes(label = symbol),
          box.padding = 1,
          label.padding = 1,
          point.padding = 0.5,
          force = 2,
          max.overlaps = Inf,
          colour = 'black',
          position = position_nudge_center(y = 2, center_y = 0)
        ) +
        guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), title = 'Class')) +
        xlab('Mean of normalised counts') +
        ylab('Log fold change (shrunk)') +
        scale_x_log10() +
        theme_bw()
      # Save
      ggsave(plot = pl,
             filename = paste0(outdir,
                               '/',
                               gsub('\\|', '-', comparison),
                               '.pdf'))
    }
  }
}


