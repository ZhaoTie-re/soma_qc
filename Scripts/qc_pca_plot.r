# %%
#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(plotly)
library(htmlwidgets)
library(optparse)

option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to the rds file", metavar = "character"),
  make_option(c("-s", "--summary"), type = "character", default = NULL,
              help = "Path to the adat summary file", metavar = "character"), 
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character")
)

opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

output_file <- paste0(opt$output, ".qc.pca_plot.html")

# %%
rds <- readRDS(opt$rds)
ex_df <- rds[[1]]

adat_summary <- read_excel(opt$summary)
replacement <- setNames(adat_summary$adat_info, adat_summary$adat_prefix)

rownames(ex_df) <- str_replace_all(rownames(ex_df), replacement)

# %%
rds_log <- log(ex_df)
p_rds_log <- prcomp(rds_log, scale = TRUE)
pca_results <- data.frame(PC1 = p_rds_log$x[,1],
                          PC2 = p_rds_log$x[,2],
                          PC3 = p_rds_log$x[,3])

# %%
sample_groups <- substr(rownames(pca_results), 1, 3)

colors <- colorRampPalette(c("red", "blue", "green"))(length(unique(sample_groups)))
color_mapping <- setNames(colors, unique(sample_groups))
sample_colors <- color_mapping[sample_groups]

explained_variance <- summary(p_rds_log)$importance[2, ]
total_explained_variance <- sum(explained_variance[1:3])

x_label <- paste("PC1 (", round(explained_variance[1] * 100, 2), "%)", sep = "")
y_label <- paste("PC2 (", round(explained_variance[2] * 100, 2), "%)", sep = "")
z_label <- paste("PC3 (", round(explained_variance[3] * 100, 2), "%)", sep = "")

plot_title <- paste("3D PCA Plot (Total variance explained: ", round(total_explained_variance * 100, 2), "%)", sep = "")

p <- plot_ly(pca_results, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers",
             marker = list(size = 3, color = sample_colors, opacity = 0.6, symbol = 'circle', line = list(color = 'rgba(204, 204, 204, 0.5)', width = 0.5)),
             hovertext = rownames(pca_results), hoverinfo = 'text') %>%
  layout(scene = list(xaxis = list(title = x_label),
                      yaxis = list(title = y_label),
                      zaxis = list(title = z_label)),
         title = plot_title)

# print(p)

saveWidget(p, output_file)
