# %%
#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("-r", "--rds_1"), type = "character", default = NULL,
              help = "Path to the first rds file", metavar = "character"),
  make_option(c("-s", "--rds_2"), type = "character", default = NULL,
              help = "Path to the second rds file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character")
)

opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
output_file <- paste0(opt$output, ".qc.merge.rds")

# %%
rds_1 <- readRDS(opt$rds_1)
rds_2 <- readRDS(opt$rds_2)

# %%
apt_info <- c("AptName", "SeqId", "TargetFullName", "Target", "UniProt", "EntrezGeneID", "EntrezGeneSymbol", "Organism", "Type")

ex_df_1 <- rds_1[[1]]
apt_df_1 <- rds_1[[2]] %>% select(all_of(apt_info))

ex_df_2 <- rds_2[[1]]
apt_df_2 <- rds_2[[2]] %>% select(all_of(apt_info))

# %%
inter_col <- intersect(colnames(ex_df_1), colnames(ex_df_2))

ex_df_1_selected <- ex_df_1 %>% select(all_of(inter_col))
ex_df_2_selected <- ex_df_2 %>% select(all_of(inter_col))

inter_ex_df <- bind_rows(ex_df_1_selected, ex_df_2_selected)

# %%
apt_df_1_selected <- apt_df_1 %>% filter(SeqId %in% inter_col)
apt_df_2_selected <- apt_df_2 %>% filter(SeqId %in% inter_col)

# %%
if (identical(apt_df_1_selected, apt_df_2_selected)) {
  inter_apt_df <- apt_df_1_selected
} else {
  warning("apt_df_1_selected and apt_df_2_selected are not identical.")
  inter_apt_df <- apt_df_1_selected
}

# %%
result <- list(inter_ex_df, inter_apt_df)
saveRDS(result, output_file)
