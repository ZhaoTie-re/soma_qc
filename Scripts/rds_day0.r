# %%
#!/usr/bin/env Rscript
library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to the rds file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character")
)

opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
output_file <- paste0(opt$output, ".qc.day0.rds")

# ex_df_csv <- paste0(opt$output, ".ex.day0.csv")
# apt_df_csv <- paste0(opt$output, ".apt.day0.csv")

# %%
rds <- readRDS(opt$rds)

apt_info <- c("AptName", "SeqId", "TargetFullName", "Target", "UniProt", "EntrezGeneID", "EntrezGeneSymbol", "Organism", "Type")

ex_df <- rds[[1]]
apt_df <- rds[[2]] %>% select(all_of(apt_info))

ex_df_day0 <- ex_df %>% filter(grepl("day0|QC", rownames(.)))

result <- list(ex_df_day0, apt_df)
saveRDS(result, output_file)

# write.csv(ex_df_day0, ex_df_csv, row.names = TRUE, col.names = TRUE)
# write.csv(apt_df, apt_df_csv, row.names = FALSE, col.names = TRUE)