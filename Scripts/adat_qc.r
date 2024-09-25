# %%
#!/usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("-d", "--data"), type = "character", default = NULL,
              help = "Path to ADAT file", metavar = "character"),
  make_option(c("-s", "--sample"), type = "character", default = NULL,
              help = "Path to the sample file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character")
)

opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

adat_data <- opt$data
sample_data <- opt$sample
output_file <- paste0(opt$output, ".qc.rds")

# %%
library(SomaDataIO)
library(tidyverse)

adat_qc <- function(adat_path, sample_path) {


  adat_data <- read_adat(adat_path)
  sample_df <- read.csv(sample_path, sep = "\t", header = TRUE)

  aptamer_qc <- function(my_adat) {
    aptamer_info <- getAnalyteInfo(my_adat)
    aptamer_info <- aptamer_info %>%
      filter(ColCheck == "PASS",
             UniProt != "",
             !is.na(UniProt),
             Organism == "Human",
             Type == "Protein")
    return(aptamer_info)
  }

  sample_qc <- function(my_adat) {
    r_m <- getMeta(my_adat)
    sample_info <- tibble(my_adat)[, r_m]
    sample_info <- sample_info %>%
      filter(RowCheck == "PASS",
             SampleType %in% c("Sample", "QC"))
    return(sample_info)
  }


  aptamer_after_qc <- aptamer_qc(adat_data)
  sample_after_qc <- sample_qc(adat_data)

  ls <- sample_df %>%
    select(ExtIdentifier, SampleCode)


  ex <- tibble(adat_data) %>%
    filter(SampleType %in% c("QC", "Sample") & RowCheck == "PASS") %>%
    select(which(!grepl("^seq.", names(.)) | names(.) %in% aptamer_after_qc$AptName)) %>%
    rename_with(~ str_replace_all(., c("^seq." = "X", "\\." = "-"))) %>%
    left_join(ls, by = "ExtIdentifier") %>%
    select(ExtIdentifier, SampleCode, everything()) %>%
    mutate(SampleCode = ifelse(SampleType == "QC", paste0("QC_", ExtIdentifier, "_in_", opt$output), SampleCode))



  ex_df <- ex %>%
    select(SampleCode, starts_with("X")) %>%
    column_to_rownames("SampleCode")


  apt_df <- aptamer_after_qc %>%
    mutate(SeqId = paste0("X", SeqId))

  return(list(ex_df, apt_df))
}

# %%
result <- adat_qc(adat_data, sample_data)
saveRDS(result, output_file)
