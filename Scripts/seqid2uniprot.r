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

output_file <- paste0(opt$output, ".uniprot.rds")

# %%
rds <- readRDS(opt$rds)
ex_df <- rds[[1]]
apt_df <- rds[[2]]

# Remove rows where the row name contains "QC"
ex_df <- ex_df[!grepl("QC", rownames(ex_df)), ]

# %%
# Recognize one SeqId corresponds to multiple UniProt IDs

seqid_1_multi_list <- list() # List to store SeqId corresponds to multiple UniProt IDs # nolint

for (i in seq_len(nrow(apt_df))) {
  # If the UniProt value contains "|", add the SeqId to the list
  if (str_detect(apt_df$UniProt[i], "\\|")) {
    seqid_1_multi_list[[length(seqid_1_multi_list) + 1]] <- apt_df$SeqId[i]
  }
}

# %%
# Recognize multiple SeqId corresponds to one UniProt ID

multi_uniprot_df <- apt_df %>%
  filter(!str_detect(UniProt, "\\|")) %>%
  group_by(UniProt) %>%
  summarise(SeqId = paste(SeqId, collapse = ",")) %>%
  filter(str_detect(SeqId, ","))

max_seqid_list <- list() # List to store the SeqId with the highest average expression value (multipy SeqId to one UniProtID) # nolint
other_seqid_list <- list() # List to store the other SeqId (multipy SeqId to one UniProtID) # nolint

for (i in seq_len(nrow(multi_uniprot_df))) {

  seqids <- str_split(multi_uniprot_df$SeqId[i], ",")[[1]]

  avg_values <- sapply(seqids, function(seqid) mean(ex_df[, seqid], na.rm = TRUE))
  
  max_seqid <- seqids[which.max(avg_values)]
  max_seqid_list <- append(max_seqid_list, max_seqid)
  
  other_seqids <- seqids[seqids != max_seqid]
  other_seqid_list <- append(other_seqid_list, other_seqids)
}

# %%
uniProt_groups <- apt_df %>% group_by(UniProt) %>% summarise(n = n())
uniProt_single_seqId <- uniProt_groups %>% filter(n == 1, !str_detect(UniProt, "\\|"))
seqid_1_1_list <- apt_df %>% filter(UniProt %in% uniProt_single_seqId$UniProt) %>% pull(SeqId)

# %%
seqid_1_multi_list <- unlist(seqid_1_multi_list) # 74 in cteph_day0
max_seqid_list <- unlist(max_seqid_list) # 761 in cteph_day0
other_seqid_list <- unlist(other_seqid_list) # 861  in cteph_day0
seqid_1_1_list <- unlist(seqid_1_1_list) # 5493 in cteph_day0

# %%
seqid_res_list <- union(seqid_1_1_list, max_seqid_list)
seqid_del_list <- union(seqid_1_multi_list, other_seqid_list)

# %%
# SOMAmers Matrix to UniProt Matrix
apt_df_uniprot <- apt_df %>% 
  filter(SeqId %in% seqid_res_list)

ex_df_uniprot <- ex_df %>%
  select(all_of(seqid_res_list)) %>%
  setNames(apt_df_uniprot$UniProt[match(seqid_res_list, apt_df_uniprot$SeqId)])

# %%
# Write SOMAmers prcossing summary
apt_summary <- apt_df %>% mutate(category = case_when(
  SeqId %in% seqid_1_1_list ~ "1 SOMAmer to 1 UniProt",
  SeqId %in% other_seqid_list ~ "Multiple SOMAmers to 1 UniProt (Not Max Mean)",   
  SeqId %in% max_seqid_list ~ "Multiple SOMAmers to 1 UniProt (Max mean)", 
  SeqId %in% seqid_1_multi_list ~ "1 SOMAmer to multiple UniProt",
  TRUE ~ "Other" # Ideal state doesn't exist.
))

apt_summary <- apt_summary %>%mutate(status = case_when(
  SeqId %in% seqid_del_list ~ "Deleted",
  SeqId %in% seqid_res_list ~ "Reserved",
  TRUE ~ "Other" # Ideal state doesn't exist.
))

# %%
result <- list(ex_df_uniprot, apt_df_uniprot, apt_summary)
saveRDS(result, output_file)
