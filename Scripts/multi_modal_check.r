# %%
library(readr)
library(tidyverse)
library(SomaDataIO)
library(SomaScan.db)
library(optparse)

option_list <- list(
  make_option(c("-r", "--rds"), type = "character", default = NULL,
              help = "Path to the rds file", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Name of the output file", metavar = "character"),
  make_option(c("-d", "--dth"), type = "numeric", default = "0.85",
              help = "dth threshold", metavar = "numeric"),
  make_option(c("-t", "--tth"), type = "numeric", default = "0.1",
              help = "tth threshold", metavar = "numeric")
)

opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

rds_path <- opt$rds
data <- readRDS(rds_path)
ex_df <- data[[1]]
apt_df <- data[[2]]

output_file <- paste0(opt$output, ".qc.multimodal.rds")

# %%
phom_ex_df <- ex_df[grep("^PHOM", rownames(ex_df)), ]
naga_ex_df <- ex_df[grep("^NAG", rownames(ex_df)), ]

#---------------
uD <- naga_ex_df
log2_uD <- log2(uD)

#---------------
puD <- prcomp(phom_ex_df,scale=T)

k <- ncol(log2_uD)
cuD <- log2_uD

#-------------
lk <- {} 

pcuD <- prcomp(cuD,scale=T)


dth <- opt$dth
tth <- opt$tth

bimp <- {}

for(t in 1:k){
  pd <- density(cuD[,t],n=512)
  
  rdiff <- pd$y[-1]-pd$y[-512]
  ldiff <- pd$y[-512]-pd$y[-1]
  
  topd <- {}
  btmd <- {}
  for(i in 1:510){
    if(rdiff[i+1]<0 && ldiff[i]<0){
      topd <- c(topd,i)
    }
    if(rdiff[i+1]>0 && ldiff[i]>0){
      btmd <- c(btmd,i)
    } 
  }
  
  ti <- length(topd)
  bi <- length(btmd)
  
  hit <- FALSE
  maxy <- which.max(pd$y[topd])
  i <- 0
  while(i < (maxy-1)){
    i <- i+1   
    if(pd$y[topd[i]] > tth*pd$y[topd[maxy]]){
      for(j in 1:bi){     
        if(topd[i] < btmd[j] && btmd[j] < topd[maxy]){
          if(dth*pd$y[topd[i]] > pd$y[btmd[j]]){
            hit <- TRUE	   
            break
          }
        }
      }
    }
  }   
  
  if(!hit){
    i <- maxy+1
    while(i <= ti){
      if(pd$y[topd[i]] > tth*pd$y[topd[maxy]]){
        for(j in 1:bi){     
          if(topd[maxy] < btmd[j] && btmd[j] < topd[i]){
            if(dth*pd$y[topd[i]] > pd$y[btmd[j]]){
              hit <- TRUE	   
              break
            }
          }
        }
      }
      i <- i+1
    }
  }
  
  if(hit){
    bimp <- c(bimp,t)
  }
  
}

bimp_data <- cuD[, bimp]

soma_v4_apt <- somascan_menu$v4.0

apt_df <- apt_df %>%
  mutate(Is_SOMA_v4 = ifelse(SeqId %in% paste0("X", soma_v4_apt), TRUE, FALSE),
         multi_modal = ifelse(SeqId %in% colnames(bimp_data), TRUE, FALSE))

result <- list(ex_df, apt_df)
saveRDS(result, output_file)
