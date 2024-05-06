#!/usr/bin/env Rscript
library(readr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

#This script renames all NGS samples according to the meta info sample table. Mapping via sample-id
#If sample-id is not listed in the sample table, the file gets deleted.

meta_file <- read_tsv(file=args[2], col_names = T, skip_empty_rows = T)
sample_name_dict <- paste0(meta_file$sample_id, "_", meta_file$condition, ".fastq.gz")
names(sample_name_dict) <- meta_file$sample_id

setwd(args[1])

samples <- list.files(pattern="fastq.gz")

for(s in samples){
  sample_id <- s %>% str_split(pattern="_") %>% unlist 
  sample_id<-sample_id[1]
  print(sample_id)
  
  sample_name <- sample_name_dict[sample_id]
  
  if(!is.na(sample_name)){
    file.rename(s,sample_name)
  }else{
    file.remove(s)
  }
}
