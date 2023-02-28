### args:
# 1. isoform table
# 2. gene bed annotation
# 3. output name
options(warn = -1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(data.table)
})

source("./spliceob/exon_corres.R")

args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], header = TRUE)
gene_bed <- readRDS(args[2])

exonTable= createExonList(data,gene_bed)

write.table(exonTable,file = args[3], row.names = FALSE,col.names = TRUE,quote = FALSE, sep = "\t")
