### args:
# 1. exon table
# 2. path

options(warn = -1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(data.table)
})

source("./spliceob/exon_corres.R")

args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], header = TRUE)

cache = saveExonList(data,path = args[2],isoform_col = "exons")
