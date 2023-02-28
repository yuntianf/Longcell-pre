### args:
# 1. isoform table
# 3. output name

options(warn = -1)
library(tidyverse)
library(dplyr)
library(data.table)

source("./splice_object.R")

args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1], header = TRUE)

ob = exonList2spliceOb(data)

saveRDS(ob,file = args[2])
