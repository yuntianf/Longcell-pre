### args
# 1. input cell exons
# 2. needleman score thresh to build graph
# 3. output file
suppressPackageStartupMessages({
  library(stats4)
  library(Rcpp)
  library(dbscan)
  library(e1071)
  library(dplyr)
  library(reshape2)
  library(tidyr)
  library(transport)
  library(RcppHungarian)
  library(NameNeedle)
  library(igraph)
  library(MASS)
  library(argparse)
})

sourceCpp("./UmiDedup/umi_dist.cpp")
source("./UmiDedup/umi_cluster.R")
source("./UmiDedup/splice_site_correct.R")
source("./UmiDedup/reads_filter.R")
source("./UmiDedup/umi_count_functions.R")

argparser <- function(){
  parser <- ArgumentParser("UMI deduplication")
  
  parser$add_argument("-c", "--cell_exon", 
                      help="The input file to do UMI deduplication")
  parser$add_argument("-s", "--thresh", type="integer", default=5,
                      help="threshold for the similarity between UMI")
  parser$add_argument("-b", "--barlen", type="integer", default=16, 
                      help="the length of the barcode")
  parser$add_argument("-u", "--umilen", type="integer", default=10, 
                      help="the length of the UMI")
  parser$add_argument("-f", "--flank", type="integer", default=1, 
                      help="the length of flank around the UMI to be tolerant of insertions and deletions")
  parser$add_argument("-ss","--splice_site_thresh", type="integer", default=10, 
                      help = "threshold to filter out infrequent splice sites")
  parser$add_argument("-o", "--outfile", 
                      help="The output file")
  
  return(parser)
}

args <- argparser()$parse_args()

cell_exon <- read.table(args$cell_exon, header = TRUE)


colnames(cell_exon) <- c("reads","barcode","start","edit","search_seq",
                                           "gene_ID","exon_seq","state","polyA","strand")

cat("There are ",length(unique(cell_exon$barcode)),
    " cells with ",length(unique(cell_exon$gene_ID))," genes\n")

sim_thresh = args$thresh

out = umi_count(cell_exon,sim_thresh = sim_thresh,bar_len = args$barlen,
                flank = args$flank,umi_len = args$umilen,
                splice_site_thresh = args$splice_site_thresh)

write.table(out,file = args$outfile, row.names = FALSE,col.names = TRUE,quote = FALSE, sep = "\t")
