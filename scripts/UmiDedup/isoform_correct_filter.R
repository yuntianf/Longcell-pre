# input:
# 1. reads cluster for each cell for a gene
# 2. output file name

source("./splice_site_correct.R")
source("./reads_filter.R")
source("./umi_count_functions.R")
args <- commandArgs(trailingOnly = TRUE)
gene_cells_cluster = readRDS(args[1])

genes_umi_count = isoform_correct_filter(gene_cells_cluster,
                                         edit = gene_cells_cluster$edit)
genes_umi_count$gene = unique(gene_cells_cluster$gene_ID)
write.table(genes_umi_count,file = args[2], 
            row.names = FALSE,col.names = TRUE,
            quote = FALSE, sep = "\t")
