# args:
# 1. cell exon softclips
# 2. number of files out
# 3. cell exon output folder

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

cell_exon <- read.table(args[1],fill = T)

cell_exon <- as.data.frame(cell_exon)

file_num = as.numeric(args[2])
genes = unique(cell_exon$gene)

gene_num = length(genes) %/% file_num
cache <- sapply(1:(file_num+1),function(i){
  start = (i-1)*gene_num+1
  if(start > length(genes)){
    return(NULL)
  }
  end = ifelse(i*gene_num >= length(genes),length(genes),i*gene_num)
  sub_genes <- genes[start:end]
  sub_cell_exon = cell_exon[cell_exon$gene %in% sub_genes,]
  file_name = paste(c("sub_cell_exon",i,"txt"),collapse = ".")
  file_name = paste(args[5],file_name,sep = "/")
  
  write.table(sub_cell_exon,file = file_name,row.names = F,col.names = T,quote = F)
})

