# 1.exon reads table
# 2.gene bed annotation rds
# 3.output folder
args <- commandArgs(trailingOnly = TRUE)
exon_reads <- read.table(args[1])
exon_reads <- as.data.frame(exon_reads[,c(2,8,9)])
colnames(exon_reads) <- c("gene","start","end")

gene_bed <- readRDS(args[2])
gene_bed <- gene_bed[gene_bed$gene == unique(exon_reads$gene),]

map_error <- sum(exon_reads$start < gene_bed[1,"start"] | 
                    exon_reads$end > gene_bed[nrow(gene_bed),"end"])
map_error <- map_error/nrow(exon_reads)
map_error <- matrix(c(unique(exon_reads$gene),map_error,nrow(exon_reads)),nrow = 1)
write.table(map_error,file = paste(args[3],unique(exon_reads$gene),sep = "/"),row.names = F, col.names = F,quote = F)