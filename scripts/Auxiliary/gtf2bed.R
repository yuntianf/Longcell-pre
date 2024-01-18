suppressPackageStartupMessages(library(GenomicFeatures))
library(parallel)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("-g", "--gtf", type="character", 
                    help="The file path for the input gtf")
parser$add_argument("-o", "--out", type="character", 
                    help="The output directory")
parser$add_argument("-c", "--cores", type="integer", default=4,
                    help="The number of cores for parallization")

args <- parser$parse_args()

gtf_path = args$gtf
out_path = args$out
cores = args$cores

dir.create(out_path, showWarnings = FALSE)

### build bed file annotation for each gene ###
txdb <- makeTxDbFromGFF(gtf_path,format="gtf")
exons <- exonsBy(txdb,by="gene")
exons_temp = exonicParts(txdb)

temp = as.data.frame(exons_temp)
temp = temp[,c(1:5,8)]
temp[,"gene_id"] = sapply(temp[,"gene_id"], `[[`, 1)
temp[,"gene_id"] = gsub("\\..*$", "", temp[,"gene_id"])
colnames(temp)[1] <- c("chr")

core_avai = parallel::detectCores()
cores = ifelse(cores > core_avai,core_avai,cores)

cache <- mclapply(unique(temp$gene),function(x){
  sub_temp = temp[temp$gene == x,]
  filename = paste(c(out_path,x,".txt"),collapse = "")
  write.table(sub_temp[,1:5],file = filename,row.names = F,col.names = F,quote = F)
},mc.cores = getOption("mc.cores", cores))

saveRDS(temp,file = paste(out_path,"gene_bed.rds",sep = "/"))

### build exon gtf ###
gtf_data = as.data.frame(exons(txdb, columns=c("gene_id","tx_name","exon_name"), 
                               filter=NULL, use.names=FALSE))
gtf_data[,"gene_id"] = sapply(gtf_data[,"gene_id"], `[[`, 1)
gtf_data[,"tx_name"] = sapply(gtf_data[,"tx_name"], `[[`, 1)

gtf_data$gene_id = substr(gtf_data$gene_id,1,15)
gtf_data = gtf_data[,c(2:3,6:8)]
colnames(gtf_data)[4:5] = c("transname","exon_id")

saveRDS(gtf_data,file = paste(out_path,"exon_gtf.rds",sep = "/"))

