library(GenomicFeatures)


args <- commandArgs(trailingOnly = TRUE)
gtf_path = args[1]
out_path = args[2]

### build bed file annotation for each gene ###
txdb <- makeTxDbFromGFF(gtf_path,format="gtf")
exons <- exonsBy(txdb,by="gene")
exons_temp = exonicParts(txdb)

temp = as.data.frame(exons_temp)
temp = temp[,c(1:5,8)]
temp <- lapply(1:nrow(temp),function(i){
  gene = unlist(temp[i,"gene_id"])
  gene = sapply(gene,function(y) unlist(strsplit(y,split = ".",fixed = T))[1])
  x = cbind(temp[i,1:5],gene)
  return(x)
})
temp <- as.data.frame(do.call(rbind,temp))
colnames(temp)[1] <- c("chr")

cache <- sapply(unique(temp$gene),function(x){
  sub_temp = temp[temp$gene == x,]
  filename = paste(c(out_path,x,".txt"),collapse = "")
  write.table(sub_temp[,1:5],file = filename,row.names = F,col.names = F,quote = F)
})

saveRDS(temp,file = paste(out_path,"gene_bed.rds",sep = "/"))

### build exon gtf ###
gtf_data = read.table(gtf_path,sep = "\t",comment.char = "#")
gtf_data = gtf_data[,c(3,4,5,9)]
colnames(gtf_data) = c("type","start","end","info")
gtf_data = gtf_data[gtf_data$type == "exon",]
gtf_data = gtf_data[,-1]

gtf_info_split <- function(info,keys = c(),split = "; ",sep = " "){
  delimiter = paste(split,sep,sep = "|")
  info = unlist(strsplit(info,split = delimiter))
  
  value = info[seq(2,length(info),2)]
  key = info[seq(1,length(info),2)]
  
  names(value) = key
  return(value)
}
info = lapply(gtf_data$info,function(x){
  out = gtf_info_split(x)
  return(out)
})
info = as.data.frame(do.call(rbind,info))

gtf_data = cbind(gtf_data,info)

gtf_data = gtf_data[,c("start","end","gene_id",
                       "transcript_type","transcript_name","exon_id")]
gtf_data$gene_id = sapply(gtf_data$gene_id,function(x){
  temp = unlist(strsplit(x,split = ".",fixed = TRUE))
  return(temp[1])
})

colnames(gtf_data)[c(4,5)] = c("transtype","transname")
saveRDS(gtf_data,file = paste(out_path,"exon_gtf.rds",sep = "/"))

