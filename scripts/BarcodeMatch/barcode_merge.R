 # args:
# 1. barcode match
# 2. cell exon softclips
# 3. number of files out
# 4. cell exon output folder

library(dplyr)
library(knapsack)

BARCODE_LEN = 16

args <- commandArgs(trailingOnly = TRUE)


#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
#                      functions                       #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
### data distribution function for parallization
distribute <- function(data,num){
    data <- as.data.frame(data)
    data[,2] <- as.numeric(data[,2])
    if(num >= nrow(data)){
        warning("Number of output is higher than or equal to the amount of data!")
        return(as.list(data[,1]))
    }
    else{
        data <- data[order(data[,2],decreasing = TRUE),]

        size = ceiling(sum(data[,2])/num)
        top = list()
        if(size < max(data[,2])){
            warning("Some weights are higher than the capcity, will be seperated as individual files")
            top = as.list(data[data[,2] > size,1])
            data = data[-which(data[,2] > size),]
        }
        bins = binpacking(weights = data[,2],cap = size)$xbins

        ids <- lapply(unique(bins),function(x){
            id = which(bins == x)
            return(data[id,1])
        })
        ids <- c(top,ids)
        cat("There will be ",length(ids)," files output!\n")
        return(ids)
    }
}


#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+
#                      integration and filtering                      #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+

### filter erronous barcodes
bc <- read.table(args[1],stringsAsFactors = F)
colnames(bc) <- c("id","read_name","barcode","start","edit")
bc <- bc[!duplicated(bc$read_name),]

bc_edit = bc %>% group_by(barcode) %>% summarise(edit = mean(edit))
preserve = bc_edit[bc_edit$edit < 1.5,]$barcode
bc = bc[bc$barcode %in% preserve,]
bc = bc[,-1]

### merge barcode with exons
exon <- read.table(args[2],fill = FALSE)
colnames(exon) <- c("read_name","seq","gene","isoform","state","polyA","strand")
cat("There are ",nrow(exon)," reads aligned,")

cell_exon  <- left_join(bc,exon,by = "read_name")

cell_exon <- as.data.frame(na.omit(cell_exon))
cat(nrow(cell_exon)," reads identified with one or more barcodes\n")

### filter out duplicated reads
read_count = table(cell_exon$read_name)
if(sum(read_count > 1)>0){
    cell_exon_uniq = cell_exon[cell_exon$read_name %in% names(read_count[read_count == 1]),]
    cell_exon_dup = cell_exon[cell_exon$read_name %in% names(read_count[read_count > 1]),]
    cat("There are ",sum(read_count > 1),
        " reads with duplicates to ",nrow(cell_exon_dup),"reads in total\n")

    infer_barcode = substr(cell_exon_dup$seq,
                    round(cell_exon_dup$start), round(cell_exon_dup$start)+BARCODE_LEN-1)
    infer_adist = sapply(1:nrow(cell_exon_dup),function(i){
        barcode = cell_exon_dup[i,"barcode"]
        dis = adist(barcode,infer_barcode[i])
        return(dis)
    })
    cell_exon_dup = cell_exon_dup[cell_exon_dup$edit == infer_adist,]

    cell_exon_dup = cell_exon_dup[order(-cell_exon_dup$state,
                                    cell_exon_dup$edit,
                                    cell_exon_dup$gene),]
    cell_exon_dup = cell_exon_dup[!duplicated(cell_exon_dup$read_name),]
    cat("Duplicated reads are filtered, with ",nrow(cell_exon_dup)," reads left\n")
    cell_exon = rbind(cell_exon_uniq,cell_exon_dup)
}
cat("There are ",nrow(cell_exon)," reads left after deduplication, ")

### distribute genes to different files for parallelization
genes = table(cell_exon$gene)
genes = as.data.frame(cbind(names(genes),genes))
colnames(genes) = c("gene","count")
genes$count = as.numeric(genes$count)

cat("including ",nrow(genes)," genes\n")

file_num = as.numeric(args[3])
gene_list = distribute(genes,file_num)

cache <- sapply(1:length(gene_list),function(i){
  sub_genes = gene_list[[i]]
  sub_cell_exon = cell_exon[cell_exon$gene %in% sub_genes,]
  file_name = paste(c("sub_cell_exon",i,"txt"),collapse = ".")
  file_name = paste(args[4],file_name,sep = "/")

  write.table(sub_cell_exon,file = file_name,
              row.names = FALSE,col.names = TRUE,quote = FALSE)
})
