### args
# 1. reads file
# 2. gene
# 3. gene bed
# 4. bam file
# 5. groups
# 6  group annotation
# 7. output folder

args <- commandArgs(trailingOnly = TRUE)

reads = read.table(args[1],header = TRUE)
gene = args[2]

### extract reads name for target genes ###
reads = reads[reads$gene == gene & !reads$V5 %in% c("5","3","8"),]

cat("There are ",nrow(reads),"mapped for gene ",gene,"\n")

### extract reads for target genes ###
bed = readRDS(args[3])
bed = bed[bed$gene == gene,]

bam = args[4]

chr = unique(as.character(bed$chr))
start = min(bed$start)
end = max(bed$end)
gene_loc = paste(c(chr,":",start,"-",end),collapse = "")

command = paste(c("samtools view -hb",bam,gene_loc,"> ./temp.bam"),collapse = " ")
print(command)
system(command)

command = "samtools index ./temp.bam"
system(command)

### extract reads name for target cells ###
groups = unlist(strsplit(args[5],split = "|",fixed = TRUE))
target = readRDS(args[6])
path = args[7]

cache <- sapply(groups,function(i){
  cells = target[target$target == i,"barcode"]
  cat("There are ",length(cells)," cells within the group ", i,",")
  sub_reads = reads[reads$barcode %in% cells,"read_name"]
  cat(" with ",length(sub_reads)," reads captured\n")
  
  file_name = paste(c("./",i,"_",gene,"_reads.txt"),collapse = "")
  write.table(sub_reads,file = file_name, col.names = FALSE,row.names = FALSE,quote = FALSE)
  
  out = paste(c(path,"/",i,"_",gene,"_reads.bam"),collapse = "")
  command = paste(c("java -jar ~/software/picard.jar FilterSamReads -I",bam,
                    "-O",out,"-READ_LIST_FILE",file_name,"-FILTER includeReadList"),collapse = " ")
  print(command)
  system(command)
  
  command = paste(c("samtools index",out),collapse = " ")
  system(command)
  
  command = paste(c("rm",file_name),collapse = " ")
  system(command)
  
  return(i)
})

command = "rm -rf ./temp.bam*"
system(command)

