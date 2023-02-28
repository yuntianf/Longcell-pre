source("./nanosim.R")
#source("/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/scripts/Bench/nanosim.R")


args <- commandArgs(trailingOnly = TRUE)
#args = c(1,10,"./","test","test.rds")

#cells =  scan("/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/bench_cells_1000.txt",
#              character(), quote = "")
cells =  scan("~/data/bench/simulation/bench_cells_1000.txt",character(), quote = "")

qual = as.numeric(args[1])
if(qual == 1){
  trans = matrix(c(0.94,0.1,0.1,0.1,
                   0.02,0.3,0.3,0.3,
                   0.02,0.3,0.3,0.3,
                   0.02,0.3,0.3,0.3),4,4)
}else if(qual == 2){
  trans = matrix(c(0.94,0.16,0.16,0.16,
                   0.02,0.28,0.28,0.28,
                   0.02,0.28,0.28,0.28,
                   0.02,0.28,0.28,0.28),4,4)
}else if(qual == 3){
  trans = matrix(c(0.94,0.25,0.25,0.25,
                   0.02,0.25,0.25,0.25,
                   0.02,0.25,0.25,0.25,
                   0.02,0.25,0.25,0.25),4,4)
}else if(qual == 4){
  trans = matrix(c(0.97,0.1,0.1,0.1,
                   0.01,0.3,0.3,0.3,
                   0.01,0.3,0.3,0.3,
                   0.01,0.3,0.3,0.3),4,4)
}else if(qual == 5){
  trans = matrix(c(0.97,0.25,0.25,0.25,
                   0.01,0.25,0.25,0.25,
                   0.01,0.25,0.25,0.25,
                   0.01,0.25,0.25,0.25),4,4)
}

transcripts = c("GGAGCCATTACTGCAGGAAAAGGTCCCGGAGAGCTGAGCAGTCAAGATGCAGTGTGACTTCACCGAAGACCAGACCGCAGAGTTCAAGGAGGCCTTCCAGCTGTTTGACCGAACAGGTGATGGCAAGATCCTGTACAGCCAGTGTGGGGATGTGATGAGGGCCCTGGGCCAGAACCCTACCAACGCCGAGGTGCTCAAGGTCCTGGGGAACCCCAAGAGTGATGAGATGAATGTGAAGGTGCTGGACTTTGAGCACTTTCTGCCCATGCTGCAGACAGTGGCCAAGAACAAGGACCAGGGCACCTATGAGGATTATGTCGAAGGACTTCGGGTGTTTGACAAGGAAGGAAATGGCACCGTCATGGGTGCTGAAATCCGGCATGTTCTTGTCACACTGGGTGAGAAGATGACAGAGGAAGAAGTAGAGATGCTGGTGGCAGGGCATGAGGACAGCAATGGTTGTATCAACTATGAAGAGCTCGTCCGCATGGTGCTGAATGGCTGAGGACCTTCCCAGTCTCCCCAGAGTCCGTGCCTTTCCCTGTGTGAATTTTGTATCTAGCCTAAAGTTTCCCTAGGCTTTCTTGTCTCAGCAACTTTCCCATCTTGTCTCTCTTGGATGATGTTTGCCGTCAGCATTCACCAAATAAACTTGCTCTCTGGG",
                "GGAGCCATTACTGCAGGAAAAGGTCCCGGAGAGCTGAGCAGTCAAGATGCAGTGTGACTTCACCGAAGACCAGACCGCAGAGTTCAAGGAGGCCTTCCAGCTGTTTGACCGAACAGGTGATGGCAAGATCCTGTACAGCCAGTGTGGGGATGTGATGAGGGCCCTGGGCCAGAACCCTACCAACGCCGAGGTGCTCAAGGTCCTGGGGAACCCCAAGAGTGATGAGATGAATGTGAAGGTGCTGGACTTTGAGCACTTTCTGCCCATGCTGCAGACAGTGGCCAAGAACAAGGACCAGGGCACCTATGAGGATTATGTCGAAGGACTTCGGGTGTTTGACAAGGAAGGAAATGGCACCGTCATGGGTGCTGAAATCCGGCATGTTCTTGTCACACTGGGTGAGAAGATGACAGAGGAAGAAGTAGAGATGCTGGTGGCAGGGCATGAGGACAGCAATGGTTGTATCAACTATGAAGCGTTTGTGAGGCATATCCTGTCGGGGTGACGGGCCCATGGGGCGGAGCTCGTCCGCATGGTGCTGAATGGCTGAGGACCTTCCCAGTCTCCCCAGAGTCCGTGCCTTTCCCTGTGTGAATTTTGTATCTAGCCTAAAGTTTCCCTAGGCTTTCTTGTCTCAGCAACTTTCCCATCTTGTCTCTCTTGGATGATGTTTGCCGTCAGCATTCACCAAATAAACTTGCTCTCTGGGCCCTCGGTTCGGT")
transname = c("MYL6-201","MYL6-202")

exprs = as.numeric(args[2])

rep = as.numeric(args[6])

data = lapply(1:rep,function(i){
  temp = reads_simulator(cells,exprs = rep(exprs,length(cells)),
                       transcripts = transcripts, transname = transname,
                       alpha = 10,beta = 10,trans = trans,
                       pcr = 20,toolkit = 5)
  temp$read_name = paste(temp$read_name,"rep","i",sep = "_")
  return(temp)
  })
data = do.call(rbind,data)

path = args[3]
prefix = args[4]
prefix = paste(path,prefix,sep = "")

file = paste(c(prefix,"fastq"),collapse = ".")
cache = dataframe2fastq(data,file)
command = paste(c("gzip ", prefix,".fastq"),collapse = "")
system(command)


saveRDS(data,file = paste(path,args[5],sep = "/"))


