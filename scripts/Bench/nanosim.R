library(stats4)
library(Rcpp)
library(e1071)
library(dplyr)
library(reshape2)
library(tidyr)

### functions ###
hidden_mutate <- function(string,S,T){
    #set.seed(1234)
    seq <- unlist(strsplit(string,""))
    dic = c("A","G","C","T")

    S = cumsum(S)
    T = t(apply(T,1,cumsum))

    i = 1
    prob = runif(1)
    state = sum(prob > S)+1
    #print(paste(i,state,seq[i],sep = ":"))
    switch(state,
          {i = i + 1},
          {seq[1] = sample(setdiff(dic,seq[1]),1)
           i = i + 1},
          {seq <- append(seq,sample(dic,1),1)},
          {seq <- seq[-1]})

    while(i <= length(seq)){
        prob = runif(1)
        state = sum(prob > T[state,])+1
        #print(paste(i,state,seq[i],sep = ":"))
        switch(state,
          {i = i + 1},
          {seq[i] = sample(setdiff(dic,seq[i]),1)
           i = i + 1},
          {seq <- append(seq,sample(dic,1),i)},
          {seq <- seq[-i]})
    }
    return(paste(seq,collapse = ""))
}

mean2var <- function(m){
  return(2.566*m-2.566)
}

PCR_amplify <- function(GC_content,ratio = 20){
    #m = ratio*(1.2-sigmoid((GC_content - 50)/5))
    m = ratio*(sigmoid(0.5*(GC_content - 50)))+2
    fold = sapply(m,function(x){
        v = mean2var(x)
        size = x^2/(v-x)
        temp = rnbinom(n = 1,size = size,mu = x)
        return(temp)
    })
    fold[fold < 1] = 1
    return(round(fold))
}

truncate <- function(reads,trun_prob = 0.6,trunc_ratio = 0.05){
    if_trunc = runif(1) < trun_prob
    if(if_trunc){
        trunc_base = rgeom(1,trunc_ratio)
        reads = substr(reads,trunc_base,nchar(reads))
    }

    return(reads)
}

isoform_sample <- function(isoforms,exprs,alpha,beta){
    if(length(isoforms) > 2){
        warning("Only two isoforms are supported in the simulation! Only the first two will be used.")
        isoforms = isoform[1:2]
    }
    prob = rbeta(1,alpha,beta)
    isoform_1 = rbinom(1,exprs,prob)
    isoform_2 = exprs-isoform_1

    return(rep(isoforms,c(isoform_1,isoform_2)))
}


GC_ratio <- function(reads){
  ratio = sapply(reads,function(x){
    bases = unlist(strsplit(x,split = ""))
    base_count = table(bases)
    ratio = sum(base_count[c("G","C")])/nchar(x)
  })
  names(ratio) = NULL
  return(ratio*100)
}

transcript_simulator <- function(transcript,trun_prob = 0.6,trunc_ratio = 0.05){
  transcript = sapply(transcript,function(x){
    truncate(x,trun_prob,trunc_ratio)
  })
  return(transcript)
}

softclips_simulator <- function(barcodes,exprs,
                               R1 = "CTACACGACGCTCTTCCGATCT",inter = "TTTCTTATATGGG",
                               toolkit = 5){
  umi_pool <- paste(sample(c("A","G","C","T"),20000,replace = TRUE), collapse = "")
  umi_start = runif(5000,1,20000-10)
  umi <- unique(sapply(umi_start,function(x) {
    substr(umi_pool,x,x+9)
  }))

  if(length(barcodes) != length(exprs)){
    warning("The length of the barcodes is different from the length of exprs!")
  }
  barcodes <- rep(barcodes,exprs)
  umi_sample <- sample(umi,length(barcodes),replace = TRUE)

  if(toolkit == 5){
    softclips = paste(R1,barcodes,umi_sample,inter,sep = "")
  }
  else if(toolkit == 3){
    softclips = paste(umi_sample,barcodes,inter,sep = "")
  }
  else{
    stop("Invalid toolkit, only 5 or 3!")
  }

  return(softclips)
}

reads_amplify <- function(reads,ratio = 20){
  GC = sapply(reads,function(x){
    GC_ratio(x)
  })

  fold = sapply(GC,function(x){
    PCR_amplify(x,ratio)
  })

  return(round(fold))
}

qual_assign <- function(reads){
  qual = readRDS("/D/Necessary/Upenn/Thesis/data/P8640/RPL41_qual.rds")
  #qual = readRDS("./RPL41_qual.rds")
  qual = paste(qual,collapse = "")

  reads_len = nchar(reads)
  qual_start = runif(length(reads),1,nchar(qual)-max(reads_len))
  qual_seq = sapply(1:length(reads),function(i){
    qual_str = substr(qual,qual_start[i],qual_start[i]+reads_len[i]-1)
    return(qual_str)
  })
  return(qual_seq)
}

dataframe2fastq <- function(data,file,
                            name_col = "read_name",
                            read_col = "reads",
                            qual_col = "qual"){
  data$sep = "+"
  out = data[,c(name_col,read_col,"sep",qual_col)]

  out = as.data.frame(t(out))
  out = unlist(out)
  names(out) = NULL

  write.table(out,file = file,quote = FALSE,row.names = FALSE,col.names = FALSE)
  return(NULL)
}

reads_simulator <- function(barcodes,exprs,transcripts,pcr,file = NULL,
                            transname = NULL,alpha = 1,beta = 0,
                            start=NULL,trans=NULL,
                            R1 = "CTACACGACGCTCTTCCGATCT",inter = "TTTCTTATATGGG",
                            trun_prob = 0.6,trunc_ratio = 0.05,
                            toolkit = 5){
  softclips = softclips_simulator(barcodes,exprs,R1,inter,toolkit)


  if(toolkit == 5){
    barcode_seq = substr(softclips,nchar(R1)+1,nchar(R1)+16)
    umi_seq = substr(softclips,nchar(R1)+17,nchar(R1)+26)
  }
  else if(toolkit == 3){
    barcode_seq = substr(softclips,11,11+15)
    umi_seq = substr(softclips,1,10)
  }

  transcripts = lapply(transcripts,function(x){
    x = rep(x,sum(exprs))
    x = transcript_simulator(x,trun_prob,trunc_ratio)
    return(x)
  })
  if(is.null(transname)){
    transname = paste("isoform",1:length(transcripts),sep = "_")
  }
  ratio = rbeta(length(barcodes),alpha,beta)
  transcripts_sample = lapply(1:length(barcodes),function(i){
    isoform1 = sample(transcripts[[1]],round(exprs[i]*ratio[i]))
    isoform2 = sample(transcripts[[2]],exprs[i]-round(exprs[i]*ratio[i]))

    isonames = rep(transname,c(length(isoform1),length(isoform2)))
    return(cbind(isonames,c(isoform1,isoform2)))
  })
  transcripts_sample = do.call(rbind,transcripts_sample)
  rownames(transcripts_sample) = NULL
  transcripts_sample[,2] = paste(transcripts_sample[,2],
                                 paste(rep("A",20),collapse = ""),sep = "")
  iso_name = transcripts_sample[,1]

  if(toolkit == 5){
    reads = paste(softclips,transcripts_sample[,2],sep = "")
  }
  else{
    reads = paste(R1,transcripts_sample[,2],softclips,sep = "")
  }

  #fold = reads_amplify(reads,pcr)
  GC = ifelse(iso_name == transname[1],45,55)
  fold = PCR_amplify(GC_content = GC,ratio = pcr)

  data = cbind(rep(barcode_seq,fold),rep(umi_seq,fold),
               rep(iso_name,fold),rep(reads,fold))
  data = as.data.frame(data)
  colnames(data) = c("barcode","umi","isoform","reads")

  if(is.null(trans)){
    trans = matrix(c(0.94,0.15,0.15,0.15,
                     0.02,0.35,0.25,0.25,
                     0.02,0.25,0.35,0.25,
                     0.02,0.25,0.25,0.35),4,4)
  }
  if(is.null(start)){
    start = trans[1,]
  }

  data$reads = sapply(data$reads,function(x) hidden_mutate(x,start,trans))
  data$qual = qual_assign(data$reads)

  data$read_name = paste(data$barcode,data$umi,data$isoform,1:nrow(data),sep = "_")
  data$read_name = paste("@",data$read_name,sep = "")

  if(!is.null(file)){
    cache = dataframe2fastq(data,file)
  }
  return(data)
  #return(data[,c("barcode","umi","isoform")])
}


