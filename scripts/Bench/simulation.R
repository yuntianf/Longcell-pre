library(stats4)
library(Rcpp)
library(e1071)
library(dplyr)
library(reshape2)
library(tidyr)

### functions ###
normalize <- function(x){
    return(log(x/sum(na.omit(x))*10000 + 1))
}

MSE <- function(x,y,if_norm = TRUE){
    if(if_norm){
        x = normalize(x)
        y = normalize(y)
    }
    return(mean((x-y)^2))
    
}

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


PCR_amplify <- function(GC_content,ratio = 20){
    m = ratio*(1.2-sigmoid((GC_content - 50)/5))
    fold = sapply(m,function(x){
        v = mean2var(x)
        size = x^2/(v-x)
        temp = rnbinom(n = 1,size = size,mu = x)
        return(temp)
    })
    
    return(fold)
}


truncate <- function(reads,trun_prob = 0.6,trunc_ratio = 0.05,split = "|",sep = ","){
    if_trunc = runif(1) < trun_prob
    if(if_trunc){
        trunc_base = rgeom(1,trunc_ratio)
        
        exons = unlist(strsplit(reads,split = split,fixed = T))
        exons = strsplit(exons,split = sep)
        
        exons_len = sapply(exons,function(x){
            x = as.numeric(x)
            return(x[2]-x[1])
        })
        
        exons_len_sum = cumsum(exons_len)
        cut_off = trunc_base < exons_len_sum
        
        if(sum(cut_off) == 0){
            return(NA)
        }
        else{
            cut_off = min(which(cut_off))
        }
        new_start = as.numeric(exons[[cut_off]][2])-(exons_len_sum[cut_off]-trunc_base)+1
        exons[[cut_off]][1] = new_start
        
        exons = exons[cut_off:length(exons)]
        reads = sapply(exons,function(x) paste(x,collapse = sep))
        reads = paste(reads,collapse = split)               
    }
    
    return(reads)
}

isoform_len <- function(isoform,split = "|",sep = ","){
    exons = unlist(strsplit(isoform,split = split,fixed = T))
    exons = strsplit(exons,split = sep)
    
    exon_len <- unlist(lapply(exons,function(i){
        i = as.numeric(i)
        return(i[2] - i[1] + 1)
    }))
    return(sum(exon_len))
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


GC_ratio <- function(len,len_up,len_low,GC_up = 1,GC_low = 0.2){
    ratio = GC_up - (len-len_low)*(GC_up-GC_low)/(len_up-len_low)
    return(ratio*100)
}

wrong_mapping <- function(reads,prob = 0.3,len_thresh = 100,split = "|",sep = ","){
    if(runif(1) < prob){
        exons = unlist(strsplit(reads,split = split,fixed = T))
        exons = strsplit(exons,split = sep)
        
        exons_len = sapply(exons,function(x){
            x = as.numeric(x)
            return(x[2]-x[1])
        })
        
        if(length(exons_len) >= 2 & sum(exons_len <= len_thresh) > 0){
            loc = max(which(exons_len <= len_thresh))
            exons[[loc+1]][1] = as.numeric(exons[[loc+1]][1])-exons_len[loc]
            exons = exons[-loc]

            reads = sapply(exons,function(x) paste(x,collapse = sep))
            reads = paste(reads,collapse = split)  
        }
    }
    return(reads)
}

iso_exprs_matrix <- function(data,cell = "cell",isoform = "isoform",size = "size"){
    data <- as.data.frame(data)
    
    cells = unique(data[,cell])
    isoforms = unique(data[,isoform])
    
    iso_exprs_count = lapply(cells,function(x){
        sub_data = data[data[,cell] == x,]
        exprs = nrow(sub_data)
        exprs_pcr = sum(sub_data[,size])
        
        iso_exprs = table(sub_data[,isoform])
        iso_vec = rep(0,length(isoforms))
        names(iso_vec) <- isoforms
        iso_vec[names(iso_exprs)] = iso_exprs
        
        iso_pcr = table(rep(sub_data[,isoform],sub_data[,size]))
        iso_pcr_vec = rep(0,length(isoforms))
        names(iso_pcr_vec) <- isoforms
        iso_pcr_vec[names(iso_pcr)] = iso_pcr
        return(c(x,exprs,exprs_pcr,iso_vec,iso_pcr_vec))
    })

    iso_exprs_count = do.call(rbind,iso_exprs_count)
    iso_exprs_count <- as.data.frame(iso_exprs_count)
    iso_exprs_count[,2:(ncol(iso_exprs_count))] = apply(iso_exprs_count[,2:(ncol(iso_exprs_count))],2,as.numeric)
    colnames(iso_exprs_count) <- c("cell","exprs","pcr",isoforms,paste(isoforms,"pcr",sep = "-"))
    
    return(iso_exprs_count)
}

simulation <- function(barcode, exprs_size,exprs_prob,exprs,
                       exon_seq,alpha,beta,
                       pcr,R1 = "CTACACGACGCTCTTCCGATCT",inter = "TTTCTTATATGGG",
                       S = NULL, T =NULL,identifier = NULL){
    #set.seed(1234)
    umi_pool <- paste(sample(c("A","G","C","T"),20000,replace = TRUE), collapse = "")
    umi <- unique(sapply(1:2000,function(x) substr(umi_pool,x,x+9)))
    if(is.null(exprs)){
        exprs = rnbinom(length(barcode),exprs_size,exprs_prob)
    }                 
    #exprs = rnbinom(length(barcode),exprs_size,exprs_prob)
    #exprs = rep(exprs,length(barcode)) 
    b <- rep(barcode,exprs)
    umi_sample <- sample(umi,length(b),replace = FALSE)
    seq_orig_set <- paste(R1,b,umi_sample,inter,sep = "") 
                     
    exon_orig_set <- unlist(lapply(exprs,function(x) isoform_sample(exon_seq,x,alpha,beta)))
    names(exon_orig_set) <- NULL                               
    exon_set <- unlist(lapply(exon_orig_set,function(x) truncate(x,trun_prob = 0.6,trunc_ratio = 0.05)))
                       
    data <- as.data.frame(cbind(seq_orig_set,b,umi_sample,exon_orig_set,exon_set))
    data <- na.omit(data)
    rownames(data) <- NULL
    colnames(data) <- c("seq","barcode","umi","exon_orig","exon_trunc")
    
    data$len = sapply(data$exon_trunc,function(x) isoform_len(x))
    data$gc <- sapply(data$len,function(x) GC_ratio(x,len_up = 500,len_low = 400))
                       
    PCR <- PCR_amplify(data$gc,pcr)                
    PCR[is.na(PCR)] = 0
    data$pcr <- PCR+1
    
    if(is.null(S)){
        S = c(0.92,0.02,0.03,0.03)
    }
    if(is.null(T)){
    T = matrix(c(0.92,0.92,0.8,0.8,
             0.02,0.02,0.02,0.02,
             0.03,0.03,0.16,0.02,
             0.03,0.03,0.02,0.16),4,4)        
    }
    data_pcr <- apply(data,1,function(x){
        pcr = as.numeric(x["pcr"])
        
        barcode <- rep(x["barcode"],pcr)
        umi <- rep(x["umi"],pcr)
        exon_orig <- rep(x["exon_orig"],pcr)
        exon_trunc <- rep(x["exon_trunc"],pcr)
        exon_trunc2 <- sapply(exon_trunc,function(x) truncate(x,trun_prob = 0.3,trunc_ratio = 0.03))
    
        seq <- rep(x["seq"],pcr)
        seq <- sapply(seq,function(x) {
            hidden_mutate(x,S,T)
        })
                  
        out <- cbind(barcode,umi,exon_orig,exon_trunc,exon_trunc2,seq)
        return(out)
    })     
    data_pcr <- as.data.frame(do.call(rbind,data_pcr))
    data_pcr <- na.omit(data_pcr)
    rownames(data_pcr) <- NULL 
        
    data_pcr$map = sapply(data_pcr$exon_trunc2,function(x) wrong_mapping(x))    
                              
    data_pcr$id <- 1:nrow(data_pcr)
    
    seq_file = paste(c("~/Upenn/Thesis/projects/Single-cell-long-reads/scripts/Bench/search_seq",identifier,"txt"),collapse = ".")
    write.table(data_pcr[,c("id","seq")],file = seq_file,
           row.names = F,col.names = F,sep = "\t",quote = F)

    barcode_file = paste(c("~/Upenn/Thesis/projects/Single-cell-long-reads/scripts/Bench/barcodes",identifier,"txt"),collapse = ".")
    write.table(barcode,file = barcode_file,
           row.names = F,col.names = F,quote = F)
        
    match_file = paste(c("~/Upenn/Thesis/projects/Single-cell-long-reads/scripts/Bench/pos_match",identifier,"txt"),collapse = ".")
    command = paste(c("python ~/Upenn/Thesis/projects/Single-cell-long-reads/scripts/BarcodeMatch/BarcodeMatch.py -p bench -q",seq_file, "-c", 
                                         barcode_file, "-m pos -o", match_file,"-co 1"),collapse = " ")
    system(command)
    
    data_match <- read.table(match_file)
    data_match <- data_match[,-1]
    colnames(data_match) <- c("id","cell","start","edit")   
    data_match$id <- as.numeric(data_match$id)
        
    cache = sapply(c(seq_file,barcode_file,match_file),function(i){
          command = paste(c("rm",i),collapse = " ")
          system(command)
    })

    cell_exon <- left_join(data_pcr,data_match,by = "id")
    cell_exon <- na.omit(cell_exon)
    #return(list(nrow(data_pcr),cell_exon))
    return(cell_exon)
}
