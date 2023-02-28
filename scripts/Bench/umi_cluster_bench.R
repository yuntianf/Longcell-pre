### args:
# 1. data quality
# 2. out file name

library(DNABarcodes)
library(stats4)
library(Rcpp)
library(dbscan)
library(e1071)
library(dplyr)
library(reshape2)
library(tidyr)
library(transport)
library(RcppHungarian)
library(NameNeedle)
library(igraph)
library(MASS)
source("~/scripts/UMI_dedup/bench/simulation.R")
source("~/scripts/UMI_dedup/dedup/umi_cluster.R")
source("~/scripts/UMI_dedup/dedup/reads_filter.R")

args <- commandArgs(trailingOnly = TRUE)

wl = read.table("~/data/BarcodeMatch/FT/FT_softclips_exon/FT_cells.txt")
wl <- (unlist(c(wl)))
names(wl) <- NULL

qual = as.numeric(args[1])
if(qual == 1){
    T = matrix(c(0.97,0.1,0.1,0.1,
             0.01,0.4,0.25,0.25,
             0.01,0.25,0.4,0.25,
             0.01,0.25,0.25,0.4),4,4)
}else if(qual == 2){
    T = matrix(c(0.955,0.1,0.1,0.1,
             0.015,0.4,0.25,0.25,
             0.015,0.25,0.4,0.25,
             0.015,0.25,0.25,0.4),4,4)    
}else if(qual == 3){
    T = matrix(c(0.94,0.1,0.1,0.1,
             0.02,0.4,0.25,0.25,
             0.02,0.25,0.4,0.25,
             0.02,0.25,0.25,0.4),4,4)    
}else if(qual == 4){
    T = matrix(c(0.925,0.1,0.1,0.1,
             0.025,0.4,0.25,0.25,
             0.025,0.25,0.4,0.25,
             0.025,0.25,0.25,0.4),4,4)    
}else if(qual == 5){
    T = matrix(c(0.91,0.1,0.1,0.1,
             0.03,0.4,0.25,0.25,
             0.03,0.25,0.4,0.25,
             0.03,0.25,0.25,0.4),4,4)    
}
S = T[1,]


exprs = c(5,10,20,40)
pcr = seq(2,80,4)
times = 1:5
#exprs = 5
#pcr = 2
exon_seq = c("1000,1200|1400,1500|1600,1800","1000,1200|1600,1800")


size_repeat <- lapply(times,function(t){
size <- lapply(exprs,function(i){
    cat("exprs:",i,"\n")
    sub <- lapply(pcr,function(j){
        cat("pcr",j,"\n")
        error <- try({
            identifier = paste(c(qual,i,j),collapse = "-")
            sim <- simulation(barcode = wl[1:10],exprs_size = 10,exprs_prob = 0.75,
                          exprs = i,exon_seq = rep(exon_seq[1],2),alpha = 50,
                          beta = 50,pcr = j,S = S, T = T,identifier = identifier)
        })
        if(class(error) == "try-error"){
            print(c(i,j))
            return(NULL)
        }
        edit = sim$edit
        ratio = sum(edit > 1)/length(edit)
        
        cells = table(sim$cell)
        cell_diff = abs(cells - mean(cells))
        sim = sim[sim$cell == names(cells)[which(cell_diff == min(cell_diff))][1],]
        
        if(nrow(sim) == 1){
            return(NULL)
        }
        sim$umi_flank = sapply(1:nrow(sim),function(x){
            substr(sim[x,"seq"],sim[x,"start"] + 16 - 1 + 1,
                   sim[x,"start"] + 16 + 10 + 1)
        })
        sim = sim[nchar(sim$umi_flank) == 12,]
        test <- umi_sim_graph(sim$umi_flank,sim_thresh = 5)
        
        test_cluster <- lapply(decompose(test),function(x){
            cluster = louvain_iter(x,sim_thresh = 6)
            return(cluster)
        })
        test_cluster <- unlist(test_cluster,recursive = FALSE)
        
        test_corres <- sim[,c("umi_flank","umi")]
        test_corres$id <- 1:nrow(test_corres)
        test_corres <- test_corres %>% group_by(umi) %>% summarise(id = list(id))
        
        test_corres_size = sapply(test_corres$id,length)
        test_cluster_size = sapply(test_cluster,length)
        
        filter_pcr = size_filter_pcr(test_cluster_size,relation = mean2var,alpha = 0.2)
        filter_error = size_filter_error(test_cluster_size,ratio = ratio)
        filter_total = length(test_cluster_size)-
        sum(size_filter(test_cluster_size,relation = mean2var,
                                   alpha = 0.05,ratio = ratio))
        
        return(c(nrow(sim),length(test_corres_size),length(test_cluster_size),
                filter_pcr,filter_error,filter_total))
    })
    sub <- do.call(rbind,sub)
})
size <- do.call(rbind,size)
})
size_repeat <- do.call(rbind,size_repeat)
size_repeat <- as.data.frame(size_repeat)
colnames(size_repeat) <- c("raw","real","dedup","filter_pcr","filter_error","filter_total")
write.table(size_repeat,file = args[2], row.names = FALSE,col.names = TRUE,quote = FALSE)
