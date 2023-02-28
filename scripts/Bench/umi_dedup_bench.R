### args:
# 1.data quality 
# 2.expression
# 3.out file name

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
source("~/scripts/UMI_dedup/dedup/splice_site_correct.R")
source("~/scripts/UMI_dedup/dedup/umi_count_functions.R")

wl = read.table("~/data/BarcodeMatch/FT/FT_softclips_exon/FT_cells.txt")
wl <- (unlist(c(wl)))
names(wl) <- NULL

test_data <- function(exprs,pcr,S,T,exon_seq, identifier= NULL){
    cell_group1 <- simulation(barcode = wl[1:50],exprs = NULL,
                              exprs_size = exprs,exprs_prob = 0.75,
                              exon_seq,alpha = 50,beta = 50,
                              pcr = pcr,S = S, T = T,identifier= identifier)
    cell_group2 <- simulation(barcode = wl[51:100],exprs = NULL,
                              exprs_size = exprs,exprs_prob = 0.75,
                              exon_seq,alpha = 0.5,beta = 0.5,
                              pcr = pcr,S = S, T = T,identifier= identifier)
    
    sim_cells <- as.data.frame(rbind(cell_group1,cell_group2))  
    sim_cells$polyA = TRUE
    colnames(sim_cells)[2] <- "orig_umi"
    
    return(sim_cells)
}

psi_test <- function(exons,count){
    psi = sum(count[grepl("1400",exons)])/sum(count[grepl("1200",exons)])
    return(psi)
}

umi_dedup_test <- function(data){
    sim_real <- data %>% 
        group_by(barcode,exon_orig) %>% 
        summarise(size = n(),count = length(unique(orig_umi)),.groups = "drop")
        
    pcr = mean(sim_real$size/sim_real$count)
    
    sim_real_summary <- sim_real %>% 
                        group_by(barcode) %>% 
                        summarise(pcr = sum(size),
                                  exprs = sum(count),
                                  psi_pcr = psi_test(exon_orig,size),
                                  psi = psi_test(exon_orig,count))
    exprs = mean(sim_real_summary$exprs)
    
    sim_dedup <- gene_umi_count(data,data$edit,
                                bar = "cell",seq = "seq",isoform = "map")
    sim_dedup_filter = sim_dedup[sim_dedup$count > 0,]
    es_pcr = mean(sim_dedup_filter$size/sim_dedup_filter$count)
        
    sim_dedup_summary <- sim_dedup %>% 
                         group_by(cell) %>% 
                         summarise(es_cluster = sum(cluster),
                                   es_exprs = sum(count),
                                   es_cluster_psi = psi_test(isoform,cluster),
                                   es_psi = psi_test(isoform,count))
    es_exprs = mean(sim_dedup_summary$es_exprs)
    
    comp = left_join(sim_dedup_summary,sim_real_summary,by = c("cell" = "barcode"))
    comp = na.omit(comp)
    
    ### exprs evaluation
    pcr_exprs_diff = MSE(comp$pcr,comp$exprs)
    cluster_exprs_diff = MSE(comp$es_cluster,comp$exprs)
    dedup_exprs_diff = MSE(comp$es_exprs,comp$exprs)
    
    ### psi evaluation    
    pcr_psi_diff = MSE(comp$psi,comp$psi_pcr,if_norm = F)
    cluster_psi_diff = MSE(comp$psi,comp$es_cluster_psi,if_norm = F)
    dedup_psi_diff = MSE(comp$psi,comp$es_psi,if_norm = F)
    
    return(c(pcr,exprs,es_pcr,es_exprs,
             pcr_exprs_diff,cluster_exprs_diff,dedup_exprs_diff,
             pcr_psi_diff,cluster_psi_diff,dedup_psi_diff))
}

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

args <- commandArgs(trailingOnly = TRUE)

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


times = 5
PCR = seq(2,50,4)
exprs = as.numeric(args[2])

exon_seq = c("1000,1200|1400,1500|1600,1800","1000,1200|1600,1800")

error_table <- lapply(1:times,function(x){
    error = lapply(PCR,function(i){
        identifier=paste(c(qual,exprs,i),collapse = "_")
        sim = test_data(exprs = exprs,pcr = i,
                        S = S,T=T,exon_seq = exon_seq,identifier = identifier)
        sub_error = umi_dedup_test(sim)
    })
    error <- do.call(rbind,error)
    return(error)
})
error_table <- do.call(rbind,error_table)
colnames(error_table) = c("pcr","exprs","es_pcr","es_exprs","pcr_exprs_diff","cluster_exprs_diff",
"dedup_exprs_diff","pcr_psi_diff","cluster_psi_diff","dedup_psi_diff")
write.table(error_table,file = args[3], row.names = FALSE,col.names = TRUE,quote = FALSE)
