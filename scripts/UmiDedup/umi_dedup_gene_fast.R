# args:
# 1. cell exon matrix
# 2. umi length
# 3. outfile

library(dbscan)
library(stringi)
library(stats4)
library(tidyverse)
library(dplyr)
library(reshape2)
library(Rcpp)

sourceCpp(file = "./umi_adist.cpp")
sourceCpp(file = "./cluster_merge.cpp")

frNN_dis <- function(dis,n,eps = 0,if_direct = T){
    if(if_direct){
        dis <- rbind(dis,dis[,c(2,1,3)])
    }
    dis <- as.data.frame(dis)
    colnames(dis) <- c("umi1","umi2","dis")
    dis <- dis[order(dis$umi1),]
    
    out_frNN <- frNN(as.dist(1), eps = 5)
    out_frNN$dist <- lapply(1:n,function(i){
        return(dis[dis$umi1 == i,"dis"])
    })
    out_frNN$id <- lapply(1:n,function(i){
        return(dis[dis$umi1 == i,"umi2"])
    })
    out_frNN$sort <- F
    
    out_frNN <- frNN(out_frNN,eps = eps)
    
    return(out_frNN)
}

umi_cluster_dbscan <- function(exon,umi,coverage,k = 6,thresh = 2,sep = "|"){  
    umi_count = table(umi)
    umi_dist = umi_adist(names(umi_count),k=k)
    if(length(umi_dist) == 0){
        cluster <- cbind(1:length(umi),exon,1:length(umi))
        return(cluster)
    }
    umi_dist_frNN <- frNN_dis(umi_dist,length(umi_count),eps = 1)

    cluster = dbscan(umi_dist_frNN,minPts =1,weights = umi_count)$cluster
    names(cluster) <- names(umi_count)
    cluster_expand = cluster[umi]
    cluster = as.data.frame(cbind(1:length(umi),exon,cluster_expand))
    colnames(cluster) <- c("id","exon","cluster")

    exon_type <- lapply(unique(cluster$cluster),function(i){
        loc = which(cluster$cluster == i)
        sub_cluster <- cluster[loc,]
        
        exon_count = table(sub_cluster$exon)
        mode_exon <- names(exon_count[which(exon_count == max(exon_count))])
        if(length(mode_exon) > 1){
            mode_exon_split <- strsplit(mode_exon,split = sep,fixed = T)
            mode_exon_len <- sapply(mode_exon_split,function(i) length(i))
            mode_exon = mode_exon[which(mode_exon_len == max(mode_exon_len))]
            if(length(mode_exon) > 1){
                mode_exon_coverage <- sapply(mode_exon,function(i){
                    sub_coverage = strsplit(coverage[sub_cluster$exon == i],split = sep,fixed = T)
                    sub_coverage <- mean(sapply(sub_coverage,function(j) sum(as.numeric(j))))
                })
                mode_exon = mode_exon[which(mode_exon_coverage == max(mode_exon_coverage))][1]
            }
        }
        sub_cluster$exon = mode_exon
        return(sub_cluster)
    })
    exon_type <- as.data.frame(do.call(rbind,exon_type))
    colnames(exon_type) <- c("id","exon","cluster")
    exon_type$id <- as.numeric(exon_type$id)
    exon_type$cluster <- as.numeric(exon_type$cluster)

    exon_type_cluster <- lapply(unique(exon_type$exon),function(i){
        loc = which(exon_type$exon == i)
        sub_exon_type <- exon_type[loc,]
        if(nrow(sub_exon_type) == 1){
            sub_exon_type$cluster = 1
            return(sub_exon_type)
        }
        sub_umi_dist <- umi_dist[umi_dist[,1] %in% sub_exon_type$id & 
                                 umi_dist[,2] %in% sub_exon_type$id,]
        if(is.null(nrow(sub_umi_dist))){
            sub_umi_dist <- t(as.matrix(sub_umi_dist))
        }
        
        merge_cluster <- cluster_merge(as.matrix(sub_exon_type[,c("id","cluster")]),sub_umi_dist,thresh)
        sub_exon_type$cluster <- merge_cluster[,2]

        if(length(unique(sub_exon_type$cluster)) == 1){
            return(sub_exon_type)
        }else{
            len = table(sub_exon_type$cluster)
            retain = as.numeric(names(len)[len > min(mean(len),sd(len))/1.5])
            sub_exon_type[!(sub_exon_type$cluster %in% retain),"cluster"] = 0
            return(sub_exon_type)
        }  
    })
    exon_type_cluster <- as.data.frame(do.call(rbind,exon_type_cluster))
    colnames(exon_type_cluster) <- c("id","exon","cluster")
    return(exon_type_cluster)
}

coverage_summary<-function(exon,exon_adjust,exon_coverage,sep = "|"){
    exon_coverage = exon_coverage[exon == exon_adjust]
    exon_coverage <- strsplit(exon_coverage,split = sep,fixed = T)
    exon_coverage <- do.call(rbind,exon_coverage)   
    exon_coverage <- apply(exon_coverage,2,function(x) median(as.numeric(x)))
    exon_coverage_str <- paste(exon_coverage,collapse = sep)
    return(exon_coverage_str)
}

umi_count <- function(cell_exon,bar = "barcode",start = "start",gene = "gene_ID",umi = NULL,
                      seq = "search_seq",exon = "exon_seq",coverage = "coverage",
                      reads_start = "reads_s",reads_end = "reads_e",polyA = "polyA",sep = "|",
                      bar_len = 16, flank = 1,umi_len = 10){
    colnames(cell_exon)[which(colnames(cell_exon) == bar)] = "cell"
    colnames(cell_exon)[which(colnames(cell_exon) == gene)] = "gene"
    colnames(cell_exon)[which(colnames(cell_exon) == exon)] = "exon"
    colnames(cell_exon)[which(colnames(cell_exon) == coverage)] = "coverage"
    colnames(cell_exon)[which(colnames(cell_exon) == reads_start)] = "reads_s"
    colnames(cell_exon)[which(colnames(cell_exon) == reads_end)] = "reads_e"
    colnames(cell_exon)[which(colnames(cell_exon) == polyA)] = "polyA"
    
    cell_exon[,start] <- as.numeric(cell_exon[,start])
    bar_len <- as.numeric(bar_len)
    umi_len <- as.numeric(umi_len)
    flank = as.numeric(flank)
    
    if(is.null(umi)){
        umi_flank <- sapply(1:nrow(cell_exon),function(i){
            substr(cell_exon[i,seq],cell_exon[i,start] + bar_len - flank + 1,
                   cell_exon[i,start] + bar_len + umi_len + flank)
        })
        cell_exon$umi = umi_flank
    }
    
    search_len = umi_len + 2*flank
    cell_exon <- cell_exon[nchar(cell_exon$umi) == search_len,]
    
    cell <- unique(cell_exon[,"cell"])
    cell_gene_exon_count <- lapply(cell,function(i){
        cell_i = cell_exon[cell_exon[,"cell"] == i,]
        genes = unique(cell_i[,"gene"])
        cat("There are",length(genes), "genes in cell",i,"\n")

        gene_exon_count <- lapply(genes,function(j){
            print(j)
            cell_i_gene_j = cell_i[cell_i[,"gene"] == j,]
            cell_i_gene_j$cluster = 0
            if(nrow(cell_i_gene_j) != 1){
                cluster = as.data.frame(umi_cluster_dbscan(cell_i_gene_j$exon,
                                                           cell_i_gene_j$umi,
                                                           cell_i_gene_j$coverage))
                colnames(cluster) <- c("id","exon","cluster")
                cluster$id <- as.numeric(cluster$id)
                cluster$cluster <- as.numeric(cluster$cluster)
                cell_i_gene_j$cluster[cluster$id] = cluster$cluster
                cell_i_gene_j$exon_adjust[cluster$id] = cluster$exon
            }
            else{
                cell_i_gene_j$exon_adjust = cell_i_gene_j$exon
                cell_i_gene_j$cluster = 1
            }
            
            cell_i_gene_j = cell_i_gene_j[cell_i_gene_j$cluster > 0,]
            
            cell_i_gene_j_dedup <- cell_i_gene_j %>% 
                                   group_by(cell,gene,exon_adjust,cluster) %>%
                                   summarise(pcr = n(),
                                             coverage = coverage_summary(exon,exon_adjust,coverage,sep),
                                             reads_s = median(reads_s),reads_e = median(reads_e),
                                             polyA = mean(as.logical(polyA)),.groups = 'drop')
            
            return(cell_i_gene_j_dedup)
        })
        
        gene_exon_count[vapply(gene_exon_count,is.null,logical(1L))] <- NULL
        gene_exon_count = as.data.frame(do.call(rbind,gene_exon_count))
        return(as.data.frame(gene_exon_count))
    })
    
    cell_gene_exon_count[vapply(cell_gene_exon_count,is.null,logical(1L))] <- NULL
    cell_gene_exon_count = as.data.frame(do.call(rbind,cell_gene_exon_count))
    rownames(cell_gene_exon_count) <- NULL
    return(cell_gene_exon_count[,-which(colnames(cell_gene_exon_count) == "cluster")])
}


args <- commandArgs(trailingOnly = TRUE)

cell_exon <- read.table(args[1])
colnames(cell_exon) <- c("read_name","barcode","start","gene_ID","exon_seq","search_seq","polyA","strand","coverage","reads_s","reads_e")

if(nrow(cell_exon)==0){
       stop("This cell doesn't exsit!")
}

cell_exon_count <- umi_count(cell_exon, umi_len=args[2])

write.table(cell_exon_count,file = args[3],row.names = F,col.names = F,quote = F,sep = "\t")
