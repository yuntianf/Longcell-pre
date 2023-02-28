# args:
# 1. cell exon matrix
# 2. cells
# 4. outfile



library(dbscan)
library(stringi)
library(stats4)
library(tidyverse)
library(dplyr)

kmer <- function(string,k){
    n = nchar(string)
    kmers <- sapply(1:(n-k+1),function(x) substr(string,x,x+k-1))
    return(kmers)
}

umi_adist <- function(umi_flank,k = 5){
    umi_bin <- sapply(umi_flank,function(y) kmer(y,10))
    umi_dist = apply(umi_bin,2,function(y) apply(umi_bin,2,function(z) min(adist(y,z))))
    colnames(umi_dist) = NULL
    rownames(umi_dist) = NULL 
                                                  
    umi_kmer <- sapply(umi_flank,function(x) kmer(x,k))
    umi_kmer_index <- unique(c(umi_kmer))
                       
    umi_kmer_index <- sapply(umi_kmer_index,function(x) 
    which(umi_kmer==x)%/%nrow(umi_kmer) + (which(umi_kmer==x)%%nrow(umi_kmer) > 0)) 
    
    if(mode(umi_kmer_index) == "numeric"){
        umi_kmer_index <- as.list(as.data.frame(umi_kmer_index))
    }                         
                             
    umi_kmer_index <- umi_kmer_index[unlist(lapply(umi_kmer_index,function(x) length(x) > 1))]   
                                                   
    dist_filter <- matrix(0,length(umi_flank),length(umi_flank))
                                                   
    for(i in 1:length(umi_kmer_index)){
        dist_filter[umi_kmer_index[[i]],umi_kmer_index[[i]]] = 1
    } 
    umi_dist[which(dist_filter != 1)] = Inf    
    return(umi_dist)                                               
}

cluster_merge <- function(cluster,dis,thresh = 2){
    if(length(cluster) == 1){
        return(cluster)
    }
    cluster_dis <- lapply(cluster,function(x) lapply(cluster,function(y){
        if(sum(dis[x,y] == Inf) > length(dis[x,y])/5){
            return(Inf)
        }
        else{
            valid = mean(c(dis[x,y])[which(c(dis[x,y]) != Inf)])
        }
    }))
                          
    cluster_dis <- matrix(unlist(cluster_dis),length(cluster),length(cluster))  
    cluster_dis[which(is.na(cluster_dis))] <- Inf
    diag(cluster_dis) <- Inf                      
                          
    if(min(cluster_dis) < thresh){
        loc = which(cluster_dis == min(cluster_dis))[1]
        if(loc %% nrow(cluster_dis) == 0){
            cluster1 = loc %/% nrow(cluster_dis)
            cluster2 = nrow(cluster_dis)
        }
        else{
            cluster1 = loc %/% nrow(cluster_dis) + 1
            cluster2 = loc %% nrow(cluster_dis)
        }
        cluster[[cluster1]] <- c(cluster[[cluster1]],cluster[[cluster2]])
        cluster[[cluster2]] <- NULL
        return(cluster_merge(cluster,dis,thresh))
    }
    else{
        return(cluster)
    }                      
}

umi_cluster_dbscan <- function(umi,k = 6,thresh = 2){
    umi_dist = umi_adist(umi,k=k)
    
    cluster = dbscan(as.dist(umi_dist),eps = 0,minPts =1)$cluster
    cluster = sapply(1:max(cluster),function(x) which(cluster == x))
                   
    if(mode(cluster) != "list"){
        if(is.null(nrow(cluster))){
            cluster <- as.list(cluster)
        }
        else{
            cluster <- as.list(as.data.frame(cluster)) 
        }
    }
    cluster <- cluster_merge(cluster,umi_dist,thresh)                 
                     
    len = unlist(lapply(cluster,length))
    retain = len > min(mean(len),sd(len))/1.5
                     
    return(cluster[retain])
}

umi_count <- function(cell_exon,bar = "barcode",start = "start",gene = "gene_ID",umi = NULL,
                      seq = "search_seq",exon = "exon_seq",bar_len = 16, flank = 2,umi_len = 10){
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
    
    cell <- unique(cell_exon[,bar])
    cell_gene_exon_count <- lapply(cell,function(i){
        print(i)
        cell_i = cell_exon[cell_exon[,bar] == i,]
        genes = unique(cell_i[,gene])
        
        gene_exon_count <- sapply(genes,function(j){
            cell_i_gene_j = cell_i[cell_i[,gene] == j,]
            
            isoform = unique(cell_i_gene_j[,exon])
            isoform_count <- sapply(isoform,function(z){
                cell_i_gene_j_iso_z = cell_i_gene_j[cell_i_gene_j[,exon] == z,]
                if(nrow(cell_i_gene_j_iso_z) == 1){
                    dedup_count = 1             
                }
                else{
                    cluster = umi_cluster_dbscan(cell_i_gene_j_iso_z$umi)
                    dedup_count = length(cluster)
                }
                return(c(z,nrow(cell_i_gene_j_iso_z),dedup_count))
            })    
        })
        
    if(mode(gene_exon_count) == "list"){
        temp = as.data.frame(t(do.call(cbind,gene_exon_count)))
        colnames(temp) <- c("exon","pcr","dedup")
        temp$gene <- rep(names(gene_exon_count),lapply(gene_exon_count,ncol))
        rownames(gene_exon_count) <- NULL
        gene_exon_count <- temp
    }
    else{
        gene <- colnames(gene_exon_count)
        gene_exon_count <- as.data.frame(t(matrix(gene_exon_count,nrow = 3)))
        colnames(gene_exon_count) <- c("exon","pcr","dedup")
        gene_exon_count$gene <- gene
    }
    rownames(gene_exon_count) = NULL    
    return(as.data.frame(gene_exon_count))
    })
    
    cell_name <- rep(cell,unlist(lapply(cell_gene_exon_count,nrow)))
    temp = as.data.frame(do.call(rbind,cell_gene_exon_count))
    temp$cell <- cell_name
    temp$pcr <- as.numeric(temp$pcr)
    temp$dedup <- as.numeric(temp$dedup)
    
    return(temp[,c(5,4,1,2,3)])
}


args <- commandArgs(trailingOnly = TRUE)

cell_exon <- readRDS(args[1])

cell_exon_count <- umi_count(cell_exon[cell_exon$barcode == args[2],], umi_len=args[3])

write.table(cell_exon_count,file = args[4],row.names = F,col.names = F,quote = F,sep = "\t")
