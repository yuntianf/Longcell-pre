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

### UMI edit distance table ###
frNN_edit <- function(dis,n,eps = 2){
  dis <- as.data.frame(dis)
  colnames(dis) <- c("node1","node2","dis")
  dis <- dis[order(dis$node1),]
  out_frNN <- frNN(as.dist(1), eps = eps)
  out_frNN$dist <- lapply(1:n,function(i){
    return(dis[dis$node1 == i,"dis"])
  })
  out_frNN$id <- lapply(1:n,function(i){
    return(dis[dis$node1 == i,"node2"])
  })
  out_frNN$sort <- F
  out_frNN <- frNN(out_frNN,eps = eps)
  return(out_frNN)
}

umi_edit_cpp <- function(umi,iso = NULL,edit_thresh = 2, k = 10,
                            iso_thresh = 80,split = "|",sep = ","){
  if(is.null(iso)){
    iso = rep("N",length(umi))
  }
  else if(length(umi) != length(iso)){
    warning("The length of umi and isoforms don't correspond, 
                the information of isoforms won't be used!")
    iso = rep("N",length(umi))
  }
  
  iso_umi_table = as.data.frame(cbind(iso,umi))
  colnames(iso_umi_table) = c("iso","umi")
  iso_umi_table$id = 1:nrow(iso_umi_table)
  iso_umi_table = iso_umi_table %>% 
    group_by(iso,umi) %>% 
    summarise(count = n(),id = list(id),.groups = "drop")
  
  umi_edit = umi_edit_table(iso_umi_table$umi,iso_umi_table$iso,
                            edit_thresh,k,iso_thresh,split,sep);
  
  umi_edit = as.data.frame(do.call(rbind,umi_edit))
  if(nrow(umi_edit) == 0){
    frNN_obj = NULL
  }
  else{
    colnames(umi_edit) <- c("node1","node2","edit")
  
    umi_edit$node1 = umi_edit$node1+1
    umi_edit$node2 = umi_edit$node2+1
  
    frNN_obj = frNN_edit(umi_edit,nrow(iso_umi_table),edit_thresh)
  }
  
  umi_corres_size = sapply(iso_umi_table$id,length)
  umi_corres = rep(1:nrow(iso_umi_table),umi_corres_size)
  umi_corres = as.data.frame(cbind(do.call(c,iso_umi_table$id),umi_corres))
  
  colnames(umi_corres) = c("id","pair_id")
  return(list(frNN_obj,iso_umi_table$count,umi_corres))
}


umi_cluster_dbscan <- function(umi,iso = NULL,thresh = 2,k = 10,
                               iso_thresh = 80,split = "|",sep = ","){
  edit_corres <- umi_edit_cpp(umi,iso = iso,edit_thresh = thresh,
                              iso_thresh = iso_thresh,split = split,sep = sep)
  frNN_ob = edit_corres[[1]]
  umi_count = edit_corres[[2]]
  umi_corres = edit_corres[[3]]
  
  if(!is.null(frNN_ob)){
    dbscan_cluster = dbscan::dbscan(frNN_ob,minPts = 1, weights = umi_count)
    cluster = dbscan_cluster$cluster
  }
  else{
    cluster = 1:length(umi_count)
  }
  cluster = as.data.frame(cbind(1:length(umi_count),cluster))
  colnames(cluster) = c("pair_id","cluster")
  
  cluster_expand = left_join(umi_corres,cluster,by = "pair_id")
  cluster_expand = cluster_expand[order(cluster_expand$id),]
  return(cluster_expand$cluster)
}