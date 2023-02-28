### UMI similarity graph ###
umi_sim_graph_cpp <- function(umi,iso = NULL,sim_thresh = 5,
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
  #return(iso_umi_table)
  umi_ns = umi_graph_table(iso_umi_table$umi,iso_umi_table$iso,iso_umi_table$count,
                           sim_thresh,iso_thresh,split,sep);
  
  umi_ns = as.data.frame(do.call(rbind,umi_ns))
  colnames(umi_ns) <- c("node1","node2","ns","capacity")
  
  umi_ns = umi_ns[umi_ns$ns > 0,]
  umi_ns$node1 = umi_ns$node1+1
  umi_ns$node2 = umi_ns$node2+1

  graph = graph_from_data_frame(umi_ns,directed = FALSE,
                                vertices = 1:nrow(iso_umi_table))
  graph = igraph::simplify(graph,remove.loops = FALSE,edge.attr.comb = "first")
  graph = graph %>%
    set_vertex_attr("count", value = iso_umi_table$count)
  
  umi_corres_size = sapply(iso_umi_table$id,length)
  umi_corres = rep(1:nrow(iso_umi_table),umi_corres_size)
  umi_corres = as.data.frame(cbind(do.call(c,iso_umi_table$id),umi_corres))
  colnames(umi_corres) = c("id","pair_id")
  return(list(graph,umi_corres))
}

createSNN <- function(dis,count = NULL,self = FALSE){
    dis <- as.matrix(dis)
    dis <- rbind(dis,dis[,c(2,1)])
    
    dis <- as.data.frame(dis)
    colnames(dis) <- c("node1","node2")
    
    if(is.null(count)){
      count = rep(1,length(unique(dis$node1)))
      names(count) = unique(dis$node1)
    }
    
    neighbors <- as.data.frame(dis %>% group_by(node1) %>% 
                               summarise(neighbor = list(node2)))
    snn <- lapply(1:nrow(neighbors),function(x){
        x_neighbor = unlist(neighbors[x,"neighbor"])
        if(!self){
            x_neighbor = setdiff(x_neighbor,neighbors[x,"node1"])
        }
        share_nodes = sapply(x_neighbor,function(y){
            y_neighbor = unlist(neighbors[neighbors$node1 == y,"neighbor"])
            if(!self){ 
                y_neighbor = setdiff(y_neighbor,y)
            }
            share = intersect(x_neighbor,y_neighbor)
            
            if(length(share) > 0){
              share = sum(count[share])
            }
            else{
              share = 0
            }
            return(c(neighbors[x,"node1"],y,share))
        })
        share_nodes <- as.data.frame(t(share_nodes))
        return(share_nodes)
    })
    
    snn <- as.data.frame(do.call(rbind,snn))
    colnames(snn) <- c("node1","node2","share")
    return(snn)
}

SNN_graph <- function(graph){    
    vertex = vertex_attr(graph)$name
    count = vertex_attr(graph)$count
    names(count) = vertex
    
    edge_table = as_long_data_frame(graph)[,c("from","to")]
    edge_table = cbind(vertex[edge_table$from],vertex[edge_table$to])
    
    snn = createSNN(dis = edge_table,count = count,self = FALSE) 
    snn <- snn[snn$node1 <= snn$node2,]

    index <- paste(snn$node1,snn$node2,sep = "|")
    edge_attr(graph, "share", index = index) <- as.numeric(snn$share)
    edge_attr(graph)$share[is.na(edge_attr(graph)$share)] = 1
    
    return(graph)
}

### cluster_function ###
louvain_iter_stack <- function(graph,weight = "ns",share = "share",
                         resolution = 1,alpha = 3, sim_thresh = 6){
  in_graph_list = decompose(graph)
  out_graph_list = list()
  
  while(length(in_graph_list) > 0){
    temp_graph = in_graph_list[[1]]
    in_graph_list[[1]] = NULL
    
    if(length(V(temp_graph)) == 1){
      out_graph_list = append(out_graph_list,list(temp_graph))
    }
    
  else{
    min_cut <- graph.mincut(temp_graph)
    if(min_cut >= length(V(temp_graph))/alpha | 
       min(edge_attr(temp_graph)[[weight]]) > sim_thresh){
      out_graph_list = append(out_graph_list,list(temp_graph))
    }
    else{
      temp_graph = SNN_graph(temp_graph)
      
      edge_weight <- edge_attr(temp_graph)[[weight]]
      edge_share <- edge_attr(temp_graph)[[share]]
      edge_share[edge_weight > sim_thresh & edge_share == 0] <- 1
      
      edge_attr(temp_graph,"weight",index = E(temp_graph)) = edge_weight*sqrt(edge_share)
      cluster <- cluster_louvain(temp_graph,resolution = resolution)
      
      graph_cut <- delete_edges(temp_graph, E(temp_graph)[crossing(cluster,temp_graph)])
      
      sub_graphs = decompose(graph_cut)
      
      if(length(sub_graphs) == 1){
        out_graph_list = append(out_graph_list,sub_graphs)
      }
      else{
        in_graph_list = append(in_graph_list, sub_graphs)
      }
    }
  }
  }
  return(out_graph_list)
}

graph_to_cluster <- function(graph_list){
    out <- lapply(1:length(graph_list),function(i){
        cluster = graph_list[[i]]
        cluster_table <- cbind(vertex_attr(cluster)$name,i)
    })
    out <- as.data.frame(do.call(rbind,out))
    colnames(out) <- c("id","cluster")
    out$id = as.numeric(out$id)
    out$cluster = as.numeric(out$cluster)
    out <- out[order(out[,"id"]),]
    
    return(out[,"cluster"])
}

umi_cluster_cpp <- function(umi,iso = NULL,thresh = 5){
  graph_corres <- umi_sim_graph_cpp(umi,iso = iso,sim_thresh = thresh)
  graph = graph_corres[[1]]
  umi_corres = graph_corres[[2]]
  
  graph_cluster = louvain_iter_stack(graph = graph,alpha = 3,sim_thresh = thresh+2)
  cluster = graph_to_cluster(graph_cluster)
  cluster = as.data.frame(cbind(1:length(cluster),cluster))
  colnames(cluster) = c("pair_id","cluster")
  
  cluster_expand = left_join(umi_corres,cluster,by = "pair_id")
  cluster_expand = cluster_expand[order(cluster_expand$id),]
  return(cluster_expand$cluster)
}