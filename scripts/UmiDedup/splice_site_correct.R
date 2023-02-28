splice_site <- function(isoform,split = "|",sep = ","){
    exons = unlist(strsplit(isoform,split = split,fixed = T))
    splice = unlist(strsplit(exons,split = sep))
    
    return(splice)
}

splice_site_count <- function(isoforms,split = "|",sep = ","){
    isoforms = table(isoforms)
    splice_sites_dic = unique(unlist(lapply(names(isoforms),function(x){
        ss = splice_site(x,split,sep)
        if(length(ss) > 2){
            return(ss[2:(length(ss)-1)])
        }
        else{
            return(NA)
        }
    })))
    splice_sites_dic <- na.omit(splice_sites_dic)
    if(length(splice_sites_dic) == 0){
        return(NA)
    }
    splice_sites_count = lapply(names(isoforms),function(x){
        count = rep(0,length(splice_sites_dic))
        names(count) = splice_sites_dic
            
        site = splice_site(x,split,sep)
        if(length(site) > 2){
            site = site[2:(length(site)-1)]
            count[site] = count[site]+isoforms[x]
        }      
        return(count)
    })
    splice_sites_count <- do.call(rbind,splice_sites_count)
    splice_sites_count <- colSums(splice_sites_count)
    splice_sites_count <- splice_sites_count[order(names(splice_sites_count))]
    return(na.omit(splice_sites_count))
}

splice_site_table <- function(isoforms,split = "|",sep = ",",
                              splice_site_thresh = 10){
    ss_count = splice_site_count(isoforms,split,sep) 
    
    ss = NA
    if(!is.na(sum(ss_count))){
        ss = names(ss_count[ss_count >= splice_site_thresh])  
        if(length(ss) == 0){
            ss = NA
        }
    }

    ss_table <- lapply(isoforms,function(x){
        temp_ss = splice_site(x,split,sep)
        iso_start = temp_ss[1]
        iso_end = temp_ss[length(temp_ss)]
        
        if(!is.na(ss[1])){
            ss_vec <- rep(0,length(ss))
            names(ss_vec) <- ss
            if(length(temp_ss) > 2){
                if(sum(!temp_ss[2:(length(temp_ss)-1)] %in% ss) > 0){
                    return(NULL)
                }
                else{
                    ss_vec[temp_ss[2:(length(temp_ss)-1)]] = 1    
                }
            }
            return(c(as.numeric(iso_start),ss_vec,as.numeric(iso_end)))
        }
        else{
          if(length(temp_ss) == 2){
            return(as.numeric(c(iso_start,iso_end)))
          }
          else{
            return(NULL)
          }
        }
    })
    return(ss_table)
}

long2wide <- function(data,names_from,value_from,symmetric = TRUE){
    if(length(names_from) != 2){
        stop("Two columns are required to build the distance matrix!")
    }
    
    dimNames = unique(unlist(data[,names_from]))
    
    out = as.data.frame(matrix(0,nrow = length(dimNames),ncol = length(dimNames)))

    colnames(out) = dimNames
    rownames(out) = dimNames
    
    data[,names_from[1]] <- as.character(data[,names_from[1]])
    data[,names_from[2]] <- as.character(data[,names_from[2]])
    
    for(i in 1:nrow(data)){
        out[data[i,names_from[1]],data[i,names_from[2]]] = out[data[i,names_from[1]],data[i,names_from[2]]]+
                                                              data[i,value_from]
        if(symmetric){
            out[data[i,names_from[2]],data[i,names_from[1]]] = out[data[i,names_from[1]],data[i,names_from[2]]]
        }
    }
    return(out)
}


cluster_mid_count <- function(splice_table){  
    splice_table <- as.data.frame(splice_table)
    
    isoforms = apply(splice_table,1,function(x){
            paste(x,collapse = "")
            })
    isoform_count = table(isoforms)
    isoform_name = names(isoform_count)
    
    if(length(isoform_name) == 1){
        return(c(isoform_name,isoform_name,isoform_count))
    }
    else{
        mode_isoform = isoform_name[which(isoform_count == max(isoform_count))][1]
        isoform_coexist = cbind(isoform_name,mode_isoform,isoform_count)
        rownames(isoform_coexist) <- NULL
        return(isoform_coexist)
    }
}


cell_clusters_mid_count <- function(splice_table,cluster){    
    cluster_dic = unique(cluster)
    splice_table <- as.data.frame(splice_table) 
    
    if(nrow(splice_table) != length(cluster)){
        stop("The size of isoforms and clusters don't match!")
    }

    isoform_coexist <- lapply(cluster_dic,function(x){
        key = which(cluster == x)
        sub_splice_table = splice_table[key,]
        
        sub_isoform_coexist = cluster_mid_count(sub_splice_table)
        return(sub_isoform_coexist)
    }) 
    isoform_coexist <- as.data.frame(do.call(rbind,isoform_coexist))
    return(isoform_coexist)   
}

cells_mid_filter <- function(cells,splice_table,cluster){
    cell_dic = unique(cells)
    splice_table <- as.data.frame(splice_table)
    
    if(nrow(splice_table) != length(cells)){
        stop("The size of isoforms and cells don't match!")
    }
    
    isoform_coexist <- lapply(cell_dic,function(x){
        key = which(cells == x)
        sub_splice_table = splice_table[key,]
        sub_cluster = cluster[key]
        
        sub_isoform_coexist = cell_clusters_mid_count(sub_splice_table,sub_cluster)
        
        return(as.matrix(sub_isoform_coexist))
    })
    isoform_coexist <- as.data.frame(do.call(rbind,isoform_coexist))
    colnames(isoform_coexist) <- c("iso1","iso2","count")
    isoform_coexist$count <- as.numeric(isoform_coexist$count)
    isoform_coexist <- isoform_coexist %>% 
                       group_by(iso1,iso2) %>% 
                       summarise(count = sum(count),.groups = "drop")
    isoform_coexist <- as.data.frame(isoform_coexist)
    if(nrow(isoform_coexist) == 1){
        return(isoform_coexist$iso1)
    }
    
    isoform_coexist <- as.matrix(long2wide(isoform_coexist,names_from = c("iso1","iso2"),
                                           value_from = "count",symmetric = F))
    isoform_coexist_filter <- cbind(diag(isoform_coexist),
                                    rowSums(isoform_coexist)-diag(isoform_coexist))
    isoform_coexist_filter <- isoform_coexist_filter[order(isoform_coexist_filter[,1],
                                                           decreasing = T),]

    isoform_coexist_filter <- names(which(isoform_coexist_filter[,1]/isoform_coexist_filter[,2] >= 1.5))
    
    return(isoform_coexist_filter)
}

cluster_isoform_correct <- function(splice_table,isoforms,
                                    polyA,start = "start",end = "end"){
    if(nrow(splice_table) != length(polyA)){
        stop("The size of isoforms and polyA don't match!")
    }
    cluster_size = nrow(splice_table)

    splice_table <- as.data.frame(splice_table)
    
    if(!is.na(isoforms[1])){
        splice_table_mid <- splice_table[,-which(colnames(splice_table) %in% c(start,end))]
        splice_table_mid <- as.data.frame(splice_table_mid)
        splice_table_mid_isoform = apply(splice_table_mid,1,function(x){
            return(paste(x,collapse = ""))
        })
        
        preserve = which(splice_table_mid_isoform %in% isoforms)
        if(length(preserve) == 0){
            return(NULL)
        }
        
        mode_isoform = table(splice_table_mid_isoform[preserve])
        mode_isoform = names(mode_isoform[which(mode_isoform == max(mode_isoform))])
        if(length(mode_isoform) > 1){
            loc= which(isoforms %in% mode_isoform)
            mode_isoform = mode_isoform[which(loc == min(loc))]
        }
        mode_preserve = which(splice_table_mid_isoform == mode_isoform)
    }
    else{
        mode_preserve = 1:nrow(splice_table)
        mode_isoform = NA
    }
    splice_table = splice_table[mode_preserve,]
        
    mode_start = table(splice_table[,start])
    mode_start = min(as.numeric(names(mode_start[which(mode_start == max(mode_start))])))
    
    mode_end = table(splice_table[splice_table[,start] == mode_start,end])
    mode_end = max(as.numeric(names(mode_end[which(mode_end == max(mode_end))])))        
    
    polyA = as.logical(polyA)
    polyA = mean(as.numeric(polyA[mode_preserve]))
    
    return(c(mode_start,mode_isoform,mode_end,cluster_size,polyA))
}

cell_clusters_isoform_correct <- function(splice_table,cluster,isoforms,polyA,
                                          start = "start",end = "end"){
    cluster_dic = unique(cluster)
    splice_table <- as.data.frame(splice_table)
    if(nrow(splice_table) != length(cluster)){
        stop("The size of isoforms and clusters don't match!")
    }
    if(nrow(splice_table) != length(polyA)){
        stop("The size of isoforms and polyA don't match!")
    }


    splice_table_adjust <- lapply(cluster_dic,function(x){
        key = which(cluster == x)
        sub_splice_table = splice_table[key,]
        #sub_UMI = UMI[key]
        sub_polyA = polyA[key]
        
        sub_splice_table_adjust = cluster_isoform_correct(sub_splice_table,isoforms,
                                                          sub_polyA,start,end)
        return(sub_splice_table_adjust)
    })
    splice_table_adjust <- do.call(rbind,splice_table_adjust)
    
    return(splice_table_adjust)
}


isoform_correct <- function(cells,splice_table,cluster,isoforms,polyA,
                            start = "start",end = "end"){
    cell_dic = unique(cells)
    splice_table <- as.data.frame(splice_table)
    if(nrow(splice_table) != length(cells)){
        stop("The size of isoforms and cells don't match!")
    }
    if(nrow(splice_table) != length(cluster)){
        stop("The size of isoforms and clusters don't match!")
    }
    if(nrow(splice_table) != length(polyA)){
        stop("The size of isoforms and polyA don't match!")
    }
    
    splice_table_adjust <- lapply(cell_dic,function(x){
        key = which(cells == x)
        sub_splice_table = splice_table[key,]
        sub_cluster = cluster[key]
        #sub_UMI = UMI[key]
        sub_polyA = polyA[key]
        
        sub_splice_table_adjust = cell_clusters_isoform_correct(sub_splice_table,sub_cluster,
                                                                isoforms,sub_polyA,start,end)
        if(is.null(sub_splice_table_adjust)){
            return(NULL)
        }
        sub_splice_table_adjust = cbind(x,sub_splice_table_adjust)
        return(sub_splice_table_adjust)
    })

    splice_table_adjust <- as.data.frame(do.call(rbind,splice_table_adjust))
    colnames(splice_table_adjust) <- c("cell","start","mid","end","size","polyA")
    splice_table_adjust$start = as.numeric(splice_table_adjust$start)
    splice_table_adjust$end = as.numeric(splice_table_adjust$end)
    splice_table_adjust$size = as.numeric(splice_table_adjust$size)
    splice_table_adjust$polyA = as.numeric(splice_table_adjust$polyA)
    return(splice_table_adjust)
}

site_recover <- function(start,mid,end,sites,sep = ",",split = "|"){
    if(is.na(mid)){
        splice_sites = paste(start,end,sep = sep)
    }
    else{
    if(nchar(mid) != length(sites)){
        stop("The size of splicing sites and binary indicator don't match!")
    }
    mid = as.logical(as.numeric(unlist(strsplit(mid,split = ""))))
    splice_sites = c(start,sites[mid],end)
    splice_sites = t(matrix(splice_sites,nrow = 2))
    splice_sites = apply(splice_sites,1,function(x) paste(x,collapse = sep))
    splice_sites = paste(splice_sites,collapse = split)
    }
    splice_sites = gsub(" ","",splice_sites,fixed = TRUE)
    return(splice_sites)
}