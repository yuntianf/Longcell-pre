### functions ###

##isoform length distance##
isoform2_dis <- function(a,b,split = "|",sep = ","){
    exons_a = unlist(strsplit(a,split = split,fixed = TRUE))
    exons_b = unlist(strsplit(b,split = split,fixed = TRUE))
    
    base_a = lapply(exons_a,function(x){
        x = unlist(strsplit(x,split = sep,fixed = TRUE))
        left = as.numeric(x[1])
        right = as.numeric(x[2])
        return(c(left:right))
    })
    base_a = unlist(base_a)
    
    base_b = lapply(exons_b,function(x){
        x = unlist(strsplit(x,split = sep,fixed = TRUE))
        left = as.numeric(x[1])
        right = as.numeric(x[2])
        return(c(left:right))
    })
    base_b = unlist(base_b)
    
    dis = length(base_a) + length(base_b) - 2*length(intersect(base_a,base_b))
    return(dis)
}

isoforms_dis <- function(isoforms,thresh = 10,split = "|",sep = ","){
    dis <- lapply(1:(length(isoforms)-1),function(i){
        sub_dis <- lapply((i+1):length(isoforms),function(j){
            temp = isoform2_dis(isoforms[i],isoforms[j],
                               split = split,sep = sep)
            if(temp <= thresh){
                return(c(i,j,temp))
            }
            else{
                return(NULL)
            }
        })
        sub_dis <- do.call(rbind,sub_dis)
    })
    dis <- as.data.frame(do.call(rbind,dis))  
    return(dis)
}

frNN_dis <- function(dis,n,eps = 10,if_direct = FALSE){
    dis <- as.matrix(dis)
    if(!if_direct){
        dis <- rbind(dis,dis[,c(2,1,3)])
    }
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

isoform_dis_cluster <- function(isoforms,thresh = 20,eps = 20,
                            split = "|",sep = ","){
    isoforms_count = table(isoforms)
    if(length(isoforms_count) == 1){
        return(rep(1,length(isoforms)))
    }
    
    dis = isoforms_dis(names(isoforms_count),
                       thresh = thresh,split = split,sep = sep)
    if(nrow(dis) == 0){
        cluster = 1:length(isoforms_count)
    }
    else{
        #dis[,3] = 1-dis[,3]
        dis_frNN = frNN_dis(dis,length(isoforms_count),eps = eps)
        cluster = dbscan(dis_frNN,minPts =1,weights = isoforms_count)$cluster
    }
    names(cluster) <- names(isoforms_count)
    cluster_expand = cluster[isoforms]
    return(cluster_expand)
}

##reads filter##
nbinom_kurt <- function(m,v){    
    kurt = 6*(v-m)/m^2+1/v
    return(kurt+3)
}

DF <- function(m,v,n){
    kurt = nbinom_kurt(m,v)
    df = 2*n/(kurt-(n-3)/(n-1)) 
    return(df)
}

var_pi <- function(m,v,n,alpha = 0.05){
    df = DF(m,v,n)
    pi = qchisq(c(alpha/2,1-alpha/2),df = df)*v/df
    return(pi)
}

sub_vec <- function(vec,loc){
    vec = sort(vec)
    cap = vec[loc:length(vec)]
    #if(loc > 1){
    #    interest = vec[1:(loc-1)]
    #    cap = cap + cap*sum(interest)/sum(cap)  
    #}  
    return(cap)
}

sub_vec_mv <- function(vec,loc){
    cap = sub_vec(vec,loc)
    return(c(length(cap),mean(cap),var(cap)))
}

mean2var <- function(m){
    return(2.566*m-2.566)
}

size_filter_pcr <- function(size,relation,alpha = 0.2){
    size = as.numeric(size)
    if(sum(is.na(size))){
        stop("size must be integer!")
    }
    if(length(size) == 1){
        return(1)
    }
    for(i in 1:length(size)){
        para = sub_vec_mv(size,i)
        sigma = relation(para[2])
        pi = var_pi(para[2],sigma,para[1],alpha = alpha)
        if(para[3] >= pi[1] & para[3] <= pi[2] | is.na(pi[1])){
            size_sort = sub_vec(size,i)
            size_sort = sort(size_sort,decreasing = TRUE)
            break
        }
    }
    return(length(size_sort))
}

size_filter_error <- function(size,ratio = 0.2){
    size = as.numeric(size)
    if(sum(is.na(size))){
        stop("size must be integer!")
    }
    if(length(size) == 1){
        return(1)
    }
    if(max(size) == 1){
        return(length(size))
    }
    size_sort = sort(size,decreasing = TRUE)
    size_summary <- lapply(1:sum(round(size_sort) > 1),function(i){
        left = sum(size_sort[1:i])*ratio/(1-ratio)
        right = sum(size_sort[(i+1):length(size_sort)])
        return(c(left,right))
    })
    size_summary <- as.data.frame(do.call(rbind,size_summary))
    colnames(size_summary) <- c("left","right")
    #print(size_summary)
    
    id = min(which(size_summary$left >= size_summary$right))    
    if(id == Inf){
        return(round(length(size) - max(size_summary$left)))
    }
    return(id)
}

size_filter <- function(size,relation,alpha = 0.05,ratio = 0.2){
    size = as.numeric(size)
    if(sum(is.na(size))){
        stop("size must be integer!")
    }
    if(length(size) == 1){
        return(0)
    }
    if(max(size) == 1){
        return(rep(0,length(size)))
    }
    
    for(i in 1:length(size)){
        para = sub_vec_mv(size,i)
        sigma = relation(para[2])
        pi = var_pi(para[2],sigma,para[1],alpha = alpha)
        #if(para[3] >= pi[1] & para[3] <= pi[2] | is.na(pi[1])){
        if(para[3] <= pi[2] | is.na(pi[1])){
            size_sort = sub_vec(size,i)
            size_sort = sort(size_sort,decreasing = TRUE)
            break
        }
    }
    #print(size_sort)
    size_summary <- lapply(1:sum(round(size_sort) > 1),function(i){
        left = sum(size_sort[1:i])*ratio/(1-ratio)
        right = sum(size_sort[(i+1):length(size_sort)])
        return(c(left,right))
    })
    size_summary <- as.data.frame(do.call(rbind,size_summary))
    colnames(size_summary) <- c("left","right")
    size_summary[is.na(size_summary)] = 0
    #print(size_summary)
    
    id = min(which(size_summary$left >= size_summary$right))
    vote = rep(0,length(size))
    #print(id)
    if(id == Inf){
        vote[size == 1] = max(size_summary$left)/sum(size == 1)
    }
    else{
        thresh = sort(size,decreasing = TRUE)[id]
        #print(thresh)
        vote[which(size < thresh)] = 1
        vote[which(size == thresh)] = (sum(size >= thresh) - id)/sum(size == thresh)
    }
    #print(vote)
    return(vote)
}

isoforms_size_filter <- function(isoform_table,relation,
                                 alpha = 0.05,ratio = 0.2,
                                 isoform = "isoform",mid = "mid",
                                 size = "size",
                                 polyA = "polyA"){
    reads_filter <- lapply(unique(isoform_table[,mid]),function(i){
        sub = isoform_table[isoform_table[,mid] == i,]
        cluster = isoform_dis_cluster(sub[,isoform],thresh = 10,eps = 10)
        weight = lapply(unique(cluster),function(x){
            id = which(cluster == x)
            sub_weight = 1 - size_filter(sub[id,size],relation = relation,
                                     alpha = alpha, ratio = ratio)
            return(cbind(id,sub_weight))
        })
        weight = as.data.frame(do.call(rbind,weight))
        colnames(weight) = c("id","weight")
        weight = weight[order(weight$id),]
        sub$weight = weight$weight
        
        sub = sub %>% 
              group_by(across(all_of(c(isoform)))) %>%
              summarise(size = sum(across(all_of(size))),
                        ### just for benchmark ###
                        cluster = n(),
                        ### just for benchmark ###
                        count = sum(weight),
                        polyA = mean(!! rlang::sym(polyA)),
                       .groups = "drop")
        #if(sum(sub$count == 0) > 0){
        #    sub = sub[-which(sub$count == 0),]
        #}
        return(sub)
    })
    reads_filter <- as.data.frame(do.call(rbind,reads_filter))
    return(reads_filter)
}

cells_isoforms_size_filter <- function(cell_isoform_table,
                                       relation,alpha = 0.05,ratio = 0.2,
                                       cell = "cell",isoform = "isoform",
                                       mid = "mid",size = "size",
                                       polyA = "polyA"){
    reads_filter = lapply(unique(cell_isoform_table[,cell]),function(i){
        #print(i)
        sub = cell_isoform_table[cell_isoform_table[,cell] == i,]
        filter = isoforms_size_filter(isoform_table = sub,relation = relation,
                                      alpha = alpha,ratio = ratio,
                                      mid = mid,isoform = isoform,
                                      size = size,polyA = polyA)
        filter = cbind(i,filter)
        colnames(filter)[1] = cell
        return(filter)
    })
    reads_filter = do.call(rbind,reads_filter)
    return(reads_filter)
}