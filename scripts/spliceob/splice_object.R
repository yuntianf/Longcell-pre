options(warn = -1)
library(tidyverse)
library(dplyr)
library(data.table)


splice_site <- function(isoform,split = "|",sep = ","){
    exons = unlist(strsplit(isoform,split = split,fixed = T))
    splice = unlist(strsplit(exons,split = sep))
    
    return(splice)
}

transcript_start <- function(isoforms,split = "|",sep = ","){
    start = sapply(isoforms,function(x){
        sites = splice_site(x,split,sep)
        return(as.numeric(sites[1]))
    })
    names(start) = NULL
    return(start)
}

transcript_end <- function(isoforms,split = "|",sep = ","){
    end = sapply(isoforms,function(x){
        sites = splice_site(x,split,sep)
        return(as.numeric(sites[length(sites)]))
    })
    names(end) = NULL
    return(end)
}

bin2exon <- function(bin,gene_bed,exon_id = NULL,
                     start = "start",end = "end",bias = 5){
    
    if(is.null(exon_id)){
        exon_id = 1:nrow(gene_bed)
    } 
    else{   
        if("N" %in% exon_id){
            warning("N means exons which don't correspond to annotation!")
        }
    }
    if(length(bin) != 2){
        stop("The input should be a bin!")
    }
    bin = as.numeric(bin)
    ### temp code: for a few bins with end smaller than start
    if(bin[2] < bin[1]){
        warning("The end of the bin is lower than the start")
        return("N")
    }
    
    cand = which(gene_bed[,start] <= bin[2] & gene_bed[,end] >= bin[1])
    if(length(cand) == 0){
        return("N")
    }
    
    exon_seq = c()
    i = 1
    while(i <= length(cand)){
        if(gene_bed[cand[i],start] > bin[1] + bias){
            exon_seq = c(exon_seq,"N")
            bin[1] = gene_bed[cand[i],start]
        }
        else{
            cover = min(gene_bed[cand[i],end],bin[2]) - bin[1] + 1
            exon_len = gene_bed[cand[i],end] - gene_bed[cand[i],start] + 1
            if(cover > 10 | cover/exon_len >= 0.5){
                exon_seq = c(exon_seq,exon_id[cand[i]])
            }
            bin[1] = gene_bed[cand[i],end] + 1
            i = i + 1
            if(bin[1] > bin[2]){
                break
            }
        }
    }
    if(bin[1] < bin[2]-bias){
        exon_seq = c(exon_seq,"N")
    }
    return(exon_seq)
}

splice2exon <- function(splice_sites,gene_bed,exon_id = NULL,
                        split = "|",sep = ",",start = "start",end = "end",
                        strand = "strand"){
    if(is.null(exon_id)){
        if(unique(gene_bed[,strand]) == "+"){
            exon_id = 1:nrow(gene_bed)
        }
        else{
            exon_id = nrow(gene_bed):1
        }
    }
    else if(length(gene_id) != nrow(gene_bed)){
        stop("The length of exon id doesn't correspond to the size of gene bed!")
    }
    
    bins = unlist(strsplit(splice_sites, split = split, fixed = TRUE))
    exons = lapply(bins,function(x){
        sites = as.numeric(unlist(strsplit(x, split = sep)))
        ### due to the bias in softclips extraction, remove 1 base to the end of exons   
        sites[2] = sites[2] - 1
        exon = bin2exon(bin = sites,gene_bed = gene_bed,exon_id = exon_id,
                         start = start,end = end) 
        return(exon)
    })
    exons = unlist(exons)
    isoform = paste(exons,collapse = split)
    return(isoform)
}

isoforms2exon <- function(isoforms,gene,gene_bed,exon_id = NULL,
                          split = "|",sep = ",",gene_col = "gene",
                          start = "start",end = "end"){
    gene_bed = gene_bed[gene_bed[,gene_col] == gene,]
    exon_seqs <- sapply(isoforms,function(x){
        seq = splice2exon(splice_sites = x,gene_bed = gene_bed,
                           exon_id = exon_id,split = split, sep = sep,
                           start = start,end = end)
        return(seq)
    })
    names(exon_seqs) <- NULL
    return(exon_seqs)
}

createExonList <- function(data,gene_bed,cells = NULL,genes = NULL,
                               cell_col = "cell",gene_col = "gene",
                               isoform_col = "isoform",split = "|",sep = ",",
                               bed_start = "start", bed_end = "end"){
    data = as.data.frame(data)
    if(is.null(cells)){
        cells = unique(data[,cell_col])
    }
    if(is.null(genes)){
        genes = unique(data[,gene_col])
    }
    
    data = data[data[,cell_col] %in% cells & data[,gene_col] %in% genes,]
    
    gene_exons = lapply(unique(data[,gene_col]),function(x){
        sub_data = data[data[,gene_col] == x,]
        start = transcript_start(sub_data[,isoform_col],split = split,sep=sep)
        end = transcript_end(sub_data[,isoform_col],split = split,sep=sep)
        exons = isoforms2exon(sub_data[,isoform_col],gene = x,
                              gene_bed = gene_bed,exon_id = NULL, 
                              split = split, sep = sep, gene_col = gene_col, 
                              start = bed_start, end = bed_end)
        
        rest = sub_data[,-which(colnames(sub_data) %in% 
                                c(cell_col,gene_col,isoform_col))]
        out = cbind(sub_data[,cell_col],x,start,exons,end,rest)
        out = as.data.frame(out)
        colnames(out)[1:5] <- c("cell","gene","start","exons","end")
        return(out)
    })
    gene_exons = do.call(rbind,gene_exons)
    gene_exons[,setdiff(colnames(gene_exons),c("cell","gene","exons"))] = sapply(gene_exons[,setdiff(colnames(gene_exons),c("cell","gene","exons"))],as.numeric)
    return(gene_exons)
}

setClass("Splice", 
         representation(cells = "data.frame", 
                        genes = "data.frame",
                        isoforms = "list"), 
         prototype(cells = NULL, genes = NULL,isoforms = NULL))

createSpliceObject <- function(data,gene_bed,cells = NULL,genes = NULL,
                               cell_col = "cell",gene_col = "gene",
                               isoform_col = "isoform",polyA = "polyA",
                               size = "size",split = "|",sep = ",",
                               bed_start = "start", bed_end = "end"){
    data = as.data.frame(data)
    if(is.null(cells)){
        cells = unique(data[,cell_col])
    }
    if(is.null(genes)){
        genes = unique(data[,gene_col])
    }
    
    data = data[data[,cell_col] %in% cells & data[,gene_col] %in% genes,]
    
    data$cell_id = as.numeric(as.factor(data[,cell_col]))
    
    #return(data)
    gene_exons = lapply(unique(data[,gene_col]),function(x){
        print(x)
        sub_data = data[data[,gene_col] == x,]
        start = transcript_start(sub_data[,isoform_col],split = split,sep=sep)
        end = transcript_end(sub_data[,isoform_col],split = split,sep=sep)
        exons = isoforms2exon(sub_data[,isoform_col],gene = x,
                              gene_bed = gene_bed,exon_id = NULL, 
                              split = split, sep = sep, gene_col = gene_col, 
                              start = bed_start, end = bed_end)
        
        out = cbind(sub_data$cell_id,start,exons,end,
                    sub_data$polyA,sub_data$size)
        out = as.data.frame(out)
        colnames(out) <- c("cell_id","start","exons","end","polyA","size")
        out[,setdiff(colnames(out),"exons")] = sapply(out[,setdiff(colnames(out),"exons")],as.numeric)
        
        return(out)
    })
    
    cell_summary <- data %>% 
                    group_by(across(all_of(cell_col))) %>%
                    summarise(gene_count = length(unique(!!sym(gene_col))),
                             umi_count = n())
    gene_summary <- data %>% 
                    group_by(across(all_of(gene_col))) %>%
                    summarise(cell_count = length(unique(!!sym(cell_col))),
                              total_exprs = n())    
    splice_ob <- new("Splice",
                     cells = cell_summary,
                     genes = gene_summary,
                     isoforms = gene_exons)
    
    return(splice_ob)
}

exonList2spliceOb <- function(data,cells = NULL,genes = NULL,
                              cell_col = "cell",gene_col = "gene",
                              count_col = "count"){
    data = as.data.frame(data)
    if(is.null(cells)){
        cells = unique(data[,cell_col])
    }
    if(is.null(genes)){
        genes = unique(data[,gene_col])
    }

    data = data[data[,cell_col] %in% cells & data[,gene_col] %in% genes,]
        
    data = cbind(as.numeric(as.factor(data[,cell_col])),data)
    colnames(data)[1] = "cell_id"
    
    gene_exons = split(data[,setdiff(colnames(data),c(cell_col,gene_col))],
                      data[,gene_col])
    
    cell_summary <- data %>% 
                    group_by(across(all_of(cell_col))) %>%
                    summarise(gene_count = length(unique(!!sym(gene_col))),
                             umi_count = sum(!!sym(count_col)))
    gene_summary <- data %>% 
                    group_by(across(all_of(gene_col))) %>%
                    summarise(cell_count = length(unique(!!sym(cell_col))),
                              total_exprs = sum(!!sym(count_col)))    
    splice_ob <- new("Splice",
                     cells = cell_summary,
                     genes = gene_summary,
                     isoforms = gene_exons)    
    return(splice_ob)
}

exonTable <- function(exons,count,polyA, split = "|"){
    exon_uniq = unique(exons)
    
    exon_vec_dic <- lapply(exon_uniq,function(x){
        exon_seq = as.numeric(unlist(strsplit(x,split = split,fixed = TRUE)))
        na_id = which(is.na(exon_seq))
        #if(sum(na_id > 1 & na_id < length(exon_seq)) > 0)
        if(sum(is.na(exon_seq)) > 0){
            return(NULL)
        }
        exon_seq = na.omit(exon_seq)
        return(exon_seq)
    })
    names(exon_vec_dic) <- exon_uniq
    
    exon_dic = sort(unique(unlist(exon_vec_dic)))
    if(is.null(exon_dic)){
        return(NULL)
    }
    
    exon_vec <- lapply(1:length(exons),function(i){
        vec = exon_vec_dic[[exons[i]]]
        exon_count = count[[i]]
        
        out = rep(0,length(exon_dic))
        names(out) <- exon_dic
        if(length(vec) == 0){
            out[1:length(out)] <- NA 
        }
        else{
            out[as.character(vec)] = exon_count
            out[exon_dic < min(vec)] = NA
            if(polyA[i] < 0.5){
                out[exon_dic > max(vec)] = NA
            }            
        }
        return(out)
    })
    
    exon_vec <- as.data.frame(do.call(rbind,exon_vec))
    return(exon_vec)
}

annoExonTable <- function(gene_bed,gtf,gene,split = "|"){
  isoforms = gtf_bed_corres(gtf = gtf,bed = gene_bed,gene = gene)
  isoform_table = exonTable(exons = isoforms$exons,count = rep(1,nrow(isoforms)),
                            polyA = rep(TRUE,nrow(isoforms)),split = split)
  isoform_table[is.na(isoform_table)] = 0
  return(isoform_table)
}

extractGenes <- function(spliceOb){
    return(spliceOb@genes$gene)
}
extractCells <- function(spliceOb){
    return(spliceOb@cells$cell)
}
extractIsoform <- function(spliceOb,gene){    
    isoform = spliceOb@isoforms[[gene]]
    return(isoform)
}

# identify if two exons can be merged(co-existing or mutually exclusive)
ifMerge <- function(a,b){
    if(length(a) != length(b)){
        stop("Length of two input vectors should be the same!")
    }
    
    exon_corres = sum(na.omit(xor(a,b)))
    if(exon_corres == 0 | exon_corres == length(na.omit(a))){
      na_corres = xor(is.na(a),is.na(b))
      if(sum(na_corres) == 0){
        return(TRUE)
      }     
      else{
        return(NA)
      }
    }
    else{
      return(FALSE)
    }
}

# build a relation table for exons to see if they could be merged
exonTableMergeFlag <- function(exon_table){
  exons = colnames(exon_table)
  
  merge_table <- lapply(1:ncol(exon_table),function(i){
    merge_vec <- lapply(i:ncol(exon_table),function(j){
      return(c(i,j,ifMerge(exon_table[,i],exon_table[,j])))
    })
    merge_vec = do.call(rbind,merge_vec)
    return(merge_vec)
  })

  merge_table = as.data.frame(do.call(rbind,merge_table))
  colnames(merge_table) = c("exon1","exon2","if_merge")
  
  merge_matrix = pivot_wider(merge_table,names_from = "exon2",values_from = "if_merge")
  merge_matrix = as.data.frame(merge_matrix)
  
  rownames(merge_matrix) = exons[merge_matrix$exon1]
  merge_matrix = merge_matrix[,-1]
  colnames(merge_matrix) = rownames(merge_matrix)
  
  merge_matrix[lower.tri(merge_matrix)] <- t(merge_matrix)[lower.tri(merge_matrix)]
  diag(merge_matrix) = 1
  return(merge_matrix)
}

# combine the merge table from data(a) and annotation(b) 
dataAnnoCombine <- function(a,b){
  exons = colnames(a)
  if(length(setdiff(exons,colnames(b))) > 0){
    stop("The exons in the data and annotation don't correspond")
  }
  
  b = b[exons,exons]
  if(length(a) != length(b)){
    stop("The size of the matrix should be the same")
  }
  a = as.matrix(a)
  b = as.matrix(b)
  
  merge_flag = sapply(1:length(a),function(i){
    if(!is.na(a[i])){
      return(a[i])
    }
    else if(!is.na(b[i])){
      return(b[i])
    }
    else{
      stop("The value from annotation should be definite!")
    }
  })
  
  merge_flag = matrix(merge_flag,nrow(a),ncol(a))
  rownames(merge_flag) = rownames(a)
  colnames(merge_flag) = colnames(a)
  return(merge_flag)
}

# extract exon list to be merged from exon merge table
moduleExtract <- function(data){
  if(nrow(data) != ncol(data)){
    stop("The similarity matrix should be symmetric!")
  }
  if(is.null(rownames(data)) | is.null(colnames(data))){
    warning("There is no annotation of each row and col of the data, will use the id")
    rownames(data) = 1:nrow(data)
    colnames(data) = 1:ncol(data)
  }
  
  data = data[order(rowSums(data),decreasing = TRUE),
              order(rowSums(data),decreasing = TRUE)]
  exons = rownames(data)
  
  merge_list = list()
  while(length(data) > 1){
    if(sum(data) == nrow(data)){
      merge_list = c(merge_list,as.list(rownames(data)))
      break
    }
    id = which(data[1,] == 1)
    subdata = data[id,id]
    while(sum(subdata == 0) > 0){
      filter = which(rowSums(subdata) == min(rowSums(subdata)))[1]
      id = id[-filter]
      subdata = subdata[-filter,-filter]
    }
    merge_list = c(merge_list,list(rownames(data)[id]))
    data = data[-id,-id]
  }
  merge_list <- c(merge_list,as.list(setdiff(exons,unlist(merge_list))))
  return(merge_list)
}

exonMatrixMerge <- function(dataExonMergeMatrix,annoExonMergeMatrix = NULL){
  if(is.null(annoExonMergeMatrix)){
    annoExonMergeMatrix = matrix(0,nrow(dataExonMergeMatrix),ncol(dataExonMergeMatrix))
    annoExonMergeMatrix = as.data.frame(annoExonMergeMatrix)
    rownames(annoExonMergeMatrix) = rownames(dataExonMergeMatrix)
    colnames(annoExonMergeMatrix) = colnames(dataExonMergeMatrix)
  }
  
  exons = rownames(dataExonMergeMatrix)
  annoExonMergeMatrix = annoExonMergeMatrix[exons,exons]
  if(sum(is.na(annoExonMergeMatrix)) > 0){
    stop("The exons in the data and annotation don't correspond!")
  }
  
  if(nrow(dataExonMergeMatrix) != nrow(annoExonMergeMatrix) | 
     ncol(dataExonMergeMatrix) != ncol(annoExonMergeMatrix)){
    stop("The number of exons in data and annotations should be the same!")
  }
  
  exon_merge_matrix = dataAnnoCombine(dataExonMergeMatrix,annoExonMergeMatrix)
  
  merge_list = moduleExtract(exon_merge_matrix)
  return(merge_list)
}

# merge two exons
exonExonMerge <- function(a,b){
  if(length(a) != length(b)){
    stop("The size of the matrix should be the same")
  }
  
  if(sum(xor(a,b),na.rm = TRUE) > 0){
    b = !b
  }
  
  merge_flag = sapply(1:length(a),function(i){
    if(is.na(a[i])){
      return(b[i])
    }
    if(is.na(b[i])){
      return(a[i])
    }
    if(!is.na(a[i]) & !is.na(b[i])){
      if(a[i] == b[i]){
        return(a[i])
      }
      else{
        stop("There exist collision between exons to be merged!")
      }
    }
  })
  
  return(merge_flag)
}

# merge multi exons
multiExonsMerge <- function(data){
  if(is.null(ncol(data)) || ncol(data) == 1){
    return(data)
  }
  
  merged_exon = rep(NA,nrow(data))
  for(i in 1:ncol(data)){
    merged_exon = exonExonMerge(merged_exon,data[,i])
  }
  return(merged_exon)
}

# merge exons in an exon table according to to be merged list
exonTableMergeList <- function(exon_table,merge_list,sep = "|"){
  gene_count = apply(exon_table,1,function(x){
    return(max(x,na.rm = TRUE))
  })
  exon_table = apply(exon_table,2,as.logical)
  merged_exon_table <- lapply(merge_list,function(x){
    temp = multiExonsMerge(exon_table[,x])
  })

  merged_exon_table = as.data.frame(do.call(cbind,merged_exon_table))
  
  merged_exons = sapply(merge_list,function(x){
    exons = paste(x,collapse = sep)
  })
  
  colnames(merged_exon_table) = merged_exons
  merged_exon_table$gene_count = gene_count
  return(merged_exon_table)
}

exonTableMerge <- function(exon_table,
                           gene_bed = NULL,gtf = NULL,
                           gene = NULL,sep = "|"){
  if(is.null(ncol(exon_table)) || ncol(exon_table) == 1){
    exon_table = as.data.frame(exon_table)
    gene_count = apply(exon_table,1,max)
    gene_count[is.na(gene_count)] = 0
    merged_table = as.data.frame(cbind(exon_table,gene_count))
    colnames(merged_table) = c("exon","gene_count")
    return(merged_table)
  }
  exon_merge_flag = exonTableMergeFlag(exon_table)
  
  if(!is.null(gene_bed) & !is.null(gtf) & !is.null(gene)){
    anno_exon_table = annoExonTable(gene_bed = gene_bed,gtf = gtf,
                                    gene = gene,split = sep)
    anno_merge_flag = exonTableMergeFlag(anno_exon_table)
  }
  else{
    anno_merge_flag = NULL
  }
  merge_list = exonMatrixMerge(exon_merge_flag,anno_merge_flag)
  merged_table = exonTableMergeList(exon_table,merge_list)
  
  return(merged_table)
}

exon_len_filter <- function(gene,gene_bed,thresh = 10,exons = NULL,
                            exon_id = NULL,gene_id = "gene",
                            len = "width",strand = "strand",sep = "|"){
  gene_bed = gene_bed[gene_bed[,gene_id] == gene,]
  
  if(is.null(exon_id)){
    if(unique(gene_bed[,strand]) == "+"){
      exon_id = 1:nrow(gene_bed)
    }
    else{
      exon_id = nrow(gene_bed):1
    }
  }
  
  if(is.null(exons)){
    exons = exon_id
  }
  
  filtered = sapply(exons,function(x){
    e = unlist(strsplit(x,split = sep,fixed = TRUE))
    exons_len = sum(gene_bed[exon_id %in% e,len])
    if(exons_len >= thresh){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  })
  return(exons[filtered])
}

exon_search <- function(exon,exons_group,sep = "|"){
  exons_group = strsplit(exons_group,split = sep,fixed = TRUE)
  id = sapply(exons_group,function(x){
    flag = exon %in% x
    return(flag)
  })
  if(sum(id) == 0){
    stop("This exon doesn't exist!")
  }
  else if(sum(id) > 1){
    warning("The exon is duplicated!")
  }
  return(which(id)[1])
}

extractExonTable <- function(spliceOb,gene,gene_bed,gtf = NULL,
                             cells = "all",exons = "all",
                             exon_len_thresh = 10,sep = "|"){
  isoform = extractIsoform(spliceOb,gene)
  if(cells != "all"){
    isoform = isoform[isoform$cell_id %in% cells,]
  } 
  
  if(length(isoform) == 0 || nrow(isoform) == 0){
    return(NULL)
  }
  
  isoform_table = exonTable(isoform$exons,isoform$count,isoform$polyA)
  isoform_table = cbind(isoform$cell_id,isoform_table)
  colnames(isoform_table)[1] = "cell_id"
  
  row_filter = rowSums(isoform_table,na.rm = TRUE)-isoform_table[,1]
  col_filter = colSums(isoform_table,na.rm = TRUE)
  isoform_table = isoform_table[row_filter > 0,col_filter > 0]

  if(length(isoform_table) == 0 || nrow(isoform_table) == 0){
    return(NULL)
  }

  #return(isoform_table)
  all_exons = setdiff(colnames(isoform_table),"cell_id")
  isoform_table = cbind(isoform_table[,"cell_id"],
                     exonTableMerge(isoform_table[,all_exons],
                                    gene_bed = gene_bed,gtf = gtf,
                                    gene = gene,sep = sep))
  colnames(isoform_table)[1] = "cell_id"
  if(length(all_exons) == 1){
    colnames(isoform_table)[2] = all_exons
  }
  
  if(exons == "all"){
    exons = setdiff(colnames(isoform_table),c("cell_id","gene_count"))
    exons = exon_len_filter(gene = gene,gene_bed = gene_bed,
                            thresh = exon_len_thresh,exons = exons)
  }   
  else{
    if(length(setdiff(exons,colnames(isoform_table))) > 0){
      warning("Some designated exons don't exsit or be merged, will use the intersection!")
    }
    exons_id = sapply(exons,function(i){
      exon_search(exon = i,exons_group = colnames(isoform_table))
    })
    exons = colnames(isoform_table)[unique(unlist(exons_id))]
  }  
  
  isoform_table = isoform_table[,colnames(isoform_table) %in% c("cell_id",exons,"gene_count")]
  return(isoform_table)
}

spliceOb2dataFrame <- function(spliceOb){
    len = sapply(spliceOb@isoforms,nrow)
    genes = names(spliceOb@isoforms)
    isoforms = as.data.frame(do.call(rbind,spliceOb@isoforms))
    isoforms$gene = rep(genes,len)
    
    cells = spliceOb@cells
    isoforms$cell = cells[isoforms$cell_id]
    
    isoforms = isoforms[,c("cell","gene",
                         setdiff(colnames(isoforms,c("cell","gene","cell_id"))))]
    return(isoforms)
}

subsetSpliceOb <- function(spliceOb,cells = NULL,genes = NULL,count_col = "count"){
    cell_attr = spliceOb@cells
    gene_attr = spliceOb@genes
    ### temp for current version, should be removed
    cell_attr$cell_id = 1:nrow(cell_attr)    
    
    if(is.null(cells) & is.null(genes)){
        warning("No cells and genes are specified for subseting, will re turn the original spliceOb!")
        return(spliceOb)
    }
    else if(is.null(cells)){
        cells = cell_attr$cell
    }
    else if(is.null(genes)){
        genes = gene_attr$gene
    }
    cells = cells[cells %in% cell_attr$cell]
    genes = genes[genes %in% gene_attr$gene]
    
    cell_attr = cell_attr[cell_attr$cell %in% cells,]
    gene_attr = gene_attr[gene_attr$gene %in% genes,]
    
    isoforms = spliceOb@isoforms[genes]
    iso_len = sapply(isoforms,nrow)
    isoforms = as.data.frame(do.call(rbind,isoforms))
    isoforms$gene = rep(1:length(genes),iso_len)
    
    ### temp for current version, should be removed
    isoforms$count = 1
    isoforms = isoforms[isoforms$cell_id %in% cell_attr$cell_id,]

    cell_summary <- isoforms %>% 
                    group_by(cell_id) %>%
                    summarise(gene_count = length(unique(gene)),
                             umi_count = sum(!!sym(count_col)))
    
    gene_summary <- isoforms %>% 
                    group_by(gene) %>%
                    summarise(cell_count = length(unique(cell_id)),
                              total_exprs = sum(!!sym(count_col)))
    gene_summary$gene <- genes[gene_summary$gene]
    
    cell_attr = left_join(cell_attr[,-which(colnames(cell_attr) %in% c("gene_count","umi_count"))],
                         cell_summary,by = "cell_id")
    gene_attr = left_join(gene_attr[,-which(colnames(gene_attr) %in% c("cell_count","total_exprs"))],
                         gene_summary,by = "gene")
    
    if(sum(is.na(cell_attr$gene_count)) > 0){
        cell_attr <- cell_attr[-which(is.na(cell_attr$gene_count)),]
    }
    if(sum(is.na(gene_attr$cell_count)) > 0){
        gene_attr <- gene_attr[-which(is.na(gene_attr$cell_count)),]
    }    
    
    if(length(unique(isoforms$cell_id))!= length(cell_attr$cell_id)){
        warning("cell count don't match!")
    }
    cell_attr$cell_id = as.numeric(as.factor(cell_attr$cell_id))
    isoforms$cell_id = as.numeric(as.factor(isoforms$cell_id))
    
    rownames(isoforms) <- NULL
    genes = genes[unique(isoforms$gene)]
    isoforms = split(isoforms[,-which(colnames(isoforms) == "gene")],isoforms$gene)
    names(isoforms) = genes
    
    sub_splice_ob <- new("Splice",
                     cells = cell_attr,
                     genes = gene_attr,
                     isoforms = isoforms)   
    return(sub_splice_ob)
}

cellGeneExonCount <- function(exon_table){ 
  exon_table <- apply(exon_table,2,as.numeric)
  
  gene_table = apply(exon_table,1,function(x){
    x[2:(ncol(exon_table)-1)][!is.na(x[2:(ncol(exon_table)-1)])] = x[ncol(exon_table)]
    return(x[1:(ncol(exon_table)-1)])
  })
  gene_table = as.data.frame(t(gene_table))

  exons_table = apply(exon_table,1,function(x){
    x[2:(ncol(exon_table)-1)][which(x[2:(ncol(exon_table)-1)] > 0)] = x[ncol(exon_table)]
    return(x[1:(ncol(exon_table)-1)])
  })
  exons_table = as.data.frame(t(exons_table))
  
  cell_gene_count = gene_table %>% 
    group_by(cell_id) %>% 
    summarise_all(funs(sum(na.omit(.))))
  cell_exon_count = exons_table %>% 
    group_by(cell_id) %>% 
    summarise_all(funs(sum(na.omit(.))))
  
  return(list(cell_gene_count,cell_exon_count))
}

gene_exons_table <- function(splice_ob,gene,gene_bed,gtf = NULL,cells = "all",
                             exons = "all",exon_len_thresh = 10){
    exon_table = extractExonTable(splice_ob,gene,gene_bed,gtf = gtf,cells = cells,
                                  exons = exons,exon_len_thresh = exon_len_thresh)
    if(is.null(exon_table) || nrow(exon_table) == 0){
        return(NULL)
    }
    
    count_list = cellGeneExonCount(exon_table)
    return(count_list)
}

gene_exons_psi <- function(splice_ob,gene,gene_bed,gtf = NULL,cells = "all",
                           exons = "all",exon_len_thresh = 10,exprs_thresh = 10){
    count_list = gene_exons_table(splice_ob,gene,gene_bed,gtf = gtf,cells = cells,
                                  exons = exons,exon_len_thresh = exon_len_thresh)
    if(is.null(count_list)){
        return(NULL)
    }
    gene_count_table = count_list[[1]]
    exon_count_table = count_list[[2]]  
    
    exons = setdiff(colnames(exon_count_table),"cell_id")

    exon_psi_table = exon_count_table[,exons]/gene_count_table[,exons]
    exon_psi_table[gene_count_table[,exons] < exprs_thresh] = NA
    exon_psi_table = cbind(exon_count_table[,"cell_id"],exon_psi_table)
    colnames(exon_psi_table)[1] = "cell_id"
    
    row_filter = apply(exon_psi_table,1,function(x){
        return(sum(is.na(x)) < length(x) - 1)
    })
    col_filter = apply(exon_psi_table,2,function(x){
        return(sum(is.na(x)) < length(x))
    })
    exon_psi_table = exon_psi_table[row_filter,col_filter]
    
    return(exon_psi_table)
}


phi <- function(x){
  return(var(x)/(mean(x)*(1-mean(x))))
}

phi_conf <- function(x,iters = 100){
  conf <- sapply(1:iters,function(i){
    sample_x <- sample(x,length(x),replace = TRUE)
    sample_phi <- phi(sample_x)
    return(sample_phi)
  })
  out <- na.omit(conf)
  return(c(mean(out),quantile(out,c(0.025,0.975),na.rm = TRUE)))
}

gene_exons_phi <- function(splice_ob,gene,gene_bed,gtf = NULL,cells = "all",
                           exons = "all",exon_len_thresh = 10,
                           cell_count_thresh = 30,exprs_thresh = 10,
                           iters = 100){
    exon_psi_table = gene_exons_psi(splice_ob,gene,gene_bed,gtf = NULL,cells = cells,
                                    exons = exons,exon_len_thresh = exon_len_thresh,
                                    exprs_thresh = exprs_thresh)
    if(is.null(exon_psi_table) || length(exon_psi_table) == 0){
        return(NULL)
    }
    exons = colnames(exon_psi_table)[which(colSums(!is.na(exon_psi_table)) >= cell_count_thresh)]
    exons = setdiff(exons,"cell_id")
    
    if(length(exons) == 0){
        return(NULL)
    }
    else{
        exons_phi <- lapply(exons,function(i){
            psi = na.omit(unlist(exon_psi_table[,i]))
            conf = phi_conf(psi,iters = iters)
            return(c(i,mean(psi),conf,length(psi)))
        })
        exons_phi <- as.data.frame(do.call(rbind,exons_phi))
        exons_phi[,2:ncol(exons_phi)] <- apply(exons_phi[,2:ncol(exons_phi)],
                                               2,as.numeric)
        colnames(exons_phi) <- c("exon","mean_psi","phi","phi_lwr","phi_upr","count")
        exons_phi <- as.data.frame(exons_phi)
        exons_phi$gene = gene
        return(as.data.frame(exons_phi))
    }
}