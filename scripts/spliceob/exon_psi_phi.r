options(warn = -1)
library(tidyverse)
library(Matrix)
library(dplyr)
library(data.table)
library(stats)

### correpsond exons in canonical annotations and home-made bed files
exon_correspond <- function(exon_parts,exon_ann,start = "start",end = "end",flank = 10){
    exon_ann[,start] <- as.numeric(exon_ann[,start])
    exon_ann[,end] <- as.numeric(exon_ann[,end])
    
    if(nrow(exon_parts) > 1){
        exon_merge <- lapply(1:nrow(exon_parts),function(i){
            return(which(exon_ann[,start] <= (as.numeric(exon_parts[i,start])+flank) & 
                         exon_ann[,end] >= (as.numeric(exon_parts[i,end])-flank)))
        })
        exon_merge <- as.data.frame(cbind(rep(exon_parts$exon_id, 
            sapply(exon_merge, length)), unlist(exon_merge)))
    }
    else{
        exon_merge = t(as.data.frame(c(1,1)))
    }
    exon_merge <- as.data.frame(exon_merge)
    colnames(exon_merge) <- c("part","exon")
    rownames(exon_merge) <- NULL
    exon_merge <- exon_merge[order(exon_merge$exon,exon_merge$part),]
   
    return(exon_merge)
}

### split reads into exons
reads_split <- function(reads, sep = "|", blank = ""){
    if (blank != "") {
        reads <- gsub(blank, "", reads)
    }
    return(as.numeric(unlist(strsplit(reads, sep, fixed = T))))
}

### annotate canonical exons in home-made exons
transcript_part <- function(gene_gtf,gene_bed,gene){
    gene_gtf = gene_gtf[gene_gtf$gene_id == gene,]
    gene_bed = gene_bed[gene_bed$gene == gene,]

    gene_gtf_exon <- unique(gene_gtf[,c("start","end","exon_id")])
    gene_gtf_exon <- gene_gtf_exon[order(gene_gtf_exon$start,gene_gtf_exon$end),]
    gene_gtf_exon$exon_num <- 1:nrow(gene_gtf_exon)
    gene_gtf <- left_join(gene_gtf,gene_gtf_exon[,c(3,4)],by = "exon_id")

    gene_transcript <- lapply(unique(gene_gtf$transname),
                           function(x) gene_gtf[gene_gtf$transname == x,"exon_num"])
    names(gene_transcript) <- unique(gene_gtf$transname)
    gene_corres <- exon_correspond(gene_bed,gene_gtf_exon)
                         
    gene_transcript_part <- lapply(gene_transcript,function(x) 
        as.character(sort(as.numeric(gene_corres[gene_corres$exon %in% x,"part"]))))
    
    gene_transcript_part <- gene_transcript_part[order(names(gene_transcript_part))]    
    return(gene_transcript_part)
}

### define if a read is truncated
trunc_status <- function(cell_exon,genes,gene_bed,gtf){
    update_cell_exon <- lapply(genes,function(x){
        print(x)
        sub_cell_exon <- cell_exon[cell_exon$gene == x,]
        bed = gene_bed[gene_bed$gene == x,]
        if(nrow(bed) == 1){
            return(NULL)
        }
        bed$exon_id = as.numeric(bed$exon_id)
        strand = ifelse(bed$exon_id[1] > bed$exon_id[2],"reverse","forward")
        
        transcript = transcript_part(gtf,gene_bed,x)
        if(strand == "forward"){
            transcript = sapply(transcript,function(y) paste(c(y,""),collapse = "|"))
        }else{
            transcript = sapply(transcript,function(y) paste(c(rev(y),""),collapse = "|"))                     
        }
                                
        trunc <- t(apply(sub_cell_exon,1,function(y){
            trunc_state = 0
            exons = reads_split(y["exon"])
            coverage = reads_split(y["coverage"])
            if(max(exons) != max(bed$exon_id)){
                if(coverage[which(exons == max(exons))] < 0.9){
                    trunc_state = trunc_state+3
                }
            }
            if(min(exons) != min(bed$exon_id)){
                if(coverage[which(exons == min(exons))] < 0.9){
                    trunc_state = trunc_state+5
                }
            }
            
            exist_ann <- y["exon"] %in% transcript
            return(c(trunc_state,exist_ann))
        }))
        colnames(trunc) <- c("trunc","ann")
        sub_cell_exon <- cbind(sub_cell_exon,trunc)
        return(sub_cell_exon)
    })
    update_cell_exon[vapply(update_cell_exon,is.null,logical(1L))] <- NULL
    update_cell_exon <- as.data.frame(do.call(rbind,update_cell_exon))   
    return(update_cell_exon)
}

### calculate cencored psi for a gene in a single cell
cell_exon_psi_cencor <- function(cell_exons,bed,sep = "|", 
                          exon = "exon",trunc = "trunc",polyA = "polyA",thresh = 10){
    
    dic = bed$exon_id
    exon_matrix <- t(apply(cell_exons,1,function(x){
        exons = reads_split(x[exon])
        in_vec <- rep(0,length(dic))
        exist_vec <- rep(0,length(dic))
        names(in_vec) <- dic
        names(exist_vec) <- dic
        
        exist_vec[as.character(exons)] <- 1
        
        if(x[trunc] == 8){
            in_vec[as.character(min(exons):max(exons))] <- 1
        }
        else if(x[trunc] == 3){
            in_vec[as.character(min(as.numeric(dic)):max(exons))] <- 1
        }
        else if(x[trunc] == 5){
            if(x[polyA] > 0.5){
                in_vec[as.character(min(exons):max(as.numeric(dic)))] <- 1
            }
            else{
                in_vec[as.character(min(exons):max(exons))] <- 1
            }
        }
        return(c(exist_vec,in_vec))
    }))
    
    exist_matrix <- exon_matrix[,1:(ncol(exon_matrix)/2)]
    in_matrix <- exon_matrix[,(ncol(exon_matrix)/2 + 1):ncol(exon_matrix)]
    
    exist_in_matrix <- exist_matrix & in_matrix
    
    if(is.null(dim(in_matrix))){
        in_matrix <- t(as.matrix(in_matrix))
        exist_in_matrix <- t(as.matrix(exist_in_matrix))
    }
                             
    cencor_psi <- colSums(exist_in_matrix)/colSums(in_matrix)
    cencor_psi[colSums(in_matrix) < thresh] <- NA
    
    return(cencor_psi)
}

### calculate exon psi for all cells
exon_psi_cencor <- function(cell_exons,bed,sep = "|",cell = "cell", 
                          exon = "exon",trunc = "trunc",polyA = "polyA", thresh = 10){
    cells <- unique(cell_exons[,cell])
    
    out <- lapply(cells,function(x){
        cell_exon <- cell_exons[cell_exons[,cell] == x,]
        cell_exon_psi = cell_exon_psi_cencor(cell_exon,bed,
                                             sep = sep, exon = exon,trunc = trunc,
                                             polyA = polyA,thresh = 5)
    })
    
    out <- as.data.frame(do.call(rbind,out))
    if(is.null(out) | nrow(out) == 0){
        return(NULL)
    }
    
    rownames(out) <- cells
    out <- out[,as.character(sort(as.numeric(colnames(out))))]
    out <- out[rowSums(is.na(out)) != ncol(out), ]
    out <- out[,colSums(is.na(out)) != nrow(out)]
    return(as.data.frame(out))
}

### calculate exon phi a gene in single cell 
cell_exon_phi_cencor <- function(exon_psi,cell_thresh = 20){
    if(nrow(exon_psi) < cell_thresh){
        return(NULL)
    }
    else{
        exons <- colnames(exon_psi)[colSums(!is.na(exon_psi)) >= cell_thresh]
        if(length(exons) == 0){
            return(NULL)
        }
        exon_psi <- exon_psi[,colSums(!is.na(exon_psi)) >= cell_thresh]
        exon_psi <- as.data.frame(exon_psi)
        colnames(exon_psi) <- exons
        
        out <- t(apply(exon_psi,2,function(x){
            x <- na.omit(x)
            phi <- sapply(1:100,function(i){
                sample_x <- sample(x,length(x),replace = T)
                sample_mean = mean(sample_x)
                sample_phi <- var(sample_x)/(sample_mean*(1-sample_mean))
                return(sample_phi)
            })
            phi <- na.omit(phi)
            return(c(mean(x),mean(phi),quantile(phi,c(0.025,0.975),na.rm = T)))
        }))
        out <- as.data.frame(out)
        colnames(out) <- c("psi","phi","phi_left","phi_right")
        out$exon <- colnames(exon_psi)
        return(out)
    }
}

### calculate exon phi for all cells for a gene
exon_phi_cencor <- function(exon_psi,cell_thresh = 20){
    out <- lapply(exon_psi,function(i){
        cell_exon_phi_cencor(i,cell_thresh)
    })
    gene <- rep(names(exon_psi),sapply(out,function(x) ifelse(is.null(nrow(x)),0,nrow(x))))
    out[vapply(out,is.null,logical(1L))] <- NULL
    out <- as.data.frame(do.call(rbind,out))
    out$gene <- gene
    rownames(out) <- NULL
    out <- out[,c("gene","exon","psi","phi","phi_left","phi_right")]
    return(out)
}