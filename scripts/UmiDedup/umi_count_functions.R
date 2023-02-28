### functions ###
isoform_correct_filter <- function(gene_cells_cluster,edit,
                            split = "|",sep = ",",
                            splice_site_thresh = 10,
                            relation = mean2var,alpha = 0.05){

  gene_isoform = splice_site_table(gene_cells_cluster[,"isoform"],
                                   split,sep,splice_site_thresh)

  gene_isoform_filter = vapply(gene_isoform,is.null,logical(1L))
  if(sum(gene_isoform_filter) > 0){
    gene_cells_cluster = gene_cells_cluster[-which(gene_isoform_filter),]
  }

  gene_isoform <- as.data.frame(do.call(rbind,gene_isoform))

  if(length(gene_isoform) == 0){
    return(NULL)
  }
  colnames(gene_isoform)[1] <- "start"
  colnames(gene_isoform)[ncol(gene_isoform)] <- "end"

  if(ncol(gene_isoform) > 2){
    filter = as.data.frame(gene_isoform[,2:(ncol(gene_isoform)-1)] == "1")
    filter = colSums(filter) > 0
    gene_isoform = gene_isoform[,c(TRUE,filter,TRUE)]
  }

  if(ncol(gene_isoform) > 2){
    splice_sites = colnames(gene_isoform)[2:(ncol(gene_isoform)-1)]
    isoform_filtered = cells_mid_filter(gene_cells_cluster$cell,
                                        gene_isoform[,2:(ncol(gene_isoform)-1)],
                                        gene_cells_cluster$cluster)
  }
  else{
    splice_sites = NA
    isoform_filtered = NA
  }

  gene_isoform = isoform_correct(gene_cells_cluster$cell,gene_isoform,
                                 gene_cells_cluster$cluster,
                                 isoform_filtered,gene_cells_cluster$polyA)

  gene_isoform$isoform <- apply(gene_isoform,1,function(x){
    site_recover(x["start"],x["mid"],x["end"],splice_sites)
  })

  gene_isoform <- gene_isoform[,c("cell","isoform","mid","size","polyA")]
  gene_isoform$mid[is.na(gene_isoform$mid)] = "null"

  ratio = sum(edit > 2)/length(edit)
  gene_isoform = cells_isoforms_size_filter(cell_isoform_table = gene_isoform,
                                            relation = relation, alpha = alpha,
                                            ratio = ratio)
  return(gene_isoform)
}

gene_umi_count <- function(cell_exon,edit,bar = "barcode",start = "start",umi = NULL,
                      seq = "search_seq",isoform = "exon_seq",polyA = "polyA",
                      sim_thresh = 5,split = "|",sep = ",",bar_len = 16, flank = 1,
                      umi_len = 10,splice_site_thresh = 10,
                      relation = mean2var,alpha = 0.05,verbose = FALSE){
    colnames(cell_exon)[which(colnames(cell_exon) == bar)] = "cell"
    colnames(cell_exon)[which(colnames(cell_exon) == isoform)] = "isoform"
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
        cell_exon$umi_flank = umi_flank
    }

    search_len = umi_len + 2*flank
    cell_exon <- cell_exon[nchar(cell_exon$umi_flank) == search_len,]
    if(nrow(cell_exon) == 0){
        return(NULL)
    }

    cells = unique(cell_exon[,"cell"])
    gene_cells_cluster <- lapply(cells,function(i){
        cell_i = cell_exon[cell_exon[,"cell"] == i,]
        if(verbose){
          cat(nrow(cell_i), " reads in cell ",i,"\n")
        }

        cell_i$cluster = 0
        if(nrow(cell_i) != 1){
          cell_i$cluster = umi_cluster_cpp(cell_i$umi_flank,
                           iso = cell_i$isoform,thresh = sim_thresh)
          #cell_i$cluster = umi_cluster_dbscan(cell_i$umi_flank,
          #                 iso = cell_i$isoform,thresh = sim_thresh)
        }
        else{
          cell_i$cluster = 1
        }

        return(cell_i)

    })
    gene_cells_cluster = as.data.frame(do.call(rbind,gene_cells_cluster))

    #return(gene_cells_cluster)
    gene_isoform = isoform_correct_filter(gene_cells_cluster,edit,
                                   split = split,sep = sep,
                                   splice_site_thresh = splice_site_thresh,
                                   relation = relation,alpha = alpha)

    return(gene_isoform)
}


umi_count <- function(cell_exon,edit = "edit",bar = "barcode",start = "start",gene = "gene_ID",
                      umi = NULL,seq = "search_seq",isoform = "exon_seq",polyA = "polyA",
                      sim_thresh = 5, split = "|",sep = ",",bar_len = 16,
                      flank = 1,umi_len = 10,
                      splice_site_thresh = 10,
                      relation = mean2var,alpha = 0.05){
    genes <- unique(cell_exon[,gene])
    edit = cell_exon[,edit]
    genes_umi_count <- lapply(genes,function(i){
        cat(i,"\n")
        sub_cell_exon = cell_exon[cell_exon[,gene] == i,]

        if(nrow(sub_cell_exon) < splice_site_thresh){
            cat("too few reads, will be filtered out\n")
            return(NULL)
        }
        sub_umi_count = gene_umi_count(sub_cell_exon,edit = edit,
                                       bar = bar,start = start,
                                       umi = umi,seq = seq,isoform = isoform,
                                       polyA = polyA,sim_thresh = sim_thresh,
                                       split = split,sep = sep,
                                       bar_len = bar_len, flank = flank,
                                       umi_len = umi_len,
                                       splice_site_thresh = splice_site_thresh,
                                       relation = relation,alpha = alpha)
        if(is.null(sub_umi_count)){
            return(NULL)
        }
        sub_umi_count$gene = i
        return(sub_umi_count)
    })

    genes_umi_count <- do.call(rbind,genes_umi_count)
    return(genes_umi_count)
}
