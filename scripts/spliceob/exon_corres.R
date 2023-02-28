#' @title splice site extraction
#' @description  Extract splice sites from a splice site sequence (used to represent an isoform).
#' @details Input a splice site sequence and return the splice sites within it.
#' @param isoform A string representing isoform.
#' @param split A character to seperate exons within the isoform.
#' @param sep A character to seperate the start and end site for an exon.
#' @return A numeric vector of splice sites.
#' @examples
#' splice_site("1000,2000|2500,3000")
splice_site <- function(isoform,split = "|",sep = ","){
  exons = unlist(strsplit(isoform,split = split,fixed = T))
  splice = unlist(strsplit(exons,split = sep))

  return(splice)
}

#' @title read start sites extraction
#' @description  Extract start sites from multiple splice site sequences (used to represent an isoform).
#' @details Input a vector of splice site sequences and return the start sites for each of them.
#' @param isoforms A vector of strings representing isoforms.
#' @param split A character to seperate exons within the isoform.
#' @param sep A character to seperate the start and end site for an exon.
#' @return A numeric vector of start sites.
#' @examples
#' transcript_start("1000,2000|2500,3000")
transcript_start <- function(isoforms,split = "|",sep = ","){
  start = sapply(isoforms,function(x){
    sites = splice_site(x,split,sep)
    return(as.numeric(sites[1]))
  })
  names(start) = NULL
  return(start)
}

#' @title read end sites extraction
#' @description  Extract end sites from multiple splice site sequences (used to represent an isoform).
#' @details Input a vector of splice site sequences and return the end sites for each of them.
#' @param isoforms A vector of strings representing isoforms.
#' @param split A character to seperate exons within the isoform.
#' @param sep A character to seperate the start and end site for an exon.
#' @return A numeric vector of end sites.
#' @examples
#' transcript_end("1000,2000|2500,3000")
transcript_end <- function(isoforms,split = "|",sep = ","){
  end = sapply(isoforms,function(x){
    sites = splice_site(x,split,sep)
    return(as.numeric(sites[length(sites)]))
  })
  names(end) = NULL
  return(end)
}

#' @title transform exon bins to exon id
#' @description  correspond an exon bin to the gene bed annotation.
#' @details Input a vector of start and end for an exon and return the correponding id for this exon in the gene bed annotation.
#' @param bin A vector of start and end for an exon.
#' @param gene_bed A dataframe for the gene bed annotation.
#' @param exon_id A vector of strings to self marked the exon name in the gene bed.
#' @param start The name of the column in the gene bed indicating the start sites of each exon.
#' @param end The name of the column in the gene bed indicating the end sites of each exon.
#' @param bias The maximum tolerance for the bias of start and end sites correspondence for an exon.
#' @return A string vector of exon ids.
#' @examples
#' bin2exon("1000,2000",gene_bed=gene_bed)
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

#' @title transform a splice site sequence to exon id sequence
#' @description  correspond exon bins in the read to the gene bed annotation.
#' @details Input a splice site sequence and return the correponding id for exons in this read according to the gene bed annotation.
#' @param splice_sites A string representing isoform.
#' @param gene_bed A dataframe for the gene bed annotation.
#' @param exon_id A vector of strings to self marked the exon name in the gene bed.
#' @param split A character to separate exons within the isoform.
#' @param sep A character to separate the start and end site for an exon.
#' @param start The name of the column in the gene bed indicating the start sites of each exon.
#' @param end The name of the column in the gene bed indicating the end sites of each exon.
#' @param strand The name for the column indicating the strand for a gene in the gene bed.
#' @param bias The maximum tolerance for the bias of start and end sites correspondence for an exon.
#' @return A string of exon id sequence.
#' @examples
#' splice2exon("1000,2000|2500,3000",gene_bed=gene_bed)
splice2exon <- function(splice_sites,gene_bed,exon_id = NULL,
                        split = "|",sep = ",",start = "start",end = "end",
                        strand = "strand",bias = 5){
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
                    start = start,end = end,bias = bias)
    return(exon)
  })
  exons = unlist(exons)
  isoform = paste(exons,collapse = split)
  return(isoform)
}

#' @title transform multiple splice site sequences to exon id sequences
#' @description  correspond exon bins in the read to the gene bed annotation.
#' @details Input a splice site sequence and return the correponding id for exons in this read according to the gene bed annotation.
#' @param splice_sites A string representing isoform.
#' @param gene_bed A dataframe for the gene bed annotation.
#' @param exon_id A vector of strings to self marked the exon name in the gene bed.
#' @param split A character to seperate exons within the isoform.
#' @param sep A character to seperate the start and end site for an exon.
#' @param start The name of the column in the gene bed indicating the start sites of each exon.
#' @param end The name of the column in the gene bed indicating the end sites of each exon.
#' @param strand The name for the column indicating the strand for a gene in the gene bed
#' @param bias The maximum tolerance for the bias of start and end sites correspondence for an exon.
#' @return A string vector of exon ids.
#' @export
#' @examples
#' isoforms2exon(c("1000,2000|2500,3000","1100,2000|2500,3200"),gene_bed=gene_bed)
isoforms2exon <- function(isoforms,gene,gene_bed,exon_id = NULL,
                          split = "|",sep = ",",gene_col = "gene",
                          start = "start",end = "end",strand = "strand",
                          bias = 5){
  gene_bed = gene_bed[gene_bed[,gene_col] == gene,]
  exon_seqs <- sapply(isoforms,function(x){
    seq = splice2exon(splice_sites = x,gene_bed = gene_bed,
                      exon_id = exon_id,split = split, sep = sep,
                      start = start,end = end,strand = strand,bias = bias)
    return(seq)
  })
  names(exon_seqs) <- NULL
  return(exon_seqs)
}

#' @title transform splice site sequences to exon id sequences to represent isoform
#' @description  correspond exon bins in the read to the gene bed annotation.
#' @details Input a data frame for single cell isoform count and transform the isoform from splice sites format into exon id format.
#' @param data A dataframe , each row represents an isoform in a gene in a cell.
#' @param cell Selected cells from the data.
#' @param gene Selected genes from the data
#' @param cell_col The name of the column indicating cells in the data
#' @param gene_col The name of the column indicating genes in the data
#' @param isoform_col The name of the column indicating isoforms in the data
#' @param gene_bed A dataframe for the gene bed annotation.
#' @param split A character to seperate exons within the isoform.
#' @param sep A character to seperate the start and end site for an exon.
#' @param bed_start The name of the column in the gene bed indicating the start sites of each exon.
#' @param bed_end The name of the column in the gene bed indicating the end sites of each exon.
#' @param bed_strand The name for the column indicating the strand for a gene in the gene bed
#' @param bias The maximum tolerance for the bias of start and end sites correspondence for an exon.
#' @return A string vector of exon ids.
#' @export
createExonList <- function(data,gene_bed,cells = NULL,genes = NULL,
                           cell_col = "cell",gene_col = "gene",
                           isoform_col = "isoform",split = "|",sep = ",",
                           bed_start = "start", bed_end = "end",
                           bed_strand = "strand",bias = 5){
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
                          start = bed_start, end = bed_end,strand = bed_strand,
                          bias = bias)

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

#' @title save the ExonList
#' @description  Save the ExonList dataframe as sparse matrix in designated path.
#' @details save the ExonList dataframe as sparse matrix in designated path.
#' @param data An ExonList dataframe created by creatExonList().
#' @param path The path to save the data in
#' @param cell_col The name of the column indicating cells in the data
#' @param gene_col The name of the column indicating genes in the data
#' @param isoform_col The name of the column indicating isoforms in the data
#' @param count_col The name of the column indicating count for the isoform in the data
#' @return NULL
#' @export
#' @import dplyr
saveExonList <- function(data,path = "./",cell_col = "cell",gene_col = "gene",
                         isoform_col = "isoform",count_col = "count"){
  data = data[order(data[,cell_col],data[,gene_col],data[,isoform_col]),]

  gene_count = data %>%
               group_by(!!as.name(cell_col),!!as.name(gene_col)) %>%
               summarise(count = sum(get(count_col)),.groups = "drop")
  gene_count = as.data.frame(gene_count)

  cells = names(table(gene_count[,cell_col]))
  genes = names(table(gene_count[,gene_col]))

  gene_count[,cell_col] = as.numeric(as.factor(gene_count[,cell_col]))
  gene_count[,gene_col] = as.numeric(as.factor(gene_count[,gene_col]))

  write.table(cells,file = paste(path,"barcodes.tsv",sep = "/"),
              row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(genes,file = paste(path,"features.tsv",sep = "/"),
              row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(gene_count,file = paste(path,"matrix.mtx",sep = "/"),
              row.names = FALSE,col.names = FALSE,quote = FALSE)

  isoform = names(table(data[,isoform_col]))
  write.table(isoform,file = paste(path,"isoforms.tsv",sep = "/"),
              row.names = FALSE,col.names = FALSE,quote = FALSE)

  data[,cell_col] = as.numeric(as.factor(data[,cell_col]))
  data[,gene_col] = as.numeric(as.factor(data[,gene_col]))
  data[,isoform_col] = as.numeric(as.factor(data[,isoform_col]))
  write.table(data,file = paste(path,"isoform_count.tsv",sep = "/"),
              row.names = FALSE,col.names = TRUE,quote = FALSE)
  return(NULL)
}
