library(Matrix)
library(dplyr)
library(geomtextpath)

normalize <- function(x){
  return(log(x/sum(x,na.rm = TRUE)*10000 + 1))
}

MSE <- function(x,y,if_norm = TRUE){
  if(if_norm){
    x = normalize(x)
    y = normalize(y)
  }
  return(mean((x-y)^2))
  
}

exon_corres <- function(start,end,bed,exon_id = NULL,
                        strand_col = "strand",start_col = "start",end_col = "end"){
  if(is.null(exon_id)){
    if(unique(bed[,strand_col]) == "+"){
      exon_id = 1:nrow(bed)
    }
    else{
      exon_id = nrow(bed):1
    }
  }
  else{
    if(length(exon_id) != nrow(bed)){
      stop("The exon id annotation should have the same length of gene bed!")
    }
  }
  
  exons = exon_id[which(bed[,start_col] >= start & bed[,end_col] <= end)]
  return(exons)
}

transcript_corres <- function(starts,ends,bed,exon_id = NULL,
                              strand_col = "strand",start_col = "start",end_col = "end"){
  if(length(starts) != length(ends)){
    stop("The start and end point of exons don't correspond!")
  }
  
  exons <- lapply(1:length(starts),function(i){
    exon = exon_corres(starts[i],ends[i],bed,
                       exon_id,strand_col,start_col,end_col)
    return(exon)
  })
  exons <- unlist(exons)
  names(exons) <- NULL
  
  if(unique(bed[,strand_col]) == "+"){
    exons = sort(exons)
  }
  else{
    exons = sort(exons,decreasing = TRUE)
  }
  exons = paste(exons,collapse = "|")
  return(exons)
}

gtf_bed_corres <- function(gtf,bed,gene,
                           gtf_gene_col = "gene_id",bed_gene_col = "gene",
                           transcript_col = "transname",
                           exon_id = NULL,strand_col = "strand",
                           gtf_start_col = "start",gtf_end_col = "end",
                           bed_start_col = "start",bed_end_col = "end"){
  gene_gtf = gtf[gtf[,gtf_gene_col] == gene,]
  gene_bed = bed[bed[,bed_gene_col] == gene,]
  
  transcript_exons <- gene_gtf %>% 
    group_by(across(all_of(transcript_col))) %>% 
    summarise(exons = transcript_corres(starts = !!sym(gtf_start_col),
                                        ends = !!sym(gtf_end_col),
                                        bed = gene_bed,
                                        exon_id = NULL,
                                        strand_col = strand_col,
                                        start_col = bed_start_col,
                                        end_col = bed_end_col))
  
  return(transcript_exons)
}

beta_mean <- function(alpha,beta){
    return(alpha/(alpha+beta))
}

beta_var <- function(alpha,beta){
    return(alpha*beta/((alpha+beta)^2*(alpha+beta+1)))
}

beta_phi <- function(alpha,beta){
    temp_var = beta_var(alpha,beta)
    temp_mean = beta_mean(alpha,beta)
    
    return(temp_var/(temp_mean*(1-temp_mean)))
}

isoform2bases <- function(isoform,count = 1,
                          sep = ",",split = "|"){
  exons = unlist(strsplit(isoform,split = split,fixed = TRUE))
  
  bases = lapply(exons,function(x){
    x = unlist(strsplit(x,split = sep,fixed = TRUE))
    left = as.numeric(x[1])
    right = as.numeric(x[2])
    return(c(left:right))
  })
  bases = rep(unlist(bases),count)
  return(bases)
}

junction_count <- function(isoform,count = 1,polyA = 0,
                           sep = ",",split = "|"){
  junctions = unlist(strsplit(isoform,split = sep,fixed = TRUE))
  end = junctions[length(junctions)]
  if(length(junctions) > 2){
    junctions = junctions[2:(length(junctions)-1)]
    junc_count = lapply(junctions,function(x){
      x = unlist(strsplit(x,split = split,fixed = TRUE))
      return(x)
    })
    junc_count = do.call(rbind,junc_count)
    junc_count = cbind(junc_count,count)
  }
  else{
    junc_count = NULL
  }
  if(polyA > 0.5){
    end_count = c(end,"polyA",count)
  }
  else{
    end_count = NULL
  }
  junc_count = rbind(junc_count,end_count)
  return(junc_count)
}

sashimi_plot_data <- function(isoforms,counts,polyA = NULL,
                              sep = ",",split = "|"){
  if(length(isoforms) != length(counts)){
    stop("The size of isoforms and their counts don't match!")
  }
  if(!is.null(polyA)){
    if(length(polyA) != length(isoforms)){
      stop("The size of isoforms and their polyA don't match!")
    }
  }
  else{
    polyA = rep(0,length(isoforms))
  }
  
  iso_count = as.data.frame(cbind(isoforms,counts))
  colnames(iso_count) = c("isoform","count")
  iso_count$count = as.numeric(iso_count$count)
  iso_count = iso_count %>% 
              group_by(isoform) %>%
              summarise(count = sum(count))
  bases = lapply(1:nrow(iso_count),function(i){
    isoform = iso_count$isoform[i]
    count = iso_count$count[i]
    sub_base = isoform2bases(isoform = isoform,count = count,
                             sep = sep,split = split)
    return(sub_base)
  })
  bases = as.data.frame(table(unlist(bases)))
  colnames(bases) = c("chr","coverage")
  bases$chr = as.numeric(as.character(bases$chr))
  
  junctions = lapply(1:length(isoforms),function(i){
    sub_junc = junction_count(isoform = isoforms[i],count = counts[i],
                              polyA = polyA[i],sep = sep,split = split)
    return(sub_junc)
  })
  junctions = as.data.frame(do.call(rbind,junctions))
  colnames(junctions) = c("start","end","count")
  junctions$count = as.numeric(junctions$count)
  junctions = junctions %>% 
              group_by(start,end) %>% 
              summarise(count = sum(count),.groups = "drop")
  
  return(list(bases,junctions))
}

range_scale = function(data,lwr = 1,upr = 5){
  diff = max(data)-min(data)
  unit = (upr-lwr)/diff
  scale_data = (data-min(data))*unit+lwr
  return(scale_data)
}

sashimi_plot <- function(coverage,junction,filter_ratio = 20,
                         color_id = 1,region = NULL,lwr = 1,upr = 5,
                         color_set = NULL){
  rownames(coverage) = coverage$chr
  
  junction$y_start = coverage[junction$start,"coverage"]
  junction$y_end = coverage[junction$end,"coverage"]
  junction = na.omit(junction)
  junction$start = as.numeric(junction$start)
  junction$end = as.numeric(junction$end)
  
  filter = max(coverage$coverage)/filter_ratio
  junction_filter = junction[junction$count > filter & 
                               junction$y_start > filter & 
                               junction$y_end > filter,]
  
  if(is.null(region)){
    region = "chr"
  }
  if(is.null(color_set)){
    color_set = list(c("steelblue","LightSkyBlue"),c("coral","LightSalmon"))
  }
  
  sashimi = ggplot()+
    geom_segment(data = coverage,aes(x = chr,y = coverage,xend = chr,yend = 0),
                 color = color_set[[color_id]][1],alpha = 0.5)+
    geom_textcurve(data = junction_filter,aes(x = start,y = y_start,xend = end,yend = y_end,label = round(count),
                                              linewidth = range_scale(log(count,10),lwr = lwr,upr = upr)),
                   linecolour = color_set[[color_id]][2],lineend = "butt",curvature = -0.1,textcolour = "black",
                   alpha = 0.75)+
    theme_classic()+
    xlab(region)+ylab("reads coverage")+ylim(c(0,max(coverage$coverage)*1.3))+
    theme(text = element_text(size = 13))
  
  return(sashimi)
}

exon_include <- function(isoform,exons_in = NULL,exons_out = NULL,
                         thresh = 0,sep = ",",split = "|"){
  if(!is.null(exons_in)){
    exon_in_bases = lapply(exons_in,function(x){
      return(isoform2bases(x,sep = sep, split = split))
    })
    exon_in_bases = unlist(exon_in_bases)
  }
  else{
    exon_in_bases = c()
  }
  if(!is.null(exons_out)){
    exon_out_bases = lapply(exons_out,function(x){
      return(isoform2bases(x,sep = sep, split = split))
    })
    exon_out_bases = unlist(exon_out_bases)
  }  
  
  isoform_bases = isoform2bases(isoform)
  lack = setdiff(exon_in_bases,isoform_bases)
  with = intersect(exon_out_bases,isoform_bases)
  
  if(length(lack) <= thresh & length(with) <= thresh){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

### unfinshed
exonsLen <- function(exons,bed,start,end,exon_id = NULL,
                     strand_col = "strand",start_col = "start",end_col = "end"){
  if(is.null(exon_id)){
    if(unique(bed[,strand_col]) == "+"){
      bed$id = as.character(1:nrow(bed))
    }
    else{
      bed$id = as.character(nrow(bed):1)
    }
  }
  else{
    bed$id = exon_id
  }
  if("N" %in% exons){
    warning("There exist unannotated exons, they won't be calculated")
    exons = setdiff(exons,"N")
  }
  sub_bed = bed[bed$id %in% exons]
  min_start = min(bed[,start_col])
  max_end = max(bed[,end_col])
  
  total_len = sum(bed[,end_col]-bed[,start_col]+1)
  if(!is.null(start)){
    len = total_len - (start - min_start)
  }
  if(!is.null(end)){
    len = total_len - (max_end - end)
  }
  return(len)
}

exonSeqSource <- function(exonSeq,annoExonSeq,bed,thresh = 10,
                          start = NULL,end = NULL,polyA = 0,
                          strand_col = "strand",split = "|"){
  exonSeq = unlist(strsplit(exonSeq,split = split,fixed = TRUE))
  annoExonSeq = unlist(strsplit(annoExonSeq,split = split,fixed = TRUE))
  
  exon_dic = sort(as.numeric(unique(exonSeq,annoExonSeq)))
  exo_vec = rep(0,length(exon_dic))
  
  
}

exonSeq2isoform <- function(exonSeq,bed,gtf,gene,start,end,
                            bed_gene_col = "gene",gtf_gene_col = "gene_id"){
  gene_bed = bed[bed[,bed_gene_col] == gene,]
  corres = gtf_bed_corres(gtf,gene_bed,gene,
                          gtf_gene_col = gtf_gene_col,
                          bed_gene_col = bed_gene_col)
  
  exon_corres_list = lapply(1:length(exonSeq),function(i){
    exon_corres_vec = sapply(corres$exons,function(y){
      return(exonSeqSource(exonSeq[i],y,gene_bed,thresh = 10,start[i],end[i]))
    })
    isoforms = corres$transname[exon_corres_vec]
    return(isoforms)
  })
  
  return(exon_corres_list)
}

mem <- function(variable){
  m = format(object.size(variable), units = "Mb")
  return(m)
}