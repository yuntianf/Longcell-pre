exon_part_merge <- function(bed,len = 10,start="start",end = "end",exon = "exon_id"){
    bed$length <- bed[,end]-bed[,start] + 1
    ts <- which(bed$length < len)
    
    if(length(ts) == 0){
        if(nrow(bed) > 1){
            if(bed[1,exon] < bed[2,exon]){
                bed[,exon] <- c(1:nrow(bed))
            }
            else{
                bed[,exon] <- rev(c(1:nrow(bed)))
            }
        }
        return(bed)
    }
    else{
        i = ts[1]
        if(i == 1){
            if(bed[i+1,start] - 1 == bed[i,end]){
                bed[i+1,start] <- bed[i,start]
            }
        }
        else if(i == nrow(bed)){
            if(bed[i-1,end] + 1 == bed[i,start]){
                bed[i-1,end] <- bed[i,end]
            }
        }
        else{
            if((bed[i-1,end] + 1 == bed[i,start]) & (bed[i+1,start] - 1 == bed[i,end])){
                if(bed[i-1,]$len <= bed[i+1,]$len){
                    bed[i-1,end] <- bed[i,end]
                }
                else{
                    bed[i+1,start] <- bed[i,start]
                }
            }
            else if(bed[i-1,end] + 1 == bed[i,start]){
                bed[i-1,end] <- bed[i,end]
            }
            else if(bed[i+1,start] - 1 == bed[i,end]){
                bed[i+1,start] <- bed[i,start]
            }
        }
        
        bed <- bed[-i,]
        return(exon_part_merge(bed,len,start,end))
    }
}

args <- commandArgs(trailingOnly = TRUE)
bed <- read.table(args[1])

colnames(bed) <- c("chrom","start","end","exon_id")
merged_bed <- exon_part_merge (bed)
write.table(merged_bed ,quote= F,row.names = F, col.names = F,file = args[2])