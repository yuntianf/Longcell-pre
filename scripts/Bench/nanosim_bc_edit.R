source("/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/scripts/Bench/nanosim.R")

wl = scan("~/data/BarcodeMatch/FT/FT_cells.txt",character(), quote = "")
#wl = scan("/D/Necessary/Upenn/Thesis/data/FT_gp5/FT_cells.txt",character(), quote = "")

cells =  scan("/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/bench_cells.txt",character(), quote = "")

trans = matrix(c(0.94,0.25,0.25,0.25,
                 0.02,0.25,0.25,0.25,
                 0.02,0.25,0.25,0.25,
                 0.02,0.25,0.25,0.25),4,4)

data = reads_simulator(cells,exprs = rep(10,length(cells)),
                       transcripts = c("",""),pcr = 20,trans = trans)
write.table(data[,c("read_name","reads")],
            file = "/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/softclips.txt",
            quote = FALSE,col.names = FALSE, row.names = FALSE,sep = "\t")

command = "/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/BarcodeMatch.sh"
system(command)

cos = read.table("/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/cos_bc.txt")
pos = read.table("/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/pos_bc.txt")

colnames(cos) = c("id","read_name","cos","start","edit","score")
colnames(pos) = c("id","read_name","pos","start","edit")

bc = inner_join(cos[,c("read_name","cos","edit")],pos[,c("read_name","pos")],by = "read_name")
bc = bc[bc$pos == bc$cos,]
bc$real = substr(bc$read_name,2,17)
sum(bc$real == bc$pos)/nrow(bc)

hist(bc$edit)
