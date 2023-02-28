library(ggplot2)
library(dplyr)

bc_judge = function(name,barcode,start = 2,end = 17){
  real = substr(name,start,end)
  return(real == barcode)
}

setwd("C:/D/Necessary/Upenn/Thesis/projects/Single-cell-long-reads/data/preprocess/simulation/MYL6_201_202_bc/")

prefix = "MYL6_201_202_"
suffix = "_name.txt"

vars = c(1:5)
##### valid reads #####
file = paste("./Longcell/",prefix,vars,"/softclips/softclips.txt",sep = "")
softclips = lapply(1:length(vars),function(i){
  tmp = read.table(file[i])
  colnames(tmp) = c("name","softclips")
  tmp$name = paste("@",tmp$name,"_",vars[i],sep = "")
  tmp$vars = vars[i]
  return(tmp)
})
softclips = do.call(rbind,softclips)

##### Longcell #####
file = paste("./Longcell/",prefix,vars,"/barcode_match/pos_bc.txt",sep = "")
pos = lapply(1:length(vars),function(i){
  tmp = read.table(file[i])[,2:5]
  colnames(tmp) = c("name","barcode","start","edit")
  tmp$name = paste("@",tmp$name,"_",vars[i],sep = "")
  tmp$vars = vars[i]
  return(tmp)
})
pos = do.call(rbind,pos)
file = paste("./Longcell/",prefix,vars,"/barcode_match/cos_bc.txt",sep = "")
cos = lapply(1:length(vars),function(i){
  tmp = read.table(file[i])[,2:5]
  colnames(tmp) = c("name","barcode","start","edit")
  tmp$name = paste("@",tmp$name,"_",vars[i],sep = "")
  tmp$vars = vars[i]
  return(tmp)
})
cos = do.call(rbind,cos)
longcell = left_join(pos,cos[,c(1,2,5)],by = c("name","vars"))
longcell = na.omit(longcell)
longcell = longcell[longcell$barcode.x == longcell$barcode.y,]
colnames(longcell)[2] = "barcode"
longcell = longcell[,-6]
##### sicelore2 #####
file = paste("./sicelore2/",prefix,vars,suffix,sep = "")
sicelore2 = lapply(1:length(vars),function(i){
  tmp = read.table(file[i])
  colnames(tmp) = c("name","barcode")
  tmp$name = sapply(tmp$name,function(y){
    strsplit(y,split = "_FWD_PS",fixed = TRUE)[[1]][1]
  })
  tmp$name = paste(tmp$name,"_",vars[i],sep = "")
  tmp$barcode = substr(tmp$barcode,8,23)
  tmp$vars = vars[i]
  return(tmp)
})
sicelore2 = do.call(rbind,sicelore2)

##### FLAMES #####
file = paste("./FLAMES/",prefix,vars,suffix,sep = "")
FLAMES = lapply(1:length(vars),function(i){
  tmp = read.table(file[i],sep = "#",comment.char = "/")
  colnames(tmp) = c("barcode","name")
  tmp$barcode = substr(tmp$barcode,2,17)
  tmp$name = paste('@',tmp$name,"_",vars[i],sep = "")
  tmp$vars = vars[i]
  return(tmp)
})
FLAMES = do.call(rbind,FLAMES)




##### mean edit distance #####
adapter = left_join(softclips,longcell[longcell$edit == 0,c("name","start")],by = "name")
adapter = na.omit(adapter)
adapter$adapter = substr(adpater$softclips,adpater$start-15,adpater$start)
adapter = adapter[nchar(adapter$adapter) == 16,]
adapter = adapter[,c("vars","adapter")]
adapter$edit = unlist(adist(adapter$adapter,"GACGCTCTTCCGATCT"))
adapter_med = adapter %>% group_by(vars) %>% summarise(edit = mean(edit))



##### comparison #####
softclips_count = table(softclips$vars)

longcell = longcell[longcell$name %in% softclips$name,]
sicelore2 = sicelore2[sicelore2$name %in% softclips$name,]
FLAMES = FLAMES[FLAMES$name %in% softclips$name,]

longcell$correct = bc_judge(longcell$name,longcell$barcode)
sicelore2$correct = bc_judge(sicelore2$name,sicelore2$barcode)
FLAMES$correct = bc_judge(FLAMES$name,FLAMES$barcode)

longcell_acc = longcell %>% group_by(vars) %>% summarise(acc = sum(correct)/n())
sicelore2_acc = sicelore2 %>% group_by(vars) %>% summarise(acc = sum(correct)/n())
FLAMES_acc = FLAMES %>% group_by(vars) %>% summarise(acc = sum(correct)/n())

longcell_pre = table(longcell$vars)/softclips_count
sicelore2_pre = table(sicelore2$vars)/softclips_count
FLAMES_pre = table(FLAMES$vars)/softclips_count

longcell_acc$pre = as.numeric(longcell_pre)
sicelore2_acc$pre = as.numeric(sicelore2_pre)
FLAMES_acc$pre = as.numeric(FLAMES_pre)

longcell_acc$f1 = longcell_acc$acc*longcell_acc$pre/(longcell_acc$acc+longcell_acc$pre)
sicelore2_acc$f1 = sicelore2_acc$acc*sicelore2_acc$pre/(sicelore2_acc$acc+sicelore2_acc$pre)
FLAMES_acc$f1 = FLAMES_acc$acc*FLAMES_acc$pre/(FLAMES_acc$acc+FLAMES_acc$pre)

sub_cos_acc$f1 = sub_cos_acc$acc*sub_cos_acc$pre/(sub_cos_acc$acc+sub_cos_acc$pre)
sub_pos_acc$f1 = sub_pos_acc$acc*sub_pos_acc$pre/(sub_pos_acc$acc+sub_pos_acc$pre)

comp = rbind(sicelore2_acc,FLAMES_acc,sub_cos_acc[,-3])
comp$tool = rep(c("sicelore2","FLAMES","Longcell"),each = 5)

comp = left_join(comp,adapter_med,by = "vars")

ggplot(pos_acc[pos_acc$edit > 1,])+
  geom_point(aes(x = diff,y = edit,color = acc,size = count))+
  theme_classic()

bc_acc_comp = ggplot(comp)+
  geom_point(aes(x = edit,y = acc,color = tool,group = tool))+
  geom_line(aes(x = edit,y = acc,color = tool,group = tool))+
  theme_classic()+xlab("mean edit distance")+ylab("accuracy")+
  labs(color='method') +theme(text = element_text(size = 18))
ggsave(plot = bc_acc_comp, filename = "../../../../../../results/bench/bc_acc_comp.png",
       device = "png", width = 6, height = 5)

bc_pre_comp = ggplot(comp)+
  geom_point(aes(x = edit,y = pre,color = tool,group = tool))+
  geom_line(aes(x = edit,y = pre,color = tool,group = tool))+
  theme_classic()+xlab("mean edit distance")+ylab("preserve ratio")+
  labs(color='method')+theme(text = element_text(size = 18))
ggsave(plot = bc_pre_comp, filename = "../../../../../../results/bench/bc_pre_comp.png",
       device = "png", width = 6, height = 5)