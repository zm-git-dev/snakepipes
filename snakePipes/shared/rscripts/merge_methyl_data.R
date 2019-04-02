library(dplyr)

args = commandArgs(TRUE)

FilePath = "custom_stats.dataset3b/"

FilePath <- args[1]
outFile <- args[2]
print(getwd())

## mean_methyl files arew from MethylDackel and bedtools map with target regions
files = Sys.glob(paste0(FilePath,"/","*mean_methyl_per_region.tsv"))

all = data_frame("chr"=character(0),"start"=integer(0),"end"=integer(0))
 
for (f in files){
 dat=read.table(f,stringsAsFactors = F)
 name = gsub(".*/(.*).mean_methyl.*",replacement = "\\1",f)
 print(name)
 colnames(dat) = c("chr","start","end",name)
 all=merge(all,dat,all=T)
}

write.table(format(all, digits=2, scientific=F),file = outFile, quote = F,row.names = F,col.names = T,sep = "\t")

dat = read.table(paste0(FilePath,"/on_target_stats.per_region.mapq20.tsv"),header = T,comment.char = "")
dat_mod = apply(dat[,-c(1:3)],2,function(x){x/sum(x)})*100
final= cbind(dat[,1:3],dat_mod)

write.table(format(final, digits=2, scientific=F),file = "custom_stats/on_target_stats.per_region.mapq20.perc.tsv", quote = F,row.names = F,col.names = T,sep = "\t")

dat = read.table(paste0(FilePath,"/on_target_stats.per_region.tsv"),header = T,comment.char = "")
dat_mod = apply(dat[,-c(1:3)],2,function(x){x/sum(x)})*100
final= cbind(dat[,1:3],dat_mod)

write.table(format(final, digits=2, scientific=F),file = "custom_stats/on_target_stats.per_region.perc.tsv", quote = F,row.names = F,col.names = T,sep = "\t")

#pheatmap::pheatmap((apply(dat[,-c(1:3)],2,function(x){x/sum(x)})*100),scale = 'none')

####################
tcpgs = read.table(paste0(FilePath,"targets.CpG.bed"),header = F)

colnames(tcpgs) = c("chr","start","end","name","score","strand","context")
tcpgs$end = tcpgs$end+1
  
files = Sys.glob(paste0(FilePath,"*_CpG.bedGraph"))

all = data_frame("chr"=character(0),"start"=integer(0),"end"=integer(0))

for (f in files){
  dat=read.table(f,stringsAsFactors = F,skip = 1)
  name = gsub(".*/(.*)_CpG.bedGraph",replacement = "\\1",f)
  print(name)
  colnames(dat) = c("chr","start","end",name,"M","U")
  all=merge(all,dat[,1:4],all=T)
}

write.table(all,file = paste0(FilePath,"Percent_Methylation_by_CpG.tsv"), quote = F,row.names = F,col.names = T,sep = "\t")


