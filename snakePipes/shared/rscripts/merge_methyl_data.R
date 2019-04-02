library(dplyr)

args = commandArgs(TRUE)

FilePath <- args[1]
outFile <- args[2]
print(getwd())
files = Sys.glob(paste0(FilePath,"/","*mean_methyl_per_region.tsv"))

all = data_frame("chr"=character(0),"start"=integer(0),"end"=integer(0))
 
for (f in files){
 dat=read.table(f,stringsAsFactors = F)
 name = gsub(".*/(.*).mean_methyl.*",replacement = "\\1",f)
 print(name)
 colnames(dat) = c("chr","start","end",name)
 all=merge(all,dat,all=T)
}

write.table(all,file = outFile, quote = F,row.names = F,col.names = T,sep = "\t")

dat = read.table(paste0(FilePath,"/on_target_stats.per_region.mapq20.tsv"),header = T,comment.char = "")
dat_mod = apply(dat[,-c(1:3)],2,function(x){x/sum(x)})
final= cbind(dat[,1:3],dat_mod)

write.table(final,file = "custom_stats/on_target_stats.per_region.mapq20.perc.tsv", quote = F,row.names = F,col.names = T,sep = "\t")

dat = read.table(paste0(FilePath,"/on_target_stats.per_region.tsv"),header = T,comment.char = "")
dat_mod = apply(dat[,-c(1:3)],2,function(x){x/sum(x)})
final= cbind(dat[,1:3],dat_mod)

write.table(final,file = "custom_stats/on_target_stats.per_region.perc.tsv", quote = F,row.names = F,col.names = T,sep = "\t")

