library(dplyr)

files = Sys.glob("custom_stats/*mean_methyl_per_region.tsv")

all = data_frame("chr"=character(0),"start"=integer(0),"end"=integer(0))
 
for (f in files){
 dat=read.table(f,stringsAsFactors = F)
 name = gsub("custom_stats/(.*).mean_methyl.*",replacement = "\\1",f)
 print(name)
 colnames(dat) = c("chr","start","end",name)
 all=merge(all,dat,all=T)
}

write.table(all,file = "mean_methyl_per_region.tsv",quote = F,row.names = F,col.names = T,sep = "\t")
