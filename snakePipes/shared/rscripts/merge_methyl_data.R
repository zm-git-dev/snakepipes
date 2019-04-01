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
