library(data.table)
library(stringr)
cpgs<-fread("analyses/2020-10-15_cpgs_to_validate.csv")

cpgs[,start:=as.integer(median(pos)-250),by="primerID"]
cpgs[,end:=as.integer(median(pos)+250),by="primerID"]
fwrite(unique(cpgs[,.(chr,start,end,primerID)]),file = "analyses/2020-10-15_cpgs_regions.bed",sep="\t",col.names = F)

cpgs[chr=="chr3"&start>200000000]

#get fasta : bedtools getfasta -fi ref/hg19/ucsc.hg19.fasta -bed analyses/2020-10-15_cpgs_regions.bed -name -tab > analyses/2020-10-15_cpgs_regions_fastas.tab

fastas<-fread("analyses/2020-10-15_cpgs_regions_fastas.tab",header = F,col.names = c("primerID","fasta"))
fastas[,primerID:=as.integer(str_extract(primerID,"^[0-9]+"))]
cpgs_fastas<-merge(cpgs,fastas)
fwrite(cpgs_fastas[,-"name"],"analyses/2020-10-15_cpgs_fastas_to_validate.csv",sep = ";")
