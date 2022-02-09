

library(data.table)
#get fasta of target regions (non bisulf convert)
fastas<-fread("outputs/bisulf_pcr/2020-10-15_cpgs_regions_fastas.tab",header = F,col.names =c("id","sequence") )
fastas[,region:=str_remove(id,"^[0-9]+::")]

cpg_info<-fread("outputs/bisulf_pcr/2020-10-15_cpgs_fastas_to_validate.csv")
cpg_info[,region:=paste(chr,paste(start,end,sep="-"),sep=":")]
fastas<-merge(fasta,unique(cpg_info[,.(region,gene)]),by="region")
fastas
fastas[,id:=paste0(">",paste(region,gene,sep="|"))]
char_vec<-c()
for(i in 1:nrow(fastas)){
  id<-fastas[i,]$id
  seq<-fastas[i,]$seq
  char_vec<-c(char_vec,id,seq)
}
writeLines(char_vec,"outputs/bisulf_pcr/2020-10-15_cpgs_regions.fa")

# bisulf convertis fasta ref with bismark
#need add bowtie2 in PATH
system("tools/bismark/Bismark-0.23.0/bismark_genome_preparation --verbose outputs/bisulf_pcr") #work on the shell

#tolerating one non-bisulfite mismatch per read

#in shell on outputs/bisulf_pcr/fastq/ directory
system("../../../tools/bismark/Bismark-0.23.0/bismark --bowtie2 -n 1 -l 50 ../ -1 *_R1.fastq -2 *_R2.fastq")


