
library(data.table)
library(stringr)
source("scripts/utils/new_utils.R")

out<-"outputs/bisulf_cbps_cpgs"
dir.create(out)
fasta_dir<-"ref/lga_cpgs_fastas/"
fastq_dir<-"~/RUN/Run_654_rna_methylation/Output/FastQ"
samples<-unique(str_extract(list.files(fastq_dir,pattern = "CBP"),"CBP[0-9]+"))


#get fasta of target regions (non bisulf convert)

# bisulf convertis fasta ref with bismark
system(paste0("tools/bismark/Bismark-0.23.0/bismark_genome_preparation --path_to_aligner /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --verbose ",fasta_dir)) 


for(s in samples){
  print(s)
  #need decrompress fastq
  dir.create(file.path(out,"",s,"fastq"),recursive = T)
  system(paste0("bzip2 -dkc  ",fastq_dir,"/",s,"_R1.fastq.bz2 > ",out,"/",s,"/fastq/",s,"_R1.fastq"))
  system(paste0("bzip2 -dkc  ",fastq_dir,"/",s,"_R2.fastq.bz2 > ",out,"/",s,"/fastq/",s,"_R2.fastq"))
  
  
  #bismark tolerating one non-bisulfite mismatch per read
  
  
  system(paste0("tools/bismark/Bismark-0.23.0/bismark --genome ",fasta_dir," --local --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o ",out,"/",s,"  -n 1 -l 20 ",out,"/",s,"/fastq/",s,"_R1.fastq"))
  system(paste0("tools/bismark/Bismark-0.23.0/bismark --genome ",fasta_dir," --local --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o ",out,"/",s,"  -n 1 -l 20 ",out,"/",s,"/fastq/",s,"_R2.fastq"))
  
  
  system(paste0("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o ",out,"/",s," ",out,"/",s,"/",s,"_R1_bismark_bt2.bam"))
  system(paste0("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o ",out,"/",s," ",out,"/",s,"/",s,"_R2_bismark_bt2.bam"))
  
  system(paste0("tools/bismark/Bismark-0.23.0/bismark2report --dir ",out,"/",s," --alignment_report ",out,"/",s,"/",s,"_R1_bismark_bt2_SE_report.txt --dedup_report NONE --splitting_report ",out,"/",s,"/",s,"_R1_bismark_bt2_splitting_report.txt --mbias_report ",out,"/",s,"/",s,"_R1_bismark_bt2.M-bias.txt"))
  system(paste0("tools/bismark/Bismark-0.23.0/bismark2report --dir ",out,"/",s," --alignment_report ",out,"/",s,"/",s,"_R2_bismark_bt2_SE_report.txt --dedup_report NONE --splitting_report ",out,"/",s,"/",s,"_R2_bismark_bt2_splitting_report.txt --mbias_report ",out,"/",s,"/",s,"_R2_bismark_bt2.M-bias.txt"))
}

#CBP519 don't work #no reads
samples<-samples[samples!="CBP519"]

res_cpgs_r1<-Reduce(rbind,lapply(samples,function(s)fread(paste0(out,"/",s,"/CpG_OT_",s,"_R1_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R1_OT"][,cpg:=T]))
res_cpgs_r2<-Reduce(rbind,lapply(samples,function(s)fread(paste0(out,"/",s,"/CpG_CTOT_",s,"_R2_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R2_CTOT"][,cpg:=T]))
res_cpgs<-rbind(res_cpgs_r1,res_cpgs_r2)

res_c_r1<-Reduce(rbind,lapply(samples,function(s)fread(paste0(out,"/",s,"/Non_CpG_OT_",s,"_R1_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R1_OT"][,cpg:=F]))
res_c_r2<-Reduce(rbind,lapply(samples,function(s)fread(paste0(out,"/",s,"/Non_CpG_CTOT_",s,"_R2_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R2_CTOT"][,cpg:=F]))
res_cs<-rbind(res_c_r1,res_c_r2)

res<-rbind(res_cpgs,res_cs)
res[,depth:=.N,by=c("chr","pos","sample","cpg")]
res[,depth.glob:=.N,by=c("chr","pos","cpg")]
res[,n.read:=.N,by=c("chr","pos","sample","cpg")]
res[,n.meth:=sum(meth=="+"),by=c("chr","pos","sample","cpg")]
res[,pct.meth:=n.meth/n.read]


fwrite(res,file.path(out,"res_Cxx_meth_all_reads.csv.gz"),sep=";")

table(res$sample,res$chr)
resf<-res[n.read>50]

resf[,median.depth.sample:=median(n.read),by=c("sample")]
resf[,pct.meth.sample:=sum(meth=="+")/.N,by=c("sample","cpg")]

res_sample<-unique(resf[,.(sample,median.depth.sample,pct.meth.sample,cpg)])
res_sample[,pct.meth.cpgs:=pct.meth.sample[cpg==T],"sample"]
res_sample[,pct.meth.non_cpgs:=pct.meth.sample[cpg==F],'sample']
res_sample<-unique(res_sample[,.(sample,median.depth.sample,pct.meth.cpgs,pct.meth.non_cpgs)])
res_sample[,bisulfite_conversion_efficency:=1-pct.meth.non_cpgs]
res_sample
#     sample median.depth.sample pct.meth.cpgs pct.meth.non_cpgs bisulfite_conversion_efficency
#  1: CBP476                7446    0.04996689        0.02843382                      0.9715662
#  2: CBP481                1684    0.04333183        0.12251671                      0.8774833
#  3: CBP497                5642    0.02108000        0.05109863                      0.9489014
#  4: CBP498                 306    0.05105029        0.10739579                      0.8926042
#  5: CBP514                 564    0.02865078        0.13141411                      0.8685859
#  6: CBP525                1258    0.03361573        0.14180054                      0.8581995
#  7: CBP527                3342    0.04712366        0.12609913                      0.8739009
#  8: CBP529                1420    0.03026797        0.04041762                      0.9595824
#  9: CBP532                1348    0.03027574        0.14142437                      0.8585756
# 10: CBP536               15124    0.03877121        0.05983334                      0.9401667
# 11: CBP543                 320    0.03390937        0.10438615                      0.8956139
# 12: CBP544                 558    0.03931576        0.14035417                      0.8596458
# 13: CBP547                1166    0.06523748        0.17730876                      0.8226912
# 14: CBP554                 448    0.04303419        0.14142173                      0.8585783
fwrite(res_sample,file.path(out,"res_samples_level.csv"),sep=";")

res_pos<-unique(resf[,.(chr,pos,sample,n.meth,n.read,pct.meth,cpg,methyl_type)])
res_pos[,pct.meth.region:=mean(pct.meth),by=c("chr","sample","cpg")]
res_pos[,context:=sapply(methyl_type,function(x)ifelse(str_to_lower(x)=="h","CHH",ifelse(x=="x","CHG","CpG")))]
res_pos[,gene:=str_extract(chr,"[A-Z0-9]+$")]
mtd<-fread("datasets/cd34/bisulfite_seq/metadata_cbps_run1.tsv")
res_pos<-merge(res_pos,mtd,by="sample")
fwrite(res_pos,file.path(out,"res_Cxx_level.csv"),sep=";")


ggplot(unique(res_pos[gene=="SOCS3"&context=='CpG'],by="sample"))+geom_boxplot(aes(x=group,y=pct.meth))+
  facet_wrap("pos")

ggplot(unique(res_pos[gene=="HES1"&context=='CpG'],by="sample"))+geom_boxplot(aes(x=group,y=pct.meth))+
  facet_wrap("pos",scales = "free")


ggplot(unique(res_pos[gene=="SESN2"&context=='CpG'],by="sample"))+geom_boxplot(aes(x=group,y=pct.meth))+
  facet_wrap("pos",scales = "free")

ggplot(unique(res_pos[gene=="E2F2"&context=='CpG'],by="sample"))+geom_boxplot(aes(x=group,y=pct.meth))+
  facet_wrap("pos",scales = "free")

unique(res_pos$gene)

ggplot(unique(res_pos[gene=="ID2"&context=='CpG'],by="sample"))+geom_boxplot(aes(x=group,y=pct.meth))+
  facet_wrap("pos",scales = "free")
