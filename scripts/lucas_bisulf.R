
library(data.table)
library(stringr)
source("scripts/utils/new_utils.R")

out<-"outputs/lucas_cpgs"
dir.create(out)
#get fasta of target regions (non bisulf convert)
#fasta save in 
system("cat ref/lucas_rat/cpg_region_Pnliptp1.fa")
# bisulf convertis fasta ref with bismark
system("tools/bismark/Bismark-0.23.0/bismark_genome_preparation --path_to_aligner /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --verbose ref/lucas_rat") 


#need decrompress fastq
dir.create("outputs/lucas_cpgs/L1-1/fastq/")
system("bzip2 -dkc  ~/RUN/Run_654_rna_methylation/Output/FastQ/L1-1_R1.fastq.bz2 > outputs/lucas_cpgs/L1-1/fastq/L1-1_R1.fastq")
system("bzip2 -dkc  ~/RUN/Run_654_rna_methylation/Output/FastQ/L1-1_R2.fastq.bz2 > outputs/lucas_cpgs/L1-1/fastq/L1-1_R2.fastq")


#bismark tolerating one non-bisulfite mismatch per read

system("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --local --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1  -n 1 -l 20 -1 outputs/lucas_cpgs/L1-1/fastq/L1-1_R1.fastq -2 outputs/lucas_cpgs/L1-1/fastq/L1-1_R2.fastq")

system("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --include_overlap --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/methylation_result outputs/lucas_cpgs/L1-1/L1-1_R1_bismark_bt2_pe.bam")

res<-fread("outputs/lucas_cpgs/L1-1/methylation_result/CpG_OT_L1-1_R1_bismark_bt2_pe.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
head(res,100)
table(res$pos)


#test R1 or R2
system("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --local --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R1  -n 1 -l 20 outputs/lucas_cpgs/L1-1/fastq/L1-1_R1.fastq")
system("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --local --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R2  -n 1 -l 20 outputs/lucas_cpgs/L1-1/fastq/L1-1_R2.fastq")


system("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R1/methylation_result outputs/lucas_cpgs/L1-1/SE_R1/L1-1_R1_bismark_bt2.bam")
system("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R2/methylation_result outputs/lucas_cpgs/L1-1/SE_R2/L1-1_R2_bismark_bt2.bam")

res_r1<-fread("outputs/lucas_cpgs/L1-1/SE_R1/methylation_result/CpG_OT_L1-1_R1_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r1$pos)
#  153  192  200  318  338  354  365  371  441  492  544  576 
# 5475    2 5304 5256    2    2    4    1    3    1  245  233  

res_r1_ctot<-fread("outputs/lucas_cpgs/L1-1/SE_R1/methylation_result/CpG_CTOT_L1-1_R1_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r1_ctot$pos)
 # 153  200  318  332  338  354  365  371  376  406  415  440  441  485  492  513  544  576 
 # 395  236  241    4    8    1   10   13    2   14    6    1    5    1    2    1 1296 1139 

res_r2_ot<-fread("outputs/lucas_cpgs/L1-1/SE_R2/methylation_result/CpG_OT_L1-1_R2_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r2_ot$pos)
 # 108  153  200  318  338  354  365  441  544  576 
 #   1 1698 1528 1431    1    1    1    2   68   69 
res_r2_ctot<-fread("outputs/lucas_cpgs/L1-1/SE_R2/methylation_result/CpG_CTOT_L1-1_R2_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r2_ctot$pos)
 #  88  108  153  192  200  269  318  338  365  370  371  376  381  396  397  404  406  415  440  441  456  484  485  492  513 
 #   4    1  293    3  168    1  231    5   34    3   75    5    4    1    8    1  107   80    2   48   21    3   16   21   13 
 # 522  523  544  549  576 
 #   2   13 4448    6 4426 

#reports
system("tools/bismark/Bismark-0.23.0/bismark2report --dir outputs/lucas_cpgs/L1-1/SE_R1 --alignment_report outputs/lucas_cpgs/L1-1/SE_R1/L1-1_R1_bismark_bt2_SE_report.txt --dedup_report NONE --splitting_report outputs/lucas_cpgs/L1-1/SE_R1_nd/methylation_result/L1-1_R1_bismark_bt2_splitting_report.txt --mbias_report outputs/lucas_cpgs/L1-1/SE_R1_nd/methylation_result/L1-1_R1_bismark_bt2.M-bias.txt")


#SE non directional
system("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R1_nd  -n 1 -l 20 outputs/lucas_cpgs/L1-1/fastq/L1-1_R1.fastq")
system("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R2_nd  -n 1 -l 20 outputs/lucas_cpgs/L1-1/fastq/L1-1_R2.fastq")

system("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R1_nd/methylation_result outputs/lucas_cpgs/L1-1/SE_R1_nd/L1-1_R1_bismark_bt2.bam")
system("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/L1-1/SE_R2_nd/methylation_result outputs/lucas_cpgs/L1-1/SE_R2_nd/L1-1_R2_bismark_bt2.bam")


res_r1_nd<-fread("outputs/lucas_cpgs/L1-1/SE_R1_nd/methylation_result/CpG_OT_L1-1_R1_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r1$pos)
#  206  253  371  407  597  629 
# 5008 5000 5004    1  138  141 

res_r2_nd<-fread("outputs/lucas_cpgs/L1-1/SE_R2_nd/methylation_result/CpG_CTOT_L1-1_R2_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r2$pos)
#  206  253  371  597  629 
# 1370 1311 1249   11   12 
#change rien

#others C
res_r1_nd_nocpg<-fread("outputs/lucas_cpgs/L1-1/SE_R1_nd/methylation_result/Non_CpG_CTOT_L1-1_R1_bismark_bt2.txt",skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))
table(res_r1_nd_nocpg$pos)


#try paired end options 

#for all

samples<-unlist(lapply(paste0("L",1:3,"-"),function(x)paste0(x,1:8)))
for(s in samples){
  print(s)
  #need decrompress fastq
dir.create(file.path("outputs/lucas_cpgs",s,"fastq"),recursive = T)
system(paste0("bzip2 -dkc  ~/RUN/Run_654_rna_methylation/Output/FastQ/",s,"_R1.fastq.bz2 > outputs/lucas_cpgs/",s,"/fastq/",s,"_R1.fastq"))
system(paste0("bzip2 -dkc  ~/RUN/Run_654_rna_methylation/Output/FastQ/",s,"_R2.fastq.bz2 > outputs/lucas_cpgs/",s,"/fastq/",s,"_R2.fastq"))


#bismark tolerating one non-bisulfite mismatch per read


system(paste0("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --local --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/",s,"  -n 1 -l 20 outputs/lucas_cpgs/",s,"/fastq/",s,"_R1.fastq"))
system(paste0("tools/bismark/Bismark-0.23.0/bismark --genome ref/lucas_rat/ --local --non_directional --bowtie2 --path_to_bowtie2 /disks/DATATMP/PhD_AlexandrePelletier/methyl/tools/bowtie2/bowtie2-2.4.4 --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/",s,"  -n 1 -l 20 outputs/lucas_cpgs/",s,"/fastq/",s,"_R2.fastq"))


system(paste0("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/",s," outputs/lucas_cpgs/",s,"/",s,"_R1_bismark_bt2.bam"))
system(paste0("tools/bismark/Bismark-0.23.0/bismark_methylation_extractor --merge_non_CpG --samtools_path /disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12 -o outputs/lucas_cpgs/",s," outputs/lucas_cpgs/",s,"/",s,"_R2_bismark_bt2.bam"))

system(paste0("tools/bismark/Bismark-0.23.0/bismark2report --dir outputs/lucas_cpgs/",s," --alignment_report outputs/lucas_cpgs/",s,"/",s,"_R1_bismark_bt2_SE_report.txt --dedup_report NONE --splitting_report outputs/lucas_cpgs/",s,"/",s,"_R1_bismark_bt2_splitting_report.txt --mbias_report outputs/lucas_cpgs/",s,"/",s,"_R1_bismark_bt2.M-bias.txt"))
system(paste0("tools/bismark/Bismark-0.23.0/bismark2report --dir outputs/lucas_cpgs/",s," --alignment_report outputs/lucas_cpgs/",s,"/",s,"_R2_bismark_bt2_SE_report.txt --dedup_report NONE --splitting_report outputs/lucas_cpgs/",s,"/",s,"_R2_bismark_bt2_splitting_report.txt --mbias_report outputs/lucas_cpgs/",s,"/",s,"_R2_bismark_bt2.M-bias.txt"))


}


res_cpgs_r1<-Reduce(rbind,lapply(samples,function(s)fread(paste0("outputs/lucas_cpgs/",s,"/CpG_OT_",s,"_R1_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R1_OT"][,cpg:=T]))
res_cpgs_r2<-Reduce(rbind,lapply(samples,function(s)fread(paste0("outputs/lucas_cpgs/",s,"/CpG_CTOT_",s,"_R2_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R2_CTOT"][,cpg:=T]))
res_cpgs<-rbind(res_cpgs_r1,res_cpgs_r2)

res_c_r1<-Reduce(rbind,lapply(samples,function(s)fread(paste0("outputs/lucas_cpgs/",s,"/Non_CpG_OT_",s,"_R1_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R1_OT"][,cpg:=F]))
res_c_r2<-Reduce(rbind,lapply(samples,function(s)fread(paste0("outputs/lucas_cpgs/",s,"/Non_CpG_CTOT_",s,"_R2_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R2_CTOT"][,cpg:=F]))
res_cs<-rbind(res_c_r1,res_c_r2)

res<-rbind(res_cpgs,res_cs)
fwrite(res,file.path(out,"res_cpgs_meth_all_reads.csv"),sep=";")

res[,depth:=.N,by=c("pos","sample","cpg")]
res[,depth.glob:=.N,by=c("pos","cpg")]

table(res[methyl_type%in%c("z","Z")&depth.glob>100000]$sample,res[methyl_type%in%c("z","Z")&depth.glob>100000]$pos)

  #       153    200    318    544    576
  # L1-1   5768   5472   5487   4693   4659
  # L1-2  16855  16360  16472  13989  13954
  # L1-3   1849   1423   1437   1538   1562
  # L1-4  29007  28168  28416  22239  22062
  # L1-5  15454  14856  15001  11644  11694
  # L1-6   4807   4582   4663   3790   3763
  # L1-7   5431   5163   5105   4237   4218
  # L1-8 365016 363239 370872 296565 287992
  # L2-1  13367  12997  13110  10759  10625
  # L2-2   9818   9781   9946   8081   8014
  # L2-3   3450   3410   3423   2744   2749
  # L2-4   1643   1590   1588   1272   1237
  # L2-5  12488  12265  12598  10170  10037
  # L2-6  24942  23750  24052  18772  18481
  # L2-7   3288   3260   3298   2802   2768
  # L2-8   3452   3192   3223   2575   2641
  # L3-1   3793   3449   3545   3266   3277
  # L3-2   3855   3786   3855   3212   3202
  # L3-3   1578   1531   1571   1352   1332
  # L3-4     46     32     31     34     28
  # L3-5  14725  14684  14936  12524  12271
  # L3-6   2086   2130   2171   1717   1682
  # L3-7   7224   7268   7377   5828   5791
  # L3-8  17445  17148  17379  13080  13019

resf<-res[depth.glob>100000]
resf[,n.read:=.N,by=c("pos","sample","cpg")]
resf[,n.meth:=sum(meth=="+"),by=c("pos","sample","cpg")]

resf[,pct.meth:=n.meth/n.read,by=c("pos","sample","cpg")]


resf[,median.depth:=median(n.read),by=c("sample")]
resf[,pct.meth.glob:=sum(meth=="+")/.N,by=c("sample","cpg")]

fwrite(resf,file.path(out,"res_pos_10kdepth_glob.csv"),sep=";")

res_sample<-unique(resf[,.(sample,median.depth,pct.meth.glob,cpg)])
res_sample[,pct.meth.cpgs:=pct.meth.glob[cpg==T],"sample"]
res_sample[,pct.meth.non_cpgs:=pct.meth.glob[cpg==F],'sample']
res_sample<-unique(res_sample[,.(sample,median.depth,pct.meth.cpgs,pct.meth.non_cpgs)])
res_sample[,bisulfite_conversion_efficency:=1-pct.meth.non_cpgs]
res_sample
fwrite(res_sample,file.path(out,"res_samples_level.csv"),sep=";")

res_pos<-unique(resf[,.(pos,sample,n.meth,n.read,pct.meth,cpg,methyl_type)])
res_pos[,chr:="chr1"]
res_pos[,pos_fasta:=pos]
res_pos[,pos:=279882093+pos]
res_pos[,context:=sapply(methyl_type,function(x)ifelse(str_to_lower(x)=="h","CHH",ifelse(x=="x","CHG","CpG")))]
res_pos<-res_pos[,-"methyl_type"]
res_pos<-res_pos[order(-cpg,pos_fasta)][,.(chr,pos,sample,n.meth,n.read,pct.meth,cpg,context,pos_fasta)]

fwrite(res_pos,file.path(out,"res_Cxx_level.csv"),sep=";")

res_pos<-fread(file.path(out,"res_Cxx_level.csv"),sep=";")
res_pos<-unique(res_pos)

ggplot(res_pos[pos_fasta==544])+geom_col(aes(x=sample,y=pct.meth))


#focus on 544 576 r2
res_cpgs_r2<-Reduce(rbind,lapply(samples,function(s)fread(paste0("outputs/lucas_cpgs/",s,"/CpG_CTOT_",s,"_R2_bismark_bt2.txt"),skip = 1,col.names = c("seq_id","meth","chr","pos","methyl_type"))[,sample:=s][,read:="R2_CTOT"][,cpg:=T]))

res_cpgs_r2<-res_cpgs_r2[pos%in%c(544,576)]
res_cpgs_r2[,n.read:=.N,by=c("pos","sample","cpg")]
res_cpgs_r2[,n.meth:=sum(meth=="+"),by=c("pos","sample","cpg")]

res_cpgs_r2[,pct.meth:=n.meth/n.read,by=c("pos","sample","cpg")]
res_pos_int<-unique(res_cpgs_r2[,.(pos,sample,n.meth,n.read,pct.meth)])
ggplot(res_pos_int)+geom_col(aes(x=sample,y=pct.meth))+facet_wrap("pos")
