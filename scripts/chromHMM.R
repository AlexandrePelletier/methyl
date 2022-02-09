
#[in BASH]

# I) CHROMHMM directly with reads bed files ChIPseq data
# [in methyl/ref/Chromatin_Annot/Liver/]
#need first binarize bed :

java -Xmx4000M -jar ../ChromHMM/ChromHMM.jar BinarizeBed ../ChromHMM/CHROMSIZES/hg19.txt ChIP_Liver/ cellmarkfiletable_GSMreads.txt ChromHMM_output_GSMreads/Binarized_Bed

#with cellmarkfiletable_GSMreads.txt : 
# liver	H3K27ac	GSM1112809_BI.Adult_Liver.H3K27ac.4.bed	GSM669910_BI.Adult_Liver.Input.5.bed.gz
# liver	H3K27me3	GSM1112814_BI.Adult_Liver.H3K27me3.4.bed	GSM669910_BI.Adult_Liver.Input.5.bed.gz
# liver	H3K36me3	GSM537708_BI.Adult_Liver.H3K36me3.5.bed	GSM669910_BI.Adult_Liver.Input.5.bed.gz
# liver	H3K4me1	GSM621654_BI.Adult_Liver.H3K4me1.5.bed	GSM669910_BI.Adult_Liver.Input.5.bed.gz
# liver	H3K4me3	GSM621675_BI.Adult_Liver.H3K4me3.4.bed	GSM669910_BI.Adult_Liver.Input.5.bed.gz
# liver	H3K9me3	GSM669986_BI.Adult_Liver.H3K9me3.4.bed	GSM669910_BI.Adult_Liver.Input.5.bed.gz

#then learn model : 
java -Xmx4000M -Djava.awt.headless=true -jar ../ChromHMM/ChromHMM.jar LearnModel ChromHMM_output_GSMreads/Binarized_Bed ChromHMM_output_GSMreads/ 6 hg19 


#### II) call peaks before from chip seq data 
# [in methyl/ref/Chromatin_Annot/Liver/ChIP_Liver]

macs2 callpeak -t GSM1112809_BI.Adult_Liver.H3K27ac.4.bed -c GSM669910_BI.Adult_Liver.Input.5.bed.gz --format=BED --name Liver_H3K27ac_peaks.bed


macs2 callpeak -t GSM621654_BI.Adult_Liver.H3K4me1.5.bed -c GSM669910_BI.Adult_Liver.Input.5.bed.gz --format=BED --name Liver_H3K4me1_peaks.bed

macs2 callpeak -t GSM1112814_BI.Adult_Liver.H3K27me3.4.bed -c GSM669910_BI.Adult_Liver.Input.5.bed.gz --format=BED --name Liver_H3K27me3_peaks.bed

macs2 callpeak -t GSM669986_BI.Adult_Liver.H3K9me3.4.bed -c GSM669910_BI.Adult_Liver.Input.5.bed.gz --format=BED --name Liver_H3K9me3_peaks.bed #doesnt work  (doesnt build model*)

macs2 callpeak -t GSM537708_BI.Adult_Liver.H3K36me3.5.bed -c GSM669910_BI.Adult_Liver.Input.5.bed.gz --format=BED --name Liver_H3K36me3_peaks.bed

macs2 callpeak -t GSM621675_BI.Adult_Liver.H3K4me3.4.bed -c GSM669910_BI.Adult_Liver.Input.5.bed.gz --format=BED --name Liver_H3K4me3_peaks.bed #doesnt work too*


#*Too few paired peaks (0) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead.

#### sort

module load bedtools2/2.24.0/gcc.4.4.7


bedtools sort -i Liver_H3K36me3_peaks.bed_peaks.narrowPeak > Liver_H3K36me3_peaks.bed_peaks.narrowPeak_sorted

bedtools sort -i Liver_H3K4me3_peaks.bed_peaks.narrowPeak > Liver_H3K4me3_peaks.bed_peaks.narrowPeak_sorted

bedtools sort -i Liver_H3K4me1_peaks.bed_peaks.narrowPeak > Liver_H3K4me1_peaks.bed_peaks.narrowPeak_sorted

bedtools sort -i Liver_H3K27me3_peaks.bed_peaks.narrowPeak > Liver_H3K27me3_peaks.bed_peaks.narrowPeak_sorted

bedtools sort -i Liver_H3K27ac_peaks.bed_peaks.narrowPeak > Liver_H3K27ac_peaks.bed_peaks.narrowPeak_sorted
#[end in methyl/ref/Chromatin_Annot/Liver/ChIP_Liver]

#[in methyl/ref/Chromatin_Annot]
#make windows across the genome and sort him (bedwindow)

bedtools makewindows -g ChromHMM/CHROMSIZES/hg19.txt -w 100 > hg19_100bpWin.bed
bedtools sort -i hg19_100bpWin.bed > hg19_100bpWin_sorted.bed

#map peaks into 100bp windows
bedtools map -a hg19_100bpWin_sorted.bed -b Liver/ChIP_Liver/Liver_H3K27ac_peaks.bed_peaks.narrowPeak -c 5 -o mean -null "0" > Liver/H3K27ac_100.bed

bedtools map -a hg19_100bpWin_sorted.bed -b Liver/ChIP_Liver/Liver_H3K27me3_peaks.bed_peaks.narrowPeak -c 5 -o mean -null "0" > Liver/H3K27me3_100.bed

bedtools map -a hg19_100bpWin_sorted.bed -b Liver/ChIP_Liver/Liver_H3K4me1_peaks.bed_peaks.narrowPeak -c 5 -o mean -null "0" > Liver/H3K4me1_100.bed

bedtools map -a hg19_100bpWin_sorted.bed -b Liver/ChIP_Liver/Liver_H3K4me3_peaks.bed_peaks.narrowPeak -c 5 -o mean -null "0" > Liver/H3K4me3_100.bed

bedtools map -a hg19_100bpWin_sorted.bed -b Liver/ChIP_Liver/Liver_H3K36me3_peaks.bed_peaks.narrowPeak -c 5 -o mean -null "0" > Liver/H3K36me3_100.bed


java -Xmx4000M -jar ~/Desktop/ChromHMM/ChromHMM.jar LearnModel  ~/Desktop/ChromHMMbed/window/  ~/Desktop/ChromHMMbed/window/output 7 hg19


java -Xmx4000M -jar ~/Desktop/ChromHMM/ChromHMM.jar LearnModel  ~/Desktop/ChromHMMbed/test/  ~/Desktop/ChromHMMbed/test/output 7 hg19


java -Xmx4000M -jar ~/Desktop/ChromHMM/ChromHMM.jar LearnModel -init random ~/Desktop/ChromHMMbed/callpeak/  ~/Desktop/ChromHMMbed/callpeak/output2 7 hg19


#ChromHMM with this win
java -Xmx4000M -jar ../ChromHMM/ChromHMM.jar BinarizeBed -center ../ChromHMM/CHROMSIZES/hg19.txt . cellmarkfiletable_100bpWin.txt ChromHMM_output_100bpWin/Binarized_Bed

java -Xmx4000M -Djava.awt.headless=true -jar ../ChromHMM/ChromHMM.jar LearnModel -init random ChromHMM_output_100bpWin/Binarized_Bed ChromHMM_output_100bpWin/ 5 hg19 
#output : doesn' work

#[end in BASH]
####plot chromHMM


options(stringAsFactors=F)

FEATURE=read.table("~/Desktop/chromHMMbed/test/placenta_7_segments.bed")
distance=read.table("/R/list/distanceTSS_020216.txt")

bedTools.2in<-function(functionstring="/usr/local/bin/intersectBed",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}


inter=bedTools.2in(bed1=FEATURE,bed2=distance,opt.string="-wa -wb")

f1=subset(inter,inter[,4]=="E1")
f2=subset(inter,inter[,4]=="E2")
f3=subset(inter,inter[,4]=="E3")
f4=subset(inter,inter[,4]=="E4")
f5=subset(inter,inter[,4]=="E5")
f6=subset(inter,inter[,4]=="E6")
f7=subset(inter,inter[,4]=="E7")


plot.multi.dens <- function(s)
{
junk.x = NULL
junk.y = NULL
for(i in 1:length(s))
{
junk.x = c(junk.x, density(s[[i]])$x)
junk.y = c(junk.y, density(s[[i]])$y)
}
xr <- range(junk.x)
yr <- range(junk.y)
plot(density(s[[1]]), xlim = xr, ylim = yr, main = "",)
for(i in 1:length(s))
{
lines(density(s[[i]]), xlim = xr, ylim = yr, col = i, lwd=1)

}
}
plot.multi.dens(list(f1[,11],f2[,11],f3[,11],f4[,11],f5[,11],f6[,11],f7[,11]))

legend(x=50, y=0.08,legend=c("feature1","feature2","feature3","feature4","feature5","feature6","feature7"), col=(1:7), lwd=2, lty = 1)

cpg=read.table("/placenta/refILLUMINA450K.txt",sep="\t")


intercpg=bedTools.2in(bed1=FEATURE,bed2=cpg[,1:3],opt.string="-wa -wb")

length(which(intercpg[,4]=="E1"))
length(which(intercpg[,4]=="E2"))
length(which(intercpg[,4]=="E3"))
length(which(intercpg[,4]=="E4"))
length(which(intercpg[,4]=="E5"))
length(which(intercpg[,4]=="E6"))
length(which(intercpg[,4]=="E7"))

length(which(FEATURE[,4]=="E1"))
length(which(FEATURE[,4]=="E2"))
length(which(FEATURE[,4]=="E3"))
length(which(FEATURE[,4]=="E4"))
length(which(FEATURE[,4]=="E5"))
length(which(FEATURE[,4]=="E6"))
length(which(FEATURE[,4]=="E7"))



### I have to look by bp 

intercpg=bedTools.2in(bed1=FEATURE,bed2=cpg[,1:3],opt.string="-wo")

s1=subset(intercpg,intercpg[,4]=="E1")
s2=subset(intercpg,intercpg[,4]=="E2")
s3=subset(intercpg,intercpg[,4]=="E3")
s4=subset(intercpg,intercpg[,4]=="E4")
s5=subset(intercpg,intercpg[,4]=="E5")
s6=subset(intercpg,intercpg[,4]=="E6")
s7=subset(intercpg,intercpg[,4]=="E7")

sum(s1[,8])
sum(s2[,8])
sum(s3[,8])
sum(s4[,8])
sum(s5[,8])
sum(s6[,8])
sum(s7[,8])

FEATURE=data.frame(FEATURE,FEATURE[,3]-FEATURE[,2])

s1=subset(FEATURE,FEATURE[,4]=="E1")
s2=subset(FEATURE,FEATURE[,4]=="E2")
s3=subset(FEATURE,FEATURE[,4]=="E3")
s4=subset(FEATURE,FEATURE[,4]=="E4")
s5=subset(FEATURE,FEATURE[,4]=="E5")
s6=subset(FEATURE,FEATURE[,4]=="E6")
s7=subset(FEATURE,FEATURE[,4]=="E7")

sum(as.numeric(s1[,5]))
sum(as.numeric(s2[,5]))
sum(as.numeric(s3[,5]))
sum(as.numeric(s4[,5]))
sum(as.numeric(s5[,5]))
sum(as.numeric(s6[,5]))
sum(as.numeric(s7[,5]))






bi1=read.table("~/Desktop/chromHMMbed/window/placenta_chr1_binary.txt",skip=1,header=T)
bi2=read.table("~/Desktop/chromHMMbed/window/placenta_chr2_binary.txt",skip=1,header=T)
bi3=read.table("~/Desktop/chromHMMbed/window/placenta_chr3_binary.txt",skip=1,header=T)
bi4=read.table("~/Desktop/chromHMMbed/window/placenta_chr4_binary.txt",skip=1,header=T)
bi5=read.table("~/Desktop/chromHMMbed/window/placenta_chr5_binary.txt",skip=1,header=T)

bi6=read.table("~/Desktop/chromHMMbed/window/placenta_chr6_binary.txt",skip=1,header=T)
bi7=read.table("~/Desktop/chromHMMbed/window/placenta_chr7_binary.txt",skip=1,header=T)
bi8=read.table("~/Desktop/chromHMMbed/window/placenta_chr8_binary.txt",skip=1,header=T)
bi9=read.table("~/Desktop/chromHMMbed/window/placenta_chr9_binary.txt",skip=1,header=T)
bi10=read.table("~/Desktop/chromHMMbed/window/placenta_chr10_binary.txt",skip=1,header=T)
bi11=read.table("~/Desktop/chromHMMbed/window/placenta_chr11_binary.txt",skip=1,header=T)
bi12=read.table("~/Desktop/chromHMMbed/window/placenta_chr12_binary.txt",skip=1,header=T)
bi13=read.table("~/Desktop/chromHMMbed/window/placenta_chr13_binary.txt",skip=1,header=T)
bi14=read.table("~/Desktop/chromHMMbed/window/placenta_chr14_binary.txt",skip=1,header=T)
bi15=read.table("~/Desktop/chromHMMbed/window/placenta_chr15_binary.txt",skip=1,header=T)
bi16=read.table("~/Desktop/chromHMMbed/window/placenta_chr16_binary.txt",skip=1,header=T)
bi17=read.table("~/Desktop/chromHMMbed/window/placenta_chr17_binary.txt",skip=1,header=T)
bi18=read.table("~/Desktop/chromHMMbed/window/placenta_chr18_binary.txt",skip=1,header=T)
bi19=read.table("~/Desktop/chromHMMbed/window/placenta_chr19_binary.txt",skip=1,header=T)
bi20=read.table("~/Desktop/chromHMMbed/window/placenta_chr20_binary.txt",skip=1,header=T)
bi21=read.table("~/Desktop/chromHMMbed/window/placenta_chr21_binary.txt",skip=1,header=T)
bi22=read.table("~/Desktop/chromHMMbed/window/placenta_chr22_binary.txt",skip=1,header=T)
biX=read.table("~/Desktop/chromHMMbed/window/placenta_chrX_binary.txt",skip=1,header=T)
biY=read.table("~/Desktop/chromHMMbed/window/placenta_chrY_binary.txt",skip=1,header=T)
biM=read.table("~/Desktop/chromHMMbed/window/placenta_chrM_binary.txt",skip=1,header=T)

summary(bi1)
summary(bi2)
summary(bi3)
summary(bi4)
summary(bi5)
summary(bi6)
summary(bi7)
summary(bi8)
summary(bi9)
summary(bi10)
summary(bi11)
summary(bi12)
summary(bi13)
summary(bi14)
summary(bi15)
summary(bi16)
summary(bi17)
summary(bi18)
summary(bi19)
summary(bi20)
summary(bi21)
summary(bi22)
summary(biX)
summary(biY)
summary(biM)



