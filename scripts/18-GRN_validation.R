out<-"outputs/18-GRN_validation"
dir.create(out)
source("scripts/utils/new_utils.R")

#validate regulons/GRN with external data
#=>integr NicheNet in our regulons analysis : overlap regulons tf-targets matrix ? wich signalling upstream of this stimulated asso regulons ?
#data : https://zenodo.org/record/3260758
#data sources : https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0667-5/MediaObjects/41592_2019_667_MOESM3_ESM.xlsx
system("wget -O outputs/18-GRN_validation/gr_network.rds https://zenodo.org/record/3260758/files/gr_network.rds?download=1")

system("wget -O outputs/18-GRN_validation/tf_target_matrix.rds https://zenodo.org/record/3260758/files/tf_target_matrix.rds?download=1")


grn<-readRDS("outputs/18-GRN_validation/gr_network.rds")

head(grn)
grn<-data.table(grn)
grn #3592299 interactions

unique(grn,by=c("from","to")) #3M non redundant interactions


tfg_hscp<-fread("outputs/16-GRN_final/tf_target_interactions.csv")
intersect(tfg[tf=="KLF2"]$target,grn[from=="KLF2"]$to)

intersect(tfg[tf=="KLF4"]$target,grn[from=="KLF4"&database=='harmonizome_gr']$to)

tft_mat<-readRDS("outputs/18-GRN_validation/tf_target_matrix.rds")
head(tft_mat[,1:10])
plot(density(tft_mat))
summary(as.vector(tft_mat)) #max 1.2798671 
dim(tft_mat) #25345  4486

tft_scores<-melt(data.table(tft_mat[,intersect(unique(tfg$tf),colnames(tft_mat))],keep.rownames = "TF"),variable.name = "target",value.name = "score")
tft_scores[,in_regulon_hsc_peak_filtered:=target%in%tfg_hscp[tf==TF]$target,by="TF"]

tft_scores[,pct_regul_hsc_peak_filtered:=sum(score>0.01)/.N,.(TF,in_regulon_hsc_peak_filtered)]

tft_scores[,pct_regul_hsc_peak_filtered:=sum(score>0.01)/.N,.(TF,in_regulon_hsc_peak_filtered)]

ggplot(unique(tft_scores[TF%in%c("EGR1","KLF2","KLF4")],by=c("TF","in_regulon_hsc_peak_filtered")))+geom_col(aes(y=pct_regul_hsc_peak_filtered,x=in_regulon_hsc_peak_filtered))+facet_wrap("TF")

p1<-ggplot(unique(tft_scores,by=c("TF","in_regulon_hsc_peak_filtered")))+geom_boxplot(aes(y=pct_regul_hsc_peak_filtered,x=in_regulon_hsc_peak_filtered))


#with non filtered regulons
tfg<-fread("outputs/13-GRN_integr/tf_target_interactions.csv")
tft_scores[,in_regulon:=target%in%tfg[tf==TF]$target,by="TF"]

tft_scores[,pct_regul:=sum(score>0)/.N,.(TF,in_regulon)]


ggplot(unique(tft_scores[TF%in%c("EGR1","KLF2","KLF4")],by=c("TF","in_regulon")))+geom_col(aes(y=pct_regul,x=in_regulon))+facet_wrap("TF")

ggplot(unique(tft_scores,by=c("TF","in_regulon")))+geom_boxplot(aes(y=pct_regul,x=in_regulon))



tft_scores[,mean_regul:=mean(score),.(TF,in_regulon)]

ggplot(unique(tft_scores,by=c("TF","in_regulon")))+geom_boxplot(aes(y=mean_regul,x=in_regulon))


tft_scores
