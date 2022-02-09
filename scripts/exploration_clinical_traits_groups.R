
#Correlation Methylation avec poids de naissance

source("scripts/utils/new_utils.R")
meth<-fread("datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
mtd<-fread("datasets/cd34/cleaned_batch_CD34_and_PCs_of_modeled_methyl_data.csv")

pca<-prcomp(t(as.matrix(meth,rownames = 1)))
summary(pca)

pct_var<-pca$sdev^2/sum(pca$sdev^2)
names(pct_var)<-paste0("PC",1:length(pct_var))
round(pct_var*100)
fwrite(data.table(PC=paste0("PC",1:length(pct_var)),var_explained=pct_var),
       "datasets/cd34/var_explained_pca_modeled_data.csv",sep=";")

#
plotPCVarExplain(pca,1:50)
pvals<-CorrelCovarPCs(pca = pca,rngPCs = 1:21,batch = data.frame(mtd,row.names = "sample"),
             var_num=c("Weight.at.term..lbs.","Weight..g.","Length..cm.","HC..cm.","PI","GA..wk.","Mat.Age"),
             var_fac=c("Sex","Group_name","Group_Sex"),res="pval",return=T,plot=F)

r2s<-CorrelCovarPCs(pca = pca,rngPCs = 1:21,batch = data.frame(mtd,row.names = "sample"),
             var_num=c("Weight.at.term..lbs.","Weight..g.","Length..cm.","HC..cm.","PI","GA..wk.","Mat.Age"),
             var_fac=c("Sex","Group_name","Group_Sex"),
             res="r2",plot=F,return = T)

r2s[pvals>0.05]<-0

pheatmap(r2s,cluster_rows = F,cluster_cols = F,
           labels_col= paste0(colnames(r2s),"(",round(pct_var[colnames(r2s)]*100),"%)"),
           display_numbers = T,
           color = colorRampPalette(c("white", "red"))(10))

dt<-data.table(clinical_var=rownames(r2s),cumulative_r2=rowSums(r2s))

ggplot(dt)+
  geom_col(aes(x=clinical_var,y=cumulative_r2,fill=clinical_var))+scale_x_discrete(limits=dt[order(-cumulative_r2)]$clinical_var)

mtd[sample=="CBP041",.(PC1)]
mtd[,Group_Sex:=factor(Group_Sex,levels = c("I_F","I_M","C_F","C_M","L_F","L_M"))]
mtd[,Group_name:=factor(Group_name,levels = c("I","C","L"))]

p1<-ggplot(mtd)+geom_boxplot(aes(x=Group_name,y=Weight..g.,fill=Gender))
p2<-ggplot(mtd)+geom_boxplot(aes(x=Group_name,y=Length..cm.,fill=Gender))
p3<-ggplot(mtd)+geom_boxplot(aes(x=Group_name,y=HC..cm.,fill=Gender))
p4<-ggplot(mtd)+geom_boxplot(aes(x=Group_name,y=PI,fill=Gender))

(p1|p2)/(p3|p4)+plot_layout(guides = 'collect')

ggplot(mtd)+geom_boxplot(aes(x=Group_Sex,y=Weight..g.,fill=Group_Sex),alpha=0.5,outlier.shape = NA)+
  geom_jitter(aes(x=Group_Sex,y=Weight..g.,col=Group_Sex))

ggplot(mtd)+geom_boxplot(aes(x=Group_Sex,y=HC..cm.,fill=Group_Sex),alpha=0.5,outlier.shape = NA)+
  geom_jitter(aes(x=Group_Sex,y=HC..cm.,col=Group_Sex))

ggplot(mtd)+geom_boxplot(aes(x=Group_Sex,y=Length..cm.,fill=Group_Sex),alpha=0.5,outlier.shape = NA)+
  geom_jitter(aes(x=Group_Sex,y=Length..cm.,col=Group_Sex))

ggplot(mtd)+geom_boxplot(aes(x=Group_Sex,y=PI,fill=Group_Sex),alpha=0.5,outlier.shape = NA)+
  geom_jitter(aes(x=Group_Sex,y=PI,col=Group_Sex))


ggplot(mtd)+geom_boxplot(aes(x=Group_Sex,y=GA..wk.,fill=Group_Sex),alpha=0.5,outlier.shape = NA)+
  geom_jitter(aes(x=Group_Sex,y=GA..wk.,col=Group_Sex))


mtd[Group_name=="I"&Weight..g.>4000] #CBP 449 bad attribution

mtd[sample=="CBP449",Group:=2][sample=="CBP449",Group_name:="L"][sample=="CBP449",Group_Sex:="L_M"]
fwrite(mtd,"datasets/cd34/cleaned_batch_CD34_and_PCs_of_modeled_methyl_data.csv")
fwrite(mtd[,.SD,.SDcols=!str_detect(colnames(mtd),"PC")],"datasets/cd34/cleaned_batch_CD34_250121.csv",sep=";")
