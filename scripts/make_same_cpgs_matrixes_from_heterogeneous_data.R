#make CpG matrixes Placenta, Pancreas (and put liver in same format if necessary)
library(data.table)
library(stringr)
#placenta : 
meth<-fread("datasets/data_placenta/run44_NICHD_ctrl-bkg.txt")
dim(meth)
head(meth[,1:10])

cols<-c("TargetID",colnames(meth)[str_detect(colnames(meth),"AVG_Beta")])

meth<-meth[,.SD,.SDcols=cols]

meth2<-fread("datasets/data_placenta/run47_NICHD_ctrl_bkg.txt")
cols2<-c("TargetID",colnames(meth2)[str_detect(colnames(meth2),"AVG_Beta")])
meth2<-meth2[,.SD,.SDcols=cols]
samples1<-str_remove(cols[-1],"\\.AVG_Beta")
samples2<-str_remove(cols2[-1],"\\.AVG_Beta")

mtd<-data.table(sample=c(samples1,samples2))
mtd[sample%in%samples1,batch:=1]
mtd[sample%in%samples2,batch:=2]

meth<-merge(meth,meth2,by="TargetID")
meth[,1:10]
ncol(meth)
fwrite(meth,"datasets/data_placenta/312_samples_beta_value_matrix.csv",sep=";")

expr<-fread("datasets/data_placenta/expression_84samples_Placenta_031715.txt")
