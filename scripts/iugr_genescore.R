#IUGR Genescore

source("scripts/utils/new_utils.R")
library(limma)

out<-"outputs/iugr_genescore"
dir.create(out)


methf<-fread("datasets/cd34/meth_data_filtered.csv.gz")

mtd<-mtd[group%in%c("CTRL","IUGR")]

#pca analysis
bool_vars<-c("latino","preterm","GDM","drugs","etoh", "smoking")
categorical_vars<-c(bool_vars,"group","sex","group_sex","ethnicity","lab","batch","date","DNA.extraction","sequencing","year")

mtd[,(categorical_vars):=lapply(.SD,as.factor),.SDcols=categorical_vars]

lapply(mtd[,.SD,.SDcols=categorical_vars],function(x)sum(table(x)==1)) #exclude sequencing and date

mtd<-mtd[,-c("sequencing","date","drugs","ethnicity","year")]
library(ggrepel)
meth_mat<-as.matrix(data.frame(methf,row.names = "cpg_id")[,mtd$sample])
pca<-prcomp(t(meth_mat))
saveRDS(pca,fp(out,"pca_iugrctrl.rds"))

pc_mtd<-merge(mtd,data.table(pca$x,keep.rownames = "sample"))


#model
vars_to_include<-c("batch","mat.age","group_complexity_fac","group_sex","latino","PC2")
mtd_f<-pc_mtd[,to_keep:=rowSums(is.na(.SD))==0,.SDcols=vars_to_include][to_keep==T] #5 samples remove
fwrite(mtd_f,fp(out,"pca_mtd_iugrctrl.csv"),sep=";")

nrow(mtd_f) #72
table(mtd_f$group,mtd_f$batch)
  #       1  2
  # CTRL 18 16
  # IUGR 19 19
  # LGA   0  0

table(mtd_f$group,mtd_f$sex)
  #       F  M
  # CTRL 16 18
  # IUGR 20 18
  # LGA   0  0
formule<- ~0 + group_sex   + batch+ group_complexity_fac +mat.age  + latino + PC2

mtd_f[,group_sex:=factor(group_sex,levels = unique(mtd_f$group_sex))]

design<-model.matrix(formule,data = data.frame(mtd_f,row.names = "sample"))
fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_f$sample], design)

cont.matrix <- makeContrasts(C.I = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexIUGR_F+group_sexIUGR_M)",
                             levels=design)

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

res<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")
fwrite(res,fp(out,"res_limma.tsv.gz"),sep="\t")

res<-fread(fp(out,"res_limma.tsv.gz"),sep="\t",
           select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
           col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))
res[pval<0.001]

table(res[pval<0.001]) #707 DMCs
ggplot(res)+geom_point(aes(x=meth.change,y=-log10(pval),col=pval<0.001&abs(meth.change)>25),size=0.5)+
  scale_color_manual(values = c("grey","red"))

#by batch
table(mtd_f$batch)
#  1  2 
# 37 35 
vars_to_include<-c("mat.age","group_complexity_fac","group_sex","latino","PC2")


res_batch<-Reduce(rbind,lapply(1:2,function(b){
  mtd_fc<-mtd_f[batch==b]
  formule<- ~0 + group_sex  + group_complexity_fac +mat.age  + latino + PC2
  design<-model.matrix(formule,data = data.frame(mtd_fc,row.names = "sample"))
  fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_fc$sample], design)
  cont.matrix <- makeContrasts(C.I = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexIUGR_F+group_sexIUGR_M)",
                             CF.IF="group_sexCTRL_F-group_sexIUGR_F",
                             CM.IM="group_sexCTRL_M-group_sexIUGR_M",
                             levels=design)
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  res<-Reduce(rbind,lapply(colnames(cont.matrix), function(comp)data.table(topTable(fit2,coef = comp,n = Inf),keep.rownames = "cpg_id")[,compa:=comp]))
  return(res[,batch:=b])
}))

fwrite(res_batch,fp(out,"res_limma_cohorts.tsv.gz"),sep="\t")
res_batch<-fread(fp(out,"res_limma_cohorts.tsv.gz"),sep="\t")
res_batchf<-res_batch[compa=="C.I"]
table(res_batchf[adj.P.Val<=0.1,.(compa,batch)])
#      batch
# compa  2
#   C.I 23
ggplot(res_batchf)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>25),size=0.1)+
  facet_wrap("batch")+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()
ggsave(fp(out,"fig1A.volcanos_limma_cohorts_C.I.png"),width = 12,height = 5)


ggplot(res_batch[compa!="C.I"])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>25),size=0.1)+
  facet_wrap(compa~batch)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()

#hypomethylation, really ?
methf_mtd<-merge(melt(methf[,.SD,.SDcols=c("cpg_id",mtd_f$sample)],id.vars = "cpg_id",value.name = "methylation",variable.name = "sample"))

#genescore calc
res<-fread(fp(out,"res_limma.tsv.gz"),sep="\t",
           select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
           col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))

cpgs_score<-fread("outputs/02-gene_score_calculation_and_validation/cpgs_genes_annot_and_weight.csv.gz")



res_anno<-merge(res,cpgs_score,by="cpg_id")

res_anno[,cpg_score:=-log10(pval)*meth.change*links_weight*regul_weight] #divided by n_sample ro normalized gene score ~ n_sampleres_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/5),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
#res_anno[,n.cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/5),by="gene"] #n.cpg_weight to reduce the influence of the n_cpg by gene to the GeneScore
#res_anno[,gene_score:=sum(cpg_score)*n.cpg_weight,by="gene"] 

# +sep prom/enh
res_anno[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other")]
res_anno[is.na(tss_dist),region_type:="other"]

res_anno[,n_cpg_weight_region:=(1/sum(1/(abs(cpg_score)+1)))^(1/3.8),by=c('region_type',"gene")]
res_anno[region_type=="promoter",gene_score_region:=sum(cpg_score)*n_cpg_weight_region,by=c("gene")]
res_anno[region_type=="other",gene_score_region:=sum(abs(cpg_score))*n_cpg_weight_region,by=c("gene")]

res_anno[,gene_score_add:=sum(unique(abs(gene_score_region)),na.rm = T),by="gene"]
meth_metrics<-c(meth_metrics,"gene_score_add")
