source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(limma)
out<-"outputs/01-lga_vs_ctrl_limma_DMCs_analysis"
dir.create(out)


meth_file<-here("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt")
mtd_file<-here("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")


#Loading====

mtd<-fread(mtd_file,
            select = c("sample","Group_name","Gender","Group_Sex",
                       "Weight..g.","Weight.at.term..lbs.","GA..wk.","Length..cm.","Mat.Age","HC..cm.","PI","SeqDepth",
                       "latino","Preterm","GDM","Drugs","Etoh", "Smoking",
                       "Race","Labor","batch","date","DNA.extraction","sequencing",
                       "Library_Complexity","Group_Complexity","Group_Complexity_Fac","GroupBatch_Complexity","GroupBatch_Complexity_Fac"),
            col.names = c("sample","group","sex","group_sex",
                       "weight.g","weight_at_term.lbs","gest.age.wk","length.cm","mat.age","head_circumf.cm.","ponderal_index","seq.depth",
                       "latino","preterm","GDM","drugs","etoh", "smoking",
                       "ethnicity","lab","batch","date","DNA.extraction","sequencing",
                       "library_complexity","group_complexity","group_complexity_fac","groupbatch_complexity","groupbatch_complexity_fac"))

mtd[,group:=ifelse(group=="L","LGA",ifelse(group=="I","IUGR","CTRL"))]
mtd[,group_sex:=paste(group,sex,sep = "_")]
mtd[,year:=str_extract(date,"20[0-9]+$")]
mtd[ethnicity%in%c("Declined",""," "),ethnicity:=NA]
mtd[,ethnicity:=factor(ethnicity,levels=unique(ethnicity))]
bool_vars<-c("latino","preterm","GDM","drugs","etoh", "smoking")
mtd[,(bool_vars):=lapply(.SD,as.logical),.SDcols=bool_vars]

categorical_vars<-c(bool_vars,"group","sex","group_sex","ethnicity","lab","batch","date","DNA.extraction","sequencing","year")
mtd[,(categorical_vars):=lapply(.SD,as.factor),.SDcols=categorical_vars]

numerical_vars<-c("weight.g","weight_at_term.lbs","gest.age.wk","length.cm","mat.age","head_circumf.cm.","ponderal_index","seq.depth","library_complexity","group_complexity","group_complexity_fac","groupbatch_complexity","groupbatch_complexity_fac")

fwrite(mtd,"datasets/cd34/metadata_formated.csv")
meth<-fread(meth_file,
            select = c("id","chr","start",mtd$sample,"msp1c","confidenceScore"),
            col.names = c("cpg_id","chr","pos",mtd$sample,"msp1c","confidence_score"))
nrow(meth) #1709224

all(mtd$sample%in%colnames(meth))
nrow(mtd) #119
table(mtd[,.(group,sex)])
#       sex
# group   F  M
#   CTRL 17 20
#   IUGR 20 20
#   LGA  21 21

#data exploration
plot(density(na.omit(as.matrix(meth[,.SD,.SDcols=mtd$sample]))))
abline(v=10)


#CpG filtering====
# add useful quality metrics 

meth[,n.na:=rowSums(is.na(.SD)),.SDcols=mtd$sample]
meth[,n.not.methyl:=rowSums(.SD>10,na.rm = T),.SDcols=mtd$sample]
meth[,pct.zero:=rowSums(.SD==0,na.rm = T)/length(mtd$sample),.SDcols=mtd$sample]
meth[,n.methyl.not.zero:=rowSums(.SD>0&.SD<10,na.rm = T),.SDcols=mtd$sample]
# + msp1count (msp1c) and confidencescore, already present in the data
#   msp1c : count for reference library digest with msp1
#   confidencescore : sum of all count for hpa2 libraries and msp1 library, normalized by library size

#   see 01A-CpG_filtering to see how we determine the threshold,


quantile(meth$msp1c,0.125,na.rm = T)*76541158 #total number of msp1 count
methf<-meth[msp1c>quantile(msp1c,0.125,na.rm=T)&n.na==0]

methf<-methf[n.not.methyl>4]

quantile(meth$confidence_score,0.2,na.rm = T) *76541158

methf<-methf[confidence_score>quantile(confidence_score,0.2,na.rm=T)]
methf<-methf[!(pct.zero>0.7&n.methyl.not.zero==0)]

methf #754931 CpGs after filtering

fwrite(methf,"datasets/cd34/meth_data_filtered.csv.gz")


#FOCUS ON CTRL LGA====
methf<-fread("datasets/cd34/meth_data_filtered.csv.gz")
mtd<-mtd[group%in%c("CTRL","LGA")]
fwrite(mtd,"datasets/cd34/metadata_cl_190421.csv",sep=";")
mtd<-fread("datasets/cd34/metadata_cl_190421.csv",sep=";")


all(mtd$sample%in%colnames(methf))
nrow(mtd) #79
table(mtd[,.(group,sex)])
#       sex
# group   F  M
#   CTRL 17 20
#   LGA  21 21


#COVARIATES ANALYSIS====
# check no vars too much unique individual by levels (singleton)
bool_vars<-c("latino","preterm","GDM","drugs","etoh", "smoking")
categorical_vars<-c(bool_vars,"group","sex","group_sex","ethnicity","lab","batch","date","DNA.extraction","sequencing","year")

mtd[,(categorical_vars):=lapply(.SD,as.factor),.SDcols=categorical_vars]

lapply(mtd[,.SD,.SDcols=categorical_vars],function(x)sum(table(x)==1)) #exclude sequencing and date

mtd<-mtd[,-c("sequencing","date","drugs")]
library(ggrepel)
meth_mat<-as.matrix(data.frame(methf,row.names = "cpg_id")[,mtd$sample])
pca<-prcomp(t(meth_mat))
saveRDS(pca,"outputs/01-lga_vs_ctrl_limma_DMCs_analysis/pca_lgactrl.rds")

pc_mtd<-merge(mtd,data.table(pca$x,keep.rownames = "sample"))

fwrite(pc_mtd,"datasets/cd34/metadata_pcs_cl_190421.csv",sep=";")

ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=group))
ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=ponderal_index))
ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC3,col=GDM))
#PC2, unexplained variation so include in the model

ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC6,col=mat.age))

ggplot(pc_mtd)+geom_point(aes(x=PC2,y=PC9,col=latino))

mtd[is.na(mat.age)]
mtd[is.na(latino)]

ggplot(pc_mtd,aes(x=PC1,y=PC2,col=GDM))+geom_point()+ggrepel::geom_label_repel(aes(label=sample))
mtd[sample%in%c("CBP186","CBP253","CBP252","CBP205")]

ggplot(pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=library_complexity))


pc_mtd[is.na(ponderal_index)]
pval_mat<-CorrelCovarPCs(pca =pca ,mtd,rngPCs =1:10,res = "pval",seuilP = 1) #batch, seq depth,group_complexity,group , ethnicity
fwrite(data.table(pval_mat,keep.rownames = "covar"),fp(out,"pval_correl_covar_pcs.csv"),sep=";")


r2_mat<-CorrelCovarPCs(pca =pca ,mtd,rngPCs =1:10,res = "r2",seuilP = 1) #batch, seq depth,group_complexity,group , ethnicity
fwrite(data.table(r2_mat,keep.rownames = "covar"),fp(out,"r2_correl_covar_pcs.csv"),sep=";")

meth.influencing.vars<-rownames(pval_mat)[rowSums(pval_mat<0.01)>0]
meth.influencing.vars<-c(meth.influencing.vars,"mat.age","latino","group")

mtd2<-copy(mtd)
mtd2[,library_complexity:=group_complexity_fac]
CorrelCovarPCs(pca =pca ,mtd2[,c("sample","batch","mat.age","group","sex",
                                "latino","ponderal_index","weight.g","length.cm",
                                "head_circumf.cm.","gest.age.wk","DNA.extraction","library_complexity")],rngPCs =1:10,res = "pval",seuilP = 0.1) #batch, seq depth,group_complexity,group , ethnicity

# correl between covars and group
infl.vars.fac<-meth.influencing.vars[meth.influencing.vars%in%categorical_vars]
infl.vars.num<-meth.influencing.vars[meth.influencing.vars%in%numerical_vars]

cor_nums<-sapply(infl.vars.num,function(var1){
  r2s<-sapply(infl.vars.num,function(var2){
    f<-as.formula(paste(var1,"~",var2))
    return(summary(lm(f,mtd))$adj.r.squared)
    })
  return(r2s)
  })
cor_nums<-data.matrix(cor_nums)
pheatmap(cor_nums,display_numbers = T,cluster_rows = F,cluster_cols = F) 

pvals_num_fac<-sapply(infl.vars.num,function(var1){
  pvals<-sapply(infl.vars.fac,function(var2){
    f<-as.formula(paste(var1,"~",var2))
    return(anova(lm(f,mtd))$Pr[1])
    })
  return(pvals)
  })
pvals_num_fac<-data.matrix(pvals_num_fac)
plotPvalsHeatMap(pvals_num_fac) 
summary(lm(mtd$seq.depth~mtd$group))

pvals_fac<-sapply(infl.vars.fac,function(var1){
  pvals<-sapply(infl.vars.fac,function(var2)correl(mtd[,.SD,.SDcols=c(var1,var2)],verbose = T))
  return(pvals)
  })
pvals_fac<-data.matrix(pvals_fac)
pvals_fac[is.na(pvals_fac)]<-1
pheatmap(-log10(pvals_fac),display_numbers = T,cluster_rows = F,cluster_cols = F) #no correl
table(mtd[,.(group,GDM)])

#all this covars necessary to explain PC1?
summary(lm(PC1~group_complexity_fac+group+mat.age+latino,data = pc_mtd)) #seq.depth and mat age not sig in PC1
summary(lm(PC1~group_complexity_fac+group+mat.age+latino+seq.depth,data = pc_mtd)) #seq.depth to rm
summary(lm(PC1~group_complexity_fac+group+latino,data = pc_mtd)) #mat.age does not help to explain PC1 but explain well PC6 and PC26
ggplot(pc_mtd)+geom_point(aes(x=PC6,y=PC26,col=mat.age))
#and corralted with pct0
pc_mtd[,pct_zero:=colSums(meth_mat[,pc_mtd$sample]==0)/nrow(meth_mat)]
ggplot(pc_mtd)+geom_point(aes(x=mat.age,y=pct_zero,col=mat.age))


#include group_sex in the model instead of group ?
summary(lm(PC1~group_complexity_fac+group_sex+latino,data = pc_mtd)) #can see that effect ++ in lga F than Lga M, and CtrlM not sig



# so include batch,mat.age, group_complexity_fac,group_sex and PC2 in the model 

#DATA MODELING and DMC analysis with limma====
# exclude samples without all necessary clinical infos
vars_to_include<-c("batch","mat.age","group_complexity_fac","group_sex","latino","PC2")
mtd_f<-pc_mtd[,to_keep:=rowSums(is.na(.SD))==0,.SDcols=vars_to_include][to_keep==T]
fwrite(mtd_f,"datasets/cd34/metadata_cl_pcs_040521.csv")
mtd_f<-fread("datasets/cd34/metadata_cl_pcs_040521.csv")

nrow(mtd_f)
table(mtd_f$group,mtd_f$batch)
  #    1  2
  # CTRL 18 16
  # LGA  20 16

table(mtd_f$group,mtd_f$sex)
formule<- ~0 + group_sex   + batch+ group_complexity_fac +mat.age  + latino + PC2

mtd_f[,group_sex:=factor(group_sex,levels = unique(mtd_f$group_sex))]
design<-model.matrix(formule,data = data.frame(mtd_f,row.names = "sample"))
fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_f$sample], design)

cont.matrix <- makeContrasts(C.L = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexLGA_F+group_sexLGA_M)",
                             levels=design)

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

res<-data.table(topTable(fit2,coef = "C.L",n = Inf),keep.rownames = "cpg_id")

fwrite(res,fp(out,"res_limma.tsv.gz"),sep="\t")
res<-fread(fp(out,"res_limma.tsv.gz"),sep="\t")

ggplot(res)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=adj.P.Val<0.1&abs(logFC)>25),size=0.1)+
  scale_color_manual(values = c("grey","red"))+theme_minimal()
ggsave(fp(out,"volcano_plot.png"))


#logFC == %difference ?
methf[,avg.lga:=rowMeans(.SD),.SDcols=mtd_f[group=="LGA"]$sample]
methf[,avg.ctrl:=rowMeans(.SD),.SDcols=mtd_f[group=="CTRL"]$sample]

res_m<-merge(res,methf[,.(cpg_id,avg.lga,avg.ctrl)],by="cpg_id")
res_m[,meth_diff:=avg.lga-avg.ctrl]
res_m[,meth_fc:=avg.lga/avg.ctrl]
res_m[,meth_lfc:=log2(meth_fc)]

plot(res_m$logFC,res_m$meth_diff)
plot(res_m[adj.P.Val<0.1]$logFC,res_m[adj.P.Val<0.1]$meth_diff)

plot(res_m$logFC,res_m$meth_fc)
plot(res_m$logFC,res_m$meth_lfc)


cor(res_m$logFC,res_m$meth_diff)^2 #r2=84%


#DMR
#run 01B


#Validation cohorts
#n by cohorts
table(mtd_f$batch)
vars_to_include<-c("mat.age","group_complexity_fac","group_sex","latino","PC2")


res_batch<-Reduce(rbind,lapply(1:2,function(b){
  mtd_fc<-mtd_f[batch==b]
  formule<- ~0 + group_sex  + group_complexity_fac +mat.age  + latino + PC2
  design<-model.matrix(formule,data = data.frame(mtd_fc,row.names = "sample"))
  fit <- lmFit(data.frame(methf,row.names = "cpg_id")[,mtd_fc$sample], design)
  cont.matrix <- makeContrasts(C.L = "(group_sexCTRL_F+group_sexCTRL_M)-(group_sexLGA_F+group_sexLGA_M)",
                             CF.LF="group_sexCTRL_F-group_sexLGA_F",
                             CM.LM="group_sexCTRL_M-group_sexLGA_M",
                             levels=design)
  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  res<-Reduce(rbind,lapply(colnames(cont.matrix), function(comp)data.table(topTable(fit2,coef = comp,n = Inf),keep.rownames = "cpg_id")[,compa:=comp]))
  return(res[,batch:=b])
}))

fwrite(res_batch,fp(out,"res_limma_cohorts.tsv.gz"),sep="\t")
res_batch<-fread(fp(out,"res_limma_cohorts.tsv.gz"),sep="\t")
res_batch<-res_batch[compa=="C.L"]
table(res_batch[adj.P.Val<=0.1,.(compa,batch)])
ggplot(res_batch)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>25),size=0.1)+
  facet_wrap("batch")+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()
ggsave(fp(out,"fig1A.volcanos_limma_cohorts_C.L.png"),width = 12,height = 5)

#interesection
# renv::install("VennDiagram")
library(VennDiagram)
# Fonction d'aide pour afficher le diagramme de Venn
display_venn <- function(x, ...){
  require(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(list(batch1=res_batch[batch==1&P.Value<0.001&abs(logFC)>25]$cpg_id,
                  batch2=res_batch[batch==2&P.Value<0.001&abs(logFC)>25]$cpg_id
                  )) #only 12/2200 in common..

inter<-intersect(res_batch[batch==1&P.Value<0.001&abs(logFC)>25]$cpg_id,
                 res_batch[batch==2&P.Value<0.001&abs(logFC)>25]$cpg_id)
summary(res_batch[cpg_id%in%inter]$adj.P.Val)
summary(res_batch[P.Value<0.001&abs(logFC)>25]$adj.P.Val) #et mm pas particulirement signif

over_repr_test_simple(set1 =res_batch[batch==1&P.Value<0.001&abs(logFC)>25]$cpg_id ,
                      set2 = res_batch[batch==2&P.Value<0.001&abs(logFC)>25]$cpg_id,
                      size_universe = length(unique(res_batch$cpg_id))) #0.04

#inter au niveau du gene ?
cpgs_annot<-fread("outputs/02A-CpGs_annotations/cpgs_closest_gene_tss_linked_within_200kb_around.tsv")
res_b<-merge(res_batch,cpgs_annot)

display_venn(list(batch1=unique(res_b[batch==1&P.Value<0.001&abs(logFC)>25]$gene),
                  batch2=unique(res_b[batch==2&P.Value<0.001&abs(logFC)>25]$gene)
                  )) #197/1222

over_repr_test_simple(set1 =res_b[batch==1&P.Value<0.001&abs(logFC)>25]$gene ,
                      set2 = res_b[batch==2&P.Value<0.001&abs(logFC)>25]$gene,
                      size_universe = length(unique(res_b$gene))) #5.933862e-26

#+++ inter au niveau du gene mais pas au niveau du cpg, due a ++ cpg lié à ces gene ?
res_perm<-res_b[batch==1,.(cpg_id,gene)]
ps<-sapply(1:100, function(i){
  set.seed(i)
  return(over_repr_test_simple(res_perm[cpg_id%in%sample(cpg_id,2200)]$gene,
                      res_perm[cpg_id%in%sample(cpg_id,2200)]$gene,
                      size_universe = length(unique(res_perm$gene))))
  }
  )
  
summary(ps)
set.seed(123)
res_perm[cpg_id%in%sample(cpg_id,2200)]$gene
res_perm[cpg_id%in%sample(cpg_id,2200)]$gene

display_venn(list(batch1=unique(res_perm[cpg_id%in%sample(cpg_id,2200)]$gene),
                  batch2=unique(res_perm[cpg_id%in%sample(cpg_id,2200)]$gene)
                  )) #374/1590

#=> intersection pas signif



