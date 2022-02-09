#Test GeneScore on other datasets
#I) Le GeneScore toujours enrichit dans les DEG ? 


#2020-12-10 Mono old vs young
library(data.table)
library(ggplot2)
library(stringr)
library(limma)
source("scripts/utils/methyl_utils.R")
set.seed(1234)
options(stringsAsFactors = F)

#explor, easy to analyze ?
mo_meth<-fread("datasets/monocyte_aging/filtered_cytosines_freq.tsv.gz")
plot(density(na.omit(as.matrix(mo_meth[,-1]))))

sum(is.na(as.matrix(mo_meth[,-1])))/length(as.matrix(mo_meth[,-1]))

mtd<-data.table(sample=colnames(mo_meth)[-1],age=ifelse(str_detect(colnames(mo_meth)[-1],"O"),"Old","Young"))
table(mtd$age)
#=> oui

#filtration des NA
mo_meth<-data.frame(mo_meth,row.names = 1)
mo_meth_F<-mo_meth[!rowSums(is.na(mo_meth))==0,]
dim(mo_meth_F) #325766     40


#limma model
mtd[,age:=as.factor(age)]
design<-model.matrix(~0+age,data = mtd)
fit <- lmFit(mo_meth_F, design)

cont.matrix <- makeContrasts(young.old= "ageYoung-ageOld",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
res<-topTable(fit2,coef = "young.old",n = Inf)
res<-data.table(res,keep.rownames = "chrBase")
res<-res[,pval:=P.Value][,pval_adj:=adj.P.Val][,meth.change:=logFC][,avg_pct.meth:=AveExpr][,.(chrBase,meth.change,pval,avg_pct.meth)]
res

#volcano
ggplot(res[!is.na(pval)])+geom_point(aes(x=meth.change,y=-log10(pval),col=abs(meth.change)>20&pval<0.001))

#GeneScore
#need to make the CpG reg ref from chr pos.
#For LinksWeight, need :tss_dist, eQTR
res[,chr:=str_extract(chrBase,"chr[0-9XYM]+")]
res[,pos:=sapply(chrBase,function(x)as.numeric(strsplit(x,"\\.")[[1]][2]))]
res[,chr:=str_extract(chr,"[0-9XYM]+")]
res[,locisID:=paste(chr,pos,sep=":")]
fwrite(res[,.(chr,pos,meth.change,pval,avg_pct.meth,locisID)],"analyses/valid_gene_score/monocytes/res.csv")
cpgs<-res[,.(chr,pos)]

#need gene closest info
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
cpgs[,start:=pos][,end:=pos]

cpgs<-makeGRangesFromDataFrame(cpgs)
anno_cpgs <- annotatePeak(cpgsGR, tssRegion=c(-3000, 3000),
                       TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                       annoDb="org.Hs.eg.db")

anno_cpgs<-as.data.frame(anno_cpgs)
head(anno_cpgs,10)
anno_cpgs<-data.table(anno_cpgs)[,chr:=seqnames][,pos:=start][,tss_dist:=distanceToTSS][,gene:=SYMBOL][,.(chr,pos,tss_dist,gene)]



#For RegWeght need : feature_type_name, (chromatinFeature)
#feature_type_name

ensemblReg<-fread("ref/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz")
ensemblReg<-ensemblReg[,chr:=V1][,start:=V4][,end:=V5][,feature_type_name:=V3][,.(chr,start,end,feature_type_name)]
ensemblReg<-ensemblReg[!str_detect(chr,"\\.")]
fwrite(ensemblReg,"ref/ensembl_regulatory_hg38.bed",sep="\t")
#need 1) convert pos in hg19, 2) use bedtools interesect
#1)convert avec crossmap
#need bed3 format (“chrom”, “start”, “end”) 

fwrite(head(ensemblReg[,.(chr,start,end)]),"ref/test_hg38to19_conv.bed3",sep = "\t") 
#works with : CrossMap.py bed hg38ToHg19.over.chain.gz ensembl_regulatory_hg38.bed
#but not need header
#now for all pos
#with CrossMap.py bed hg38ToHg19.over.chain.gz temp_hg38to19_conv.bed3 converted_ensemblReg_hg38tohg19.bed3
ensemblReg19<-fread("ref/ensembl_regulatory_hg19.bed",col.names = c("chr","start","end","feature_type_name"))
#need put res in same chr format
anno_cpgs[,chr:=str_extract(chr,"[0-9XYM]+")]

#then run cmd like : python CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed3 test.hg19.bed3
#then intersect with bedtools intersect

dir.create("analyses/valid_gene_score/monocytes")
bed_cpgs<-anno_cpgs[,.(chr,pos)][,start:=pos][,end:=pos+1]
fwrite(bed_cpgs[order(chr,start)][,.(chr,start,end)],"analyses/valid_gene_score/monocytes/temp_cpgs.bed",sep = "\t",col.names = F)
#sort bed ensembl
bed_ensembl<-fread("ref/ensembl_regulatory_hg19.bed")
fwrite(bed_ensembl[order(V1,V2)],"ref/ensembl_regulatory_hg19.bed",col.names = F,sep = "\t")
#reorder 
#run in bash : bedtools intersect -wa -wb -a analyses/valid_gene_score/monocytes/temp_cpgs.bed -b ref/ensembl_regulatory_hg19.bed > analyses/valid_gene_score/monocytes/cpgs_in_ensembl_reg.bed
cpgs_ensembl<-fread("analyses/valid_gene_score/monocytes/cpgs_in_ensembl_reg.bed")
cpgs_ensembl<-fread("analyses/valid_gene_score/monocytes/cpgs_in_ensembl_reg.bed",
                   select = c(1,2,7),
                   col.names = c("chr","pos","feature_type_name"))
cpgs_ensembl[,locisID:=paste(chr,pos,sep=":")]
cpgs_ensembl[,feature_type_name:=paste(feature_type_name,collapse="/"),by="locisID"]
cpgs_ensembl<-unique(cpgs_ensembl)
anno_cpgs<-merge(anno_cpgs,cpgs_ensembl,by=c("chr","pos"),all=T)
fwrite(anno_cpgs,"analyses/valid_gene_score/monocytes/anno_cpgs.csv",sep = ";")



#eQTR 

meta_eqtr<-fread("ref/2020-06-23_meta.eQTR.csv")
#make in bed format for intersect
meta_eqtr[,chr:=str_extract(chr,"[0-9XYM]+")]
fwrite(eqtr[,.(chr,start.eQTR,end.eQTR,gene,RegScore)],"ref/meta.eQTR.bed",sep="\t")

cpgs_in_meta_eqtr<-fread("analyses/valid_gene_score/monocytes/cpgs_in_meta_eqtr.bed")
cpgs_in_meta_eqtr<-fread("analyses/valid_gene_score/monocytes/cpgs_in_meta_eqtr.bed",
                         select = c(1,2,7,8),
                         col.names = c("chr","pos","gene","avg.mlog10pval"))
cpgs_in_meta_eqtr[,in_eQTR:=TRUE]

anno_cpgs<-fread("analyses/valid_gene_score/monocytes/anno_cpgs.csv")
anno_cpgs

anno_cpgs<-merge(merge(anno_cpgs[,-"feature_type_name"],cpgs_in_meta_eqtr,by=c("chr","pos","gene"),all=T),unique(anno_cpgs[,.(chr,pos,feature_type_name)]),by=c("chr","pos"),all.x=T)
anno_cpgs[in_eQTR==T]
fwrite(anno_cpgs,"analyses/valid_gene_score/monocytes/anno_cpgs.csv",sep = ";")


blood_eqtr<-fread("ref/2020-05-31_Whole_Blood.eQTR_light.csv")
#make in bed format for intersect
blood_eqtr[,chr:=str_extract(chr,"[0-9XYM]+")]
fwrite(blood_eqtr[,.(chr,start.eQTR,end.eQTR,gene,RegScore)],"ref/blood.eQTR.bed",sep="\t")

cpgs_in_blood_eqtr<-fread("analyses/valid_gene_score/monocytes/cpgs_in_blood_eqtr.bed",
                         select = c(1,2,7,8),
                         col.names = c("chr","pos","gene","avg.mlog10pval"))
cpgs_in_blood_eqtr[,in_eQTR:=TRUE]
cpgs_in_blood_eqtr[,in_wb_eQTR:=TRUE]
cpgs_in_blood_eqtr[,chr:=as.character(chr)]
anno_cpgs<-fread("analyses/valid_gene_score/monocytes/anno_cpgs.csv")
anno_cpgs[in_eQTR==TRUE,in_meta_eQTR:=TRUE]

anno_cpgs<-merge(merge(anno_cpgs[,-"feature_type_name"],cpgs_in_blood_eqtr,all=T),unique(anno_cpgs[,.(chr,pos,feature_type_name)]),by=c("chr","pos"),all.x=T)
anno_cpgs[in_eQTR==T]
fwrite(anno_cpgs,"analyses/valid_gene_score/monocytes/anno_cpgs.csv",sep = ";")



library(data.table)
source("scripts/utils/methyl_utils.R")
anno_cpgs<-fread("analyses/valid_gene_score/monocytes/anno_cpgs.csv")
anno_cpgs[,locisID:=paste(chr,pos,sep=":")]
anno_cpgs[is.na(in_eQTR),in_eQTR:=FALSE]
anno_cpgs[in_eQTR==FALSE][duplicated(locisID)]
anno_cpgs

#ensregscore
anno_cpgs[,EnsRegScore:=sapply(feature_type_name,function(x){
  vecX<-strsplit(as.character(x),"/")[[1]]
  if(any(c("CTCF_binding_site","promoter","enhancer")%in%vecX)){
    score<-0.5
  }else if(any(c("open_chromatin_region","promoter_flanking_region")%in%vecX)){
    score<-0.25
  }else{
    score<-0
  }
  if("TF_binding_site"%in%vecX){
    score<-score+0.5
  }
  return(score)
})]

# Regulatory weight:
anno_cpgs[,RegWeight:=0.5+1.5*EnsRegScore]

# Links Weight[0.1-1] = 0.1+0.9*max(LinkScore [0-1])
#   - for CpG associated to a gene based on TSS proximity  

anno_cpgs[in_eQTR==FALSE,
          LinkScore:=ifelse(abs(tss_dist)<1000,
                            1,
                            ifelse(abs(tss_dist)<20000,
                                   0.5+0.5*sqrt(1000/abs(tss_dist)),
                                   0.5*sqrt(20000/abs(tss_dist)))),
          by="locisID"
 ]

anno_cpgs[in_eQTR==TRUE,LinkScore:=1] 

anno_cpgs[,LinksWeight:=max(LinkScore),by=c("locisID","gene")]

res<-fread("analyses/valid_gene_score/monocytes/res.csv")
res<-merge(res[,-c("chr","pos")],anno_cpgs[,-c("chr","pos")],all=T,by="locisID")
fwrite(res,"analyses/valid_gene_score/monocytes/res_annotated.csv")

res<-unique(res[order(locisID,gene,-LinkScore)],by=c("locisID","gene"))
res[,CpGScore:=(-log10(pval)/4*meth.change)*RegWeight*LinksWeight]



res[,n_cpg_gene:=.N,by=.(gene)]
res[,n_gene_cpg:=.N,by=.(locisID)]
res[n_gene_cpg>1]
res[,n_cpg_sig_gene:=sum(pval<0.001),by=.(gene)]
res[,nCpGWeight:=(1/sum(1/(abs(CpGScore)+1)))^(1/4),by="gene"]
res[,GeneScore:=sum(CpGScore)*nCpGWeight,by="gene"]
fwrite(res,"analyses/valid_gene_score/monocytes/res_annotated_with_genescore.csv")
library(ggplot2)
library(patchwork)
p1<-ggplot(unique(res,by="gene"))+
  geom_boxplot(aes(x = as.factor(n_cpg_gene),y =nCpGWeight )) 

p2<-ggplot(unique(res,by="gene"))+
  geom_boxplot(aes(x = as.factor(n_cpg_gene),y =GeneScore )) 
p1+p2
resg<-unique(res[order(pval)],by="gene")
resg[abs(GeneScore)>20]

#DEG
library(data.table)
library(DESeq2)
rna<-fread("datasets/monocyte_aging/TableS6_monocyte_rnaseq_counts2.csv")
rna<-data.frame(rna[,-"gene_id"][!is.na(gene_symbol)][order(gene_symbol,-YD18)][!duplicated(gene_symbol)],row.names = 1)
head(rna)
dds1 <- DESeqDataSetFromMatrix(rna, 
                               colData = data.frame(mtd[,cell_type:="Mo"],row.names =1 )[colnames(rna),], 
                               design = ~ age)
dds1 <- DESeq(dds1)
resultsNames(dds1)
# resultsNames(dds)
res1 <- results(dds1, 
                contrast = c("age","Old","Young"),
                alpha = 0.05)

res1 <- lfcShrink(dds1, 
                  contrast = c("age","Old","Young"),
                  res=res1)
res1<-data.table(data.frame(res1),keep.rownames = "gene")

resg<-merge(res1,resg)

fwrite(resg[order(pvalue)],
       "analyses/valid_gene_score/monocytes/res_with_deg.csv",
       sep=";")
resg[is.na(padj),padj:=1]
ggplot(resg)+geom_boxplot(aes(x=padj<0.05,y=abs(GeneScore)))

ggplot(resg[padj<0.2])+geom_point(aes(x=abs(log2FoldChange),y=abs(GeneScore)))


res_meth1[is.na(padj),padj:=1]
ggplot(res_meth1)+geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.05))+ scale_color_manual(values = c("grey","red")) + theme_minimal() +theme(legend.position = "none")
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_volcano_plot_DEG.png")

ggplot(res_meth1[padj<1])+geom_point(aes(x=log2FoldChange,y=GeneScore,col=padj<0.05))+
  scale_x_continuous(limits = c(-2.5,2.5))+scale_color_manual(values = c("grey","red")) +
  theme_minimal() +theme(legend.position = "none")
ggsave("analyses/paper/ctrlM_vs_ctrlF/2020-10-06_plot_DEG_vs_DMG.png")


resp<-fread("datasets/monocyte_aging/TableS3_somascan_values_log2RFU_in_row.csv")
resp<-resp[,gene:=EntrezGeneSymbol][,.(gene,logFC,P.Value,adj.P.Val)]

resp[str_detect(gene," "),gene:=sapply(gene,function(x)strsplit(x," ")[[1]][1])]

resp[str_detect(gene," ")]
resp<-merge(resg,resp,by="gene")
ggplot(resp)+geom_boxplot(aes(x=adj.P.Val<0.1,y=log(abs(GeneScore))))
ggplot(resp)+geom_boxplot(aes(x=adj.P.Val<0.1,y=log(abs(meth.change))))
ggplot(resp)+geom_boxplot(aes(x=adj.P.Val<0.1,y=-log10(pval)))
ggplot(resp)+geom_point(aes(x=abs(logFC),y=abs(GeneScore)))


resp[adj.P.Val<0.1]


#2) TCGA
#2020-10-19 glioma young vs old

#package useful
library(data.table)
library(stringr)
#pour la norm
# renv::install("bioc::ENmix")
#library(ENmix)

library(limma)


mtd<-merge(fread("datasets/TCGA/old_male_glioma/gdc_sample_sheet.2020-09-10.tsv")
           [`Data Type`=="Methylation Beta Value"&`Sample Type`=="Primary Tumor"]
           [,sample:=`Sample ID`]
           [,case:=`Case ID`]
           [,.(`File ID`,`File Name`,case,sample)]
           [,age:="old"],
           fread("datasets/TCGA/young_male_glioma/gdc_sample_sheet.2020-09-10.tsv")
           [`Data Type`=="Methylation Beta Value"&`Sample Type`=="Primary Tumor"]
           [,sample:=`Sample ID`]
           [,case:=`Case ID`]
           [,.(`File ID`,`File Name`,case,sample)]
           [,age:="young"],all=T)

mtd<-unique(mtd,by="case")
sample<-mtd$case[1]
meth<-fread(file.path("datasets/TCGA/young_male_glioma",
                      mtd$`File ID`,
                      mtd$`File Name`)[1])[,locisID:=`Composite Element REF`][,(sample):=Beta_value][,.SD,.SDcols=c("locisID",sample)]


for(i in 2:nrow(mtd)){
  
  sample<-mtd$case[i]
  print(sample)
  if(!sample%in%colnames(meth)){
    path<-file.path("datasets/TCGA",paste0(mtd$age[i],"_male_glioma"),mtd$`File ID`[i],mtd$`File Name`[i])
    
    meth<-merge(meth,
                fread(path)[,locisID:=`Composite Element REF`][,(sample):=Beta_value][,.SD,.SDcols=c("locisID",sample)],
                by="locisID",
                all=T)
    
  }
  
  
}
dir.create("analyses/valid_gene_score/tcga_glioma/")
fwrite(meth,"datasets/TCGA/methylation_b_value_old_vs_young_male.csv")
#QC 
meth<-meth[rowSums(is.na(as.matrix(meth[,-"locisID"])))==0,]

#rm SNP etc
anno<-fread("ref/HM450.hg38.manifest.tsv.gz")
anno[MASK_general==F]
anno[,locisID:=probeID]
meth<-merge(meth,anno[,.(locisID,MASK_general,designType)])
meth<-meth[MASK_general==F]
fwrite(meth,"datasets/TCGA/methylation_b_value_old_vs_young_glioma_male_filtered.csv",sep = ";")


plot(density(as.matrix(meth[designType=="I"][,.SD,.SDcols=-c("locisID","MASK_general","designType")])))
lines(density(as.matrix(meth[designType=="II"][,.SD,.SDcols=-c("locisID","MASK_general","designType")])),col=2)
plot(density(as.matrix(data_met[,.SD,.SDcols=mtd$sample])))

#enmix norm type I and II probe [refaire fct pour accepter tcga data]
library(ENmix)
?getB
?rcp
rcp_from_beta<-function (beta, dist = 25, quantile.grid = seq(0.001, 0.999, 
                                               by = 0.001), qcscore = NULL, nbthre = 3, detPthre = 1e-06) 
{
  if (!is(mdat, "methDataSet") & !is(mdat, "MethylSet")) {
    stop("the input needs to be of class 'methDataSet' or 'MethylSet'")
  }
  raw.M <- B2M(beta)
  if (is(mdat, "methDataSet")) {
    annotation = rowData(mdat)
    names(annotation)[which(names(annotation) == "Infinium_Design_Type")] = "Type"
    annotation = annotation[annotation$Type %in% c("I", 
                                                   "II"), ]
    rownames(annotation) = annotation$Name
  }
  else if (is(mdat, "MethylSet")) {
    annotation <- getAnnotation(mdat)
  }
  annotation = annotation[as.vector(annotation$Name) %in% 
                            rownames(beta), ]
  probe.II.Name = annotation$Name[annotation$Type == "II"]
  annotation = annotation[order(annotation$chr, annotation$pos), 
                          ]
  anno1 = annotation[1:(nrow(annotation) - 1), ]
  anno2 = annotation[2:nrow(annotation), ]
  flag = (abs(anno1$pos - anno2$pos) < dist & anno1$chr == 
            anno2$chr & anno1$Relation_to_Island == anno2$Relation_to_Island & 
            anno1$Type != anno2$Type)
  anno1 = anno1[flag, ]
  anno2 = anno2[flag, ]
  probe.I = anno1$Name
  probe.II = anno2$Name
  probe.I[anno2$Type == "I"] = anno2$Name[anno2$Type == "I"]
  probe.II[anno1$Type == "II"] = anno1$Name[anno1$Type == 
                                              "II"]
  raw.M.t = raw.M[c(probe.I, probe.II), ]
  if (is.null(qcscore)) {
  }
  else if ((sum(!(rownames(raw.M.t) %in% rownames(qcscore$detP))) + 
            sum(!(colnames(raw.M.t) %in% colnames(qcscore$detP)))) > 
           0) {
    stop("Wrong qcscore matrix, please check...\n")
  }
  else {
    temp <- qcscore$nbead < nbthre | qcscore$detP > detPthre
    temp = temp[rownames(raw.M.t), ]
    temp = temp[, colnames(raw.M.t)]
    raw.M.t[temp] = NA
  }
  M.II <- raw.M.t[probe.II, ]
  M.I <- raw.M.t[probe.I, ]
  qtl <- function(x) quantile(x, quantile.grid, na.rm = TRUE)
  M.I = apply(M.I, 2, qtl)
  M.II = apply(M.II, 2, qtl)
  beta.est <- mat.or.vec(2, ncol(beta))
  for (i in 1:ncol(beta)) {
    index <- (M.II[, i] != Inf & M.II[, i] != -Inf & M.I[, 
                                                         i] != Inf & M.I[, i] != -Inf)
    X <- cbind(rep(1, sum(index)), M.II[index, i])
    Y <- M.I[index, i]
    beta.est[, i] <- solve(t(X) %*% X) %*% t(X) %*% Y
  }
  M.II.all <- raw.M[probe.II.Name, ]
  M.II.new <- mat.or.vec(nrow(M.II.all), ncol(M.II.all))
  for (i in 1:ncol(M.II.all)) {
    M.II.new[, i] <- beta.est[1, i] + beta.est[2, i] * M.II.all[, 
                                                                i]
  }
  M.II.new[M.II.all == Inf] <- Inf
  M.II.new[M.II.all == -Inf] <- (-Inf)
  beta[probe.II.Name, ] <- M2B(M.II.new)
  beta
}
beta_mat<-as.matrix(data.frame(meth[,-c("designType","MASK_general")],row.names = "locisID"))
head(beta_mat)
beta_norm<-rcp_from_beta(beta) #doestn worh

# [end refaire fct norm pour accepter tcga data]


#for the moment, go through

#LIMMA
mtd[,age:=as.factor(age)]
formule<- ~0 + age 

design<-model.matrix(formule,data = mtd)

fit <- lmFit(beta_mat, design)
colnames(design)
cont.matrix <- makeContrasts(old.young = "ageold-ageyoung",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
results <- decideTests(fit2)
sum(abs(results)) #6809

res<-topTable(fit2,coef = "old.young",n = Inf)
res<-data.table(res,keep.rownames = "locisID")

res<-res[,pval:=P.Value][,pval_adj:=adj.P.Val][,meth.change:=logFC][,avg_pct.meth:=AveExpr][,.(locisID,meth.change,pval,pval_adj,avg_pct.meth)]
fwrite(res,"analyses/valid_gene_score/tcga_glioma/old_vs_young_glioma_meth.csv",sep=";")

res<-fread("analyses/valid_gene_score/tcga_glioma/old_vs_young_glioma_meth.csv")
#volcano
ggplot(res[!is.na(pval)])+geom_point(aes(x=meth.change,y=-log10(pval),col=abs(meth.change)>0.2&pval<0.001))

res[abs(meth.change)>0.2&pval<0.001]

res[meth.change>0.2&pval<0.001] #858
res[meth.change>0.4&pval<0.001]
summary(unlist(meth[locisID=="cg23746497",.SD,.SDcols=mtd[age=="old"]$case]))#hypermet in old
summary(unlist(meth[locisID=="cg23746497",.SD,.SDcols=mtd[age=="young"]$case]))

anno[probeID=="cg23746497",.(gene,)]


#anno cpgs
anno_genes<-readRDS("ref/HM450.hg38.manifest.gencode.v22.rds")
anno_genes<-data.table(as.data.frame(anno_genes),keep.rownames = "locisID")


anno_genes[,n.gene:=sapply(geneNames,function(x)length(strsplit(x,";")[[1]]))]
summary(anno_genes$n.gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   1.000   3.000   4.769   6.000 123.000 


anno_genes[,gene:=sapply(geneNames,function(x)strsplit(x,";")[[1]][1])]
anno_genes[,tss_dist:=sapply(distToTSS,function(x)as.numeric(strsplit(x,";")[[1]][1]))]
anno_genes[,gene_type:=sapply(transcriptTypes,function(x)strsplit(x,";")[[1]][1])]

for(i in 2:max(anno_genes$n.gene)){
  anno_genes_temp<-anno_genes[n.gene==i]
  anno_genes_temp[,gene:=sapply(geneNames,function(x)strsplit(x,";")[[1]][i])]
  anno_genes_temp[,tss_dist:=sapply(distToTSS,function(x)as.numeric(strsplit(x,";")[[1]][i]))]
  anno_genes_temp[,gene_type:=sapply(transcriptTypes,function(x)strsplit(x,";")[[1]][i])]
  
  anno_genes<-rbind(anno_genes,anno_genes_temp)
  
}
gt<-unique(anno_genes$gene_type)
major_gt<-c("protein_coding","miRNA","lincRNA","processed_transcript")
anno_genes<-anno_genes[,gene_type:=factor(gene_type,levels = c(major_gt,setdiff(gt,major_gt)))][order(locisID,gene_type,tss_dist)]
anno_genes[n.gene==3]
fwrite(anno_genes,"ref/all_anno_450K_hg38.csv",sep=";")


res<-merge(res,anno[,.(locisID,gene,tss_dist,gene_type)],by="locisID")[order(pval,gene_type,tss_dist)]


fwrite(res,"analyses/valid_gene_score/tcga_glioma/old_vs_young_glioma_meth_annotated.csv")

cpgs<-res[,.(chr,pos)]


#need gene closest info
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
cpgs[,start:=pos][,end:=pos]

cpgs<-makeGRangesFromDataFrame(cpgs)
anno_cpgs <- annotatePeak(cpgsGR, tssRegion=c(-3000, 3000),
                          TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                          annoDb="org.Hs.eg.db")

anno_cpgs<-as.data.frame(anno_cpgs)
head(anno_cpgs,10)
anno_cpgs<-data.table(anno_cpgs)[,chr:=seqnames][,pos:=start][,tss_dist:=distanceToTSS][,gene:=SYMBOL][,.(chr,pos,tss_dist,gene)]



#For RegWeght need : feature_type_name, (chromatinFeature)
#feature_type_name

ensemblReg<-fread("ref/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz")
ensemblReg<-ensemblReg[,chr:=V1][,start:=V4][,end:=V5][,feature_type_name:=V3][,.(chr,start,end,feature_type_name)]
ensemblReg<-ensemblReg[!str_detect(chr,"\\.")]
fwrite(ensemblReg,"ref/ensembl_regulatory_hg38.bed",sep="\t")
#need 1) convert pos in hg19, 2) use bedtools interesect
#1)convert avec crossmap
#need bed3 format (“chrom”, “start”, “end”) 

fwrite(head(ensemblReg[,.(chr,start,end)]),"ref/test_hg38to19_conv.bed3",sep = "\t") 
#works with : CrossMap.py bed hg38ToHg19.over.chain.gz ensembl_regulatory_hg38.bed
#but not need header
#now for all pos
#with CrossMap.py bed hg38ToHg19.over.chain.gz temp_hg38to19_conv.bed3 converted_ensemblReg_hg38tohg19.bed3
ensemblReg19<-fread("ref/ensembl_regulatory_hg19.bed",col.names = c("chr","start","end","feature_type_name"))
#need put res in same chr format
anno_cpgs[,chr:=str_extract(chr,"[0-9XYM]+")]

#then run cmd like : python CrossMap.py bed hg18ToHg19.over.chain.gz test.hg18.bed3 test.hg19.bed3
#then intersect with bedtools intersect

dir.create("analyses/valid_gene_score/monocytes")
bed_cpgs<-anno_cpgs[,.(chr,pos)][,start:=pos][,end:=pos+1]
fwrite(bed_cpgs[order(chr,start)][,.(chr,start,end)],"analyses/valid_gene_score/monocytes/temp_cpgs.bed",sep = "\t",col.names = F)
#sort bed ensembl
bed_ensembl<-fread("ref/ensembl_regulatory_hg19.bed")
fwrite(bed_ensembl[order(V1,V2)],"ref/ensembl_regulatory_hg19.bed",col.names = F,sep = "\t")
#reorder 
#run in bash : bedtools intersect -a analyses/valid_gene_score/monocytes/temp_cpgs.bed -b ref/ensembl_regulatory_hg19.bed > analyses/valid_gene_score/monocytes/cpgs_in_ensembl_reg.bed -wa -wb
cpgs_ensembl<-fread("analyses/valid_gene_score/monocytes/cpgs_in_ensembl_reg.bed")
cpgs_ensembl<-fread("analyses/valid_gene_score/monocytes/cpgs_in_ensembl_reg.bed",
                    select = c(1,2,7),
                    col.names = c("chr","pos","feature_type_name"))
cpgs_ensembl[,locisID:=paste(chr,pos,sep=":")]
cpgs_ensembl[,feature_type_name:=paste(feature_type_name,collapse="/"),by="locisID"]
cpgs_ensembl<-unique(cpgs_ensembl)
anno_cpgs<-merge(anno_cpgs,cpgs_ensembl,by=c("chr","pos"),all=T)
fwrite(anno_cpgs,"analyses/valid_gene_score/monocytes/anno_cpgs.csv",sep = ";")
#eQTR [TODO]
#II) Mieux que DMR ?
