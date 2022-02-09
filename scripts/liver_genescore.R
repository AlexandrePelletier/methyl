#Liver GeneScore
#data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5383522/#sup1

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

#[end in BASH]

library(data.table)
library(stringr)
library(limma)
library(ggplot2)
dir.create("analyses/liver/")

features<-fread("ref/Chromatin_Annot/Liver/ChromHMM_output_GSMreads/liver_6_segments.bed",col.names = c("chr","start","end","state"))
features[,state:=as.numeric(str_remove(state,"E"))]
annot_feat<-data.table(state=1:6,feature=c("Heterochromatin","Heterochromatin","Gene Body","Inactive Enhancer","Active Enhancer","Promoter"))
features<-merge(features,annot_feat,all.x=T)
fwrite(features[,.(chr,start,end,state,feature)],"ref/Chromatin_Annot/Liver/chromatin_features_reg.bed",sep = "\t")

#Get DMC for Liver
meth<-fread("datasets/liver/GSE82176_Combined_processed_HELP_tagging_data.txt")
meth
samples<-colnames(meth)[5:ncol(meth)]
plot(density(as.matrix(meth[,.SD,.SDcols=samples])))

#limma model
mtd<-data.table(sample=samples,group=str_extract(samples,"C|NT|T"))
mtd[,group:=as.factor(group)]
design<-model.matrix(~0+group,data = mtd)
fit <- lmFit(data.frame(meth,row.names = 1)[,samples], design)

colnames(design)
cont.matrix <- makeContrasts(nt.ctrl= "groupNT-groupC",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
res<-topTable(fit2,coef = "nt.ctrl",n = Inf)
res<-data.table(res,keep.rownames = "cpgID")
res<-res[,pval:=P.Value][,pval_adj:=adj.P.Val][,meth.change:=logFC][,avg_pct.meth:=AveExpr][,.(chrBase,meth.change,pval,avg_pct.meth)]
res<-res[,cpgID:=as.numeric(cpgID)]
res[abs(meth.change)>0.20&pval<0.001]
res[abs(meth.change)>0.20&pval_adj<0.05]

#volcano
p<-ggplot(res[!is.na(pval)])+geom_point(aes(x=meth.change,y=-log10(pval),col=abs(meth.change)>0.20&pval_adj<0.05))

ggsave("analyses/liver/volcano_DMC_NTvsCtrl.png",p)
meth[,cpgID:=as.numeric(X.tid)]

res<-merge(res,meth[,.(cpgID,chr,pos)])
fwrite(res[,.SD,.SDcols=c("chr","pos",colnames(res)[-c(ncol(res)-1,ncol(res))])],"analyses/liver/res_NTvsCtrl_liver.csv",sep=";")

#ANNOT CPGs
res<-fread("analyses/liver/res_NTvsCtrl_liver.csv")

cpgs<-res[,.(chr,pos,cpgID)]
rm(res)
#1)cpgs_genes_links
#by distance to tss
tss<-fread("ref/hg19/hg19_refseq_curated.txt.gz")
tss<-tss[,chr:=chrom][strand=="+",start:=txStart][strand=="-",start:=txEnd][,end:=start+1][,gene:=name2][,.(chr,start,end,gene,strand)]
fwrite(tss[order(chr,start)],"ref/hg19/tss_genes.bed",sep="\t")

#4 first genes within +/-200kb
fwrite(cpgs[,start:=pos-200e3][start<0,start:=0][,end:=pos+200e3][,.(chr,start,end,cpgID)][order(chr,start)],"analyses/liver/temp_cpgs_reg_200kb.bed",sep="\t")
#[in BASH]
bedtools intersect -a analyses/liver/temp_cpgs_reg_200kb.bed -b ref/hg19/tss_genes.bed -wa -wb > analyses/liver/genes_in_200kb_cpgs_regions.tsv
#[end in BASH]
genes_in_200kb<-fread("analyses/liver/genes_in_200kb_cpgs_regions.tsv",select = c(5,6,8,9,4),
                      col.names = c("chr","tss_pos","gene","strand","cpgID"))


cpgs_genes_tss<-merge(cpgs,genes_in_200kb,all.x=T)[,.(chr,pos,cpgID,gene,tss_pos,strand)]
cpgs_genes_tss[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]
cpgs[,ngene.cpg:=.N,by="cpgID"]
summary(cpgs$ngene.cpg)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #   1.00   16.00   28.00   33.54   46.00  135.00 

cpgs_genes_4<-cpgs_genes_tss[order(cpgID,abs(tss_dist))][,top4:=gene%in%unique(gene)[1:4],by="cpgID"][top4==T][,-"top4"]
#for >1 cpg_gene link (because multiple TSS), keep most closest cpg_gene link
cpgs_genes_4<-unique(cpgs_genes_4[order(cpgID,gene,abs(tss_dist))],by=c("cpgID","gene")) 
cpgs_genes_4
rm(cpgs_genes_tss)

fwrite(cpgs_genes_4,"analyses/liver/cpgs_4closest_gene_tss_linked_within_200kb_around.tsv",sep="\t")
file.remove("analyses/liver/temp_cpgs_reg_200kb.bed")

#by presence in eQTL region

#[in BASH]
bedtools intersect -wa -wb -a analyses/liver/temp_cpgs.bed -b ref/eQTL/liver_eQTR_hgnc_gene_hg19.bed > analyses/liver/cpgs_in_liver_eQTR_region.bed
#[end in BASH]
cpgs_eQTR<-fread("analyses/liver/cpgs_in_liver_eQTR_region.bed",select = c(4,8),col.names = c("cpgID","gene"))
cpgs_eQTR[,in_eQTR:=T]

#merge the 2 links
cpgs_genes_4[,in_eQTR:=F]
cpgs_genes<-merge(cpgs_genes_4[,.(cpgID,gene,in_eQTR,tss_dist)],cpgs_eQTR,all=T)
cpgs_genes[in_eQTR==T]

#merge with coord
cpgs_genes<-merge(cpgs[,.(cpgID,chr,pos)],cpgs_genes,all.x = T,by="cpgID")

fwrite(cpgs_genes,"analyses/liver/cpgs_4closest_gene_tss_linked_within_200kb_around_and_eQTR_linked_genes.tsv",sep="\t")

#2)cpgs_regulatory region

#ensembl regulatory reg matching
fwrite(unique(cpgs,by="cpgID")[,chr:=str_remove(chr,'chr')][,start:=pos][,end:=pos+1][,.(chr,start,end,cpgID)][order(chr,start)],"analyses/liver/temp_cpgs.bed",sep="\t",col.names = F)
#[in BASH]
bedtools intersect -wa -wb -a analyses/liver/temp_cpgs.bed -b ref/ensembl_regulatory/ensembl_regulatory_hg19.bed  > analyses/liver/cpgs_in_ensembl_regulatory_regions.bed
#[end in BASH]
cpgs_ensembl<-fread("analyses/liver/cpgs_in_ensembl_regulatory_regions.bed",
                    select = c(4,8),col.names = c("cpgID","ensembl_regulatory_domain"))

cpgs_ensembl
cpgs_ensembl[,ensembl_regulatory_domain:=paste(ensembl_regulatory_domain,collapse = "/"),by="cpgID"] #take a long time
cpgs_ensembl<-unique(cpgs_ensembl)
cpgs_ensembl

#chromatinfeature
#[in BASH]
bedtools intersect -wa -wb -a analyses/liver/temp_cpgs.bed -b ref/Chromatin_Annot/Liver/chromatin_features_reg.bed > analyses/liver/cpgs_in_liver_chromatin_features_region.bed
#[end in BASH]
cpgs_chromatin<-fread("analyses/liver/cpgs_in_liver_chromatin_features_region.bed",select = c(4,9),col.names = c("cpgID","chromatin_feature"))
cpgs_reg<-merge(cpgs_chromatin,cpgs_ensembl,all.x=T,by="cpgID")

fwrite(cpgs_reg,"analyses/liver/cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.tsv",sep="\t")

cpgs_reg<-fread("analyses/liver/cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.tsv",sep="\t")

cpgs_genes<-fread("analyses/liver/cpgs_4closest_gene_tss_linked_within_200kb_around_and_eQTR_linked_genes.tsv",sep="\t")

cpgs_anno<-merge(cpgs_reg,cpgs_genes,all=T,by="cpgID")

fwrite(unique(cpgs_anno[order(cpgID,gene)][,.(cpgID,chr,pos,in_eQTR,tss_dist,gene,chromatin_feature,ensembl_regulatory_domain)]),"analyses/liver/cpgs_annot.tsv",sep="\t")


#2020-12-02 Calculate GeneScore

library(data.table)
library(ggplot2)


#1) calculate cpg weight
cpgs<-fread("analyses/liver/cpgs_annot.tsv")
cpgs_genes<-cpgs[!is.na(gene)&gene!=""]
# a) linksWeight
cpgs_genes[in_eQTR==F,links_score:=sapply(abs(tss_dist),function(x){
  if(x<1000)return(1)
  else if(x<20000)return(0.5+0.5*sqrt(1000/x))
  else return(0.5*sqrt(20000/x))
    })]

cpgs_genes[in_eQTR==T,links_score:=1]

cpgs_genes[,links_weight:=max(links_score),by=c("cpgID","gene")]

cpgs_genes[links_score!=links_weight] #247074 match tss_dist-eQTR / 1M818k cpgs
cpgs_genes[in_eQTR==T] #1M707k CpGs-Gene linked by eQTR

#b) Regulatory Weights
cpgs<-unique(cpgs_genes,by=c("cpgID"))
unique(cpgs$chromatin_feature)
cpgs[,chromatin_score:=sapply(chromatin_feature,function(x){
      if(x%in%c("Promoter","Active Enhancer"))return(1)
      else if(x=="Inactive Enhancer")return(0.75)
      else if(x%in%"Gene Body")return(0.5)
      else return(0)
    })]



cpgs[,ensembl_reg_score:=sapply(ensembl_regulatory_domain,function(x){
    vecX<-strsplit(as.character(x),"/")[[1]]
    if(any(c("CTCF Binding Site","Promoter","Enhancer")%in%vecX)){
      score<-0.5
    }else if(any(c("Open chromatin","Promoter Flanking Region")%in%vecX)){
      score<-0.25
    }else{
      score<-0
    }
    if("TF binding site"%in%vecX){
      score<-score+0.5
    }
    return(score)
  })]

cpgs[,regul_weight:=(0.5+1.5*((chromatin_score+ensembl_reg_score)/2))]

cpgs_genes<-merge(cpgs_genes,cpgs[,.(cpgID,regul_weight)],by="cpgID")
fwrite(cpgs_genes,"analyses/liver/cpgs_genes_annot_and_weight_reference.tsv",sep="\t")

#2) merge with res and calculate CpG Score
cpgs_genes<-fread("analyses/liver/cpgs_genes_annot_and_weight_reference.tsv")
res<-fread("analyses/liver/res_NTvsCtrl_liver.csv")
res<-merge(res[,-c("chr","pos")],unique(cpgs_genes[order(cpgID,gene,-links_score)],by=c("cpgID","gene")),by=c("cpgID"))
res[,cpg_score:=(-log10(pval)/4*meth.change*100)*regul_weight*links_weight]

#3) GeneScore

res[,n_cpg.gene:=.N,by=.(gene)]
res[,n_cpg_sig.gene:=sum(pval_adj<0.05),by=.(gene)]
res[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
res[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

plot(density(unique(res,by="gene")$n_cpg.gene))
summary(unique(res,by="gene")$n_cpg.gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     1.0   139.0   213.0   256.7   319.0  2072.0 


ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()

ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=n_cpg_weight))+scale_x_log10()
ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg_sig.gene,y=gene_score))+scale_x_log10()
ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg_sig_nom3.gene,y=gene_score))+scale_x_log10()

res[order(-abs(gene_score),-abs(cpg_score))][pval_adj<0.05][gene_score>500]

plot(density(log10(abs(unique(res,by=c("gene"))$gene_score))))


unique(res[gene_score>100]$gene)

fwrite(res,"analyses/liver/res_gene_score_NTvsCtrl_liver.csv",sep=";")

#VALIDATION
#1) Calculate DE Analysis
library(data.table)
library(DESeq2)
library(stringr)
counts<-as.matrix(data.frame(fread("datasets/liver/GSE82177_Merged_Gene_Counts_HTseq.txt"),row.names = 1))
head(counts)

samples<-colnames(counts)
mtd<-data.table(sample=samples,group=str_extract(samples,"C|NT|T"))
mtd[,group:=as.factor(group)]
dds1 <- DESeqDataSetFromMatrix(counts, 
                               colData = data.frame(mtd[,tissue:="liver"],row.names =1 )[colnames(counts),], 
                               design = ~ group)
dds1 <- DESeq(dds1)
resultsNames(dds1)
# resultsNames(dds)
res1 <- results(dds1, 
                contrast = c("group","NT","C"),
                alpha = 0.05)

res1 <- lfcShrink(dds1, 
                  coef = "group_NT_vs_C" ,
                  res=res1)
res1<-data.table(data.frame(res1),keep.rownames = "gene_id")

library(biomartr)
getMarts()
getDatasets("ENSEMBL_MART_ENSEMBL")
dt<-data.table(getAttributes(mart ="ENSEMBL_MART_ENSEMBL" ,dataset = "hsapiens_gene_ensembl"))
dt[str_detect(name,"refseq")]
genes_translator <- biomart( genes      = unique(res1$gene_id), # genes that we wanted info
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = "hgnc_symbol", # attributes were selected with biomartr::getAttributes()
                                filters    = "refseq_mrna")

genes_translator<-data.table(genes_translator)
fwrite(genes_translator,"ref/refseq_mrna_to_hgnc_symbol_translator.csv",sep=";")
genes_translator<-genes_translator[,gene_id:=refseq_mrna][,gene:=hgnc_symbol][,-c("refseq_mrna","hgnc_symbol")]
res1<-merge(res1,genes_translator)
res1f<-res1[baseMean>0&!is.na(gene)]
res1f[padj<0.05]

res2<-merge(res,res1f,by="gene")



fwrite(res2[order(pvalue)],
       "analyses/liver/resNT_with_degs.csv",
       sep=";")
res2[is.na(padj),padj:=1]
ggplot(res2)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score)))+scale_y_log10()

ggplot(res2[padj<0.05])+geom_point(aes(x=abs(log2FoldChange),y=abs(gene_score)))


#2021-01-11
#Liver with n.cpg = 1 without eQTR, est ce que ça change la correl ?

library(data.table)
library(ggplot2)
library(patchwork)
res_n<-fread("analyses/liver/res_gene_score_HCCvsCtrl_liver.csv")
cpgs_score<-fread("analyses/liver/cpgs_genes_annot_and_weight_reference.tsv")

res_n_sans_4genes<-merge(unique(res_n[,.(cpgID,pval,meth.change)]),
                       unique(cpgs_score[in_eQTR==F][,is.closest:=abs(tss_dist)==min(abs(tss_dist)),
                                                     by=.(cpgID)][is.closest==T][order(cpgID,gene,-links_score)],
                              by=c("cpgID","gene")),by=c("cpgID"))

res_n_sans_4genes[,cpg_score:=(-log10(pval)/4*meth.change)*regul_weight*links_weight]

res_n_sans_4genes[,n_cpg.gene:=.N,by=.(gene)]
res_n_sans_4genes[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.25,by="gene"]
res_n_sans_4genes[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

ggplot(unique(res_n_sans_4genes,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
degs<-fread("analyses/liver/res_with_degs.csv")
degs<-unique(degs,by="gene")
degs[is.na(padj),padj:=1]
res_n_sans_4genes_degs<-merge(unique(res_n_sans_4genes[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)],by="gene")

pnew<-ggplot(res_n_sans_4genes_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(0,0.30))
pnew
res_n_sans_4genes_degs[padj<0.05]#812 DEGs
wilcox.test(res_n_sans_4genes_degs[padj<0.05]$gene_score,res_n_sans_4genes_degs[padj>=0.05]$gene_score) #0.0001757


pnew2<-ggplot(res_n_sans_4genes_degs)+geom_boxplot(aes(x=padj<0.01,y=abs(gene_score),fill=padj<0.01),outlier.shape = NA)+coord_cartesian(ylim = c(0,0.30))
pnew2
res_n_sans_4genes_degs[padj<0.01]#217 DEGs
wilcox.test(res_n_sans_4genes_degs[padj<0.01]$gene_score,res_n_sans_4genes_degs[padj>=0.01]$gene_score) #0.006575

pnew+pnew2

fwrite(res_n_sans_4genes_degs,"analyses/liver/2021-01-11_res_1gene_by_cpg_bonne_correl_with_degs.csv",sep=";")
