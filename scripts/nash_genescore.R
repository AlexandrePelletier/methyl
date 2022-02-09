#make CpG matrixes Placenta, Pancreas (and put liver in same format if necessary)
library(limma)
source("scripts/utils/new_utils.R")
#PLACENTA : 
out<-"analyses/nash"
dir.create(out)
#methylation data
#I) Get DMC 
res<-fread("../ewas_de_data/PreciNASH_EWAS_DMP_NASH_cell.csv.gz",
            select = c(1,2,5,7,9,10),col.names =c("cpg_id","meth.change","pval","pval_adj","chr","pos") )

fwrite(res,"analyses/nash/res_meth_nash_vs_not.csv.gz")
#II) functional CpGs annotation
cpgs<-res[,.(chr,pos,cpg_id)]
rm(res)

#1)annot for cpgs_genes_links
# a)by distance to tss

#first genes within +/-200kb 
cpgs<-cpgs[!is.na(chr)&chr!='']


genes_in_200kb<-bed_inter(cpgs[,start:=pos-200e3][start<0,start:=0][,end:=pos+200e3][,.(chr,start,end,cpg_id)][order(chr,start)],
          "ref/hg19/tss_genes.bed",
          select = c(5,6,8,9,4),col.names = c("chr","tss_pos","gene","strand","cpg_id"))

genes_in_200kb
cpgs_genes_tss<-merge(unique(cpgs),unique(genes_in_200kb),all.x=T,by=c("chr","cpg_id"))[,.(chr,pos,cpg_id,gene,tss_pos,strand)]
cpgs_genes_tss[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]

cpgs_genes_1<-cpgs_genes_tss[,closest.gene:=abs(tss_dist)==min(abs(tss_dist)),by=.(cpg_id)][closest.gene==T][,-"closest.gene"]
summary(abs(cpgs_genes_1$tss_dist))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #     0     566    6800   24239   31640  199983 


fwrite(cpgs_genes_1,fp(out,"cpgs_closest_gene_tss_within_200kb_around.csv"),sep=";")


#by presence in eQTL region

cpgs_eQTR<-bed_inter(a=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b="ref/eQTL/liver_eQTR_hgnc_gene_hg19.bed",
          select = c(4,8),col.names = c("cpg_id","gene"))

cpgs_eQTR[,in_eQTR:=T]

cpgs_eQTR #1M5 cpgs - snp+/-500pb match
unique(cpgs_eQTR) #800k cpg-gene match
unique(cpgs_eQTR,by="cpg_id")#300k cpgs linked to a gene
unique(cpgs_eQTR,by="gene")#13k genes linked to a cpg


#merge the 2 links
cpgs_genes_1[,in_eQTR:=F]
cpgs_genes<-merge(cpgs_genes_1[,.(cpg_id,gene,in_eQTR,tss_dist)],unique(cpgs_eQTR),all=T)
cpgs_genes[in_eQTR==T]


cpgs_genes<-merge(cpgs[,.(cpg_id,chr,pos)],cpgs_genes,all.x = T,by="cpg_id")#merge with coord


cpgs_genes<-unique(cpgs_genes[!is.na(gene)&gene!=""])
unique(cpgs_genes,by="cpg_id") #716k/~760k cpgs link to a gene

cpgs_genes[,double.linked:=any(in_eQTR==T)&any(in_eQTR==F),by=c("cpg_id","gene")]
unique(cpgs_genes[double.linked==T],by="cpg_id")#dont 36k linked by tss and eqtl

cpgs_genes[in_eQTR==F,links_score:=sapply(abs(tss_dist),function(x){
  if(x<1000)return(1)
  else if(x<20000)return(0.5+0.5*sqrt(1000/x))
  else return(0.5*sqrt(20000/x))
    })]

cpgs_genes[in_eQTR==T,links_score:=1]

cpgs_genes[,links_weight:=max(links_score),by=c("cpg_id","gene")]

fwrite(cpgs_genes,fp(out,"cpgs_closest_gene_tss_linked_within_200kb_around_and_eQTR_linked_genes.tsv"),sep="\t")

#2)annot for regulatory region
# a)ensembl regulatory domain matching

cpgs_ensembl<-bed_inter(a=unique(cpgs,by="cpg_id")[,chr:=str_remove(chr,'chr')][,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b="ref/ensembl_regulatory/ensembl_regulatory_hg19.bed",
          select = c(4,8),col.names = c("cpg_id","ensembl_regulatory_domain"))



cpgs_ensembl[,ensembl_regulatory_domain:=paste(ensembl_regulatory_domain,collapse = "/"),by="cpg_id"] 
cpgs_ensembl<-unique(cpgs_ensembl)
cpgs_ensembl

# b)chromatin regulatory feature  matching
chrine_feat<-fread("ref/Chromatin_Annot/Liver/chromatin_features_reg.bed")

cpgs_chromatin<-bed_inter(a=unique(cpgs,by="cpg_id")[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b=chrine_feat[,.(chr,start,end,feature)][order(chr,start)],
          select = c(4,8),col.names = c("cpg_id","chromatin_feature"))

cpgs_chromatin
cpgs_reg<-merge(cpgs_chromatin,cpgs_ensembl,all.x=T,by="cpg_id")
cpgs_reg<-merge(cpgs[,.(cpg_id,chr,pos)],cpgs_reg,all.x = T,by="cpg_id")#merge with coord

# c) regulatory_weight
cpgs_reg[is.na(chromatin_feature),chromatin_score:=0]
cpgs_reg[!is.na(chromatin_feature),chromatin_score:=sapply(chromatin_feature,function(x){
  
      if(x%in%c("Promoter","Active Enhancer"))return(1)
      else if(x=="Inactive Enhancer")return(0.75)
      else if(x%in%"Gene Body")return(0.5)
      else return(0)
    })]


cpgs_reg[,ensembl_reg_score:=sapply(ensembl_regulatory_domain,function(x){
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

cpgs_reg[,regul_weight:=(0.5+1.5*((chromatin_score+ensembl_reg_score)/2))]

fwrite(cpgs_reg,fp(out,"cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.tsv"),sep="\t")

#3) merge and save cpgs anno and weight
cpgs_anno<-merge(cpgs_genes,cpgs_reg[,-c("chr",'pos')],all.x=T,by=c("cpg_id"))


fwrite(cpgs_anno,fp(out,"cpgs_genes_annot_and_weight_reference.tsv"),sep="\t")

#III) calculate Gene Score
#CpG score
res<-fread("analyses/nash/res_meth_nash_vs_not.csv.gz")
cpgs_anno<-fread(fp(out,"cpgs_genes_annot_and_weight_reference.tsv"))

res<-merge(res[,.(cpg_id,pval,meth.change)],unique(cpgs_anno[order(cpg_id,gene,-links_score)],by=c("cpg_id","gene")),by=c("cpg_id"))
res[,cpg_score:=(-log10(pval)/4*(meth.change*100))*regul_weight*links_weight]

#GeneScore
res[,n_cpg.gene:=.N,by=.(gene)]
plot(density(unique(res,by="gene")$n_cpg.gene))
summary(unique(res,by="gene")$n_cpg.gene)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #   1.00   16.00   33.00   56.37   62.00 2894.00 

res[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.4,by="gene"]
res[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

plot(density(unique(res,by="gene")$gene_score))

ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
fwrite(res,fp(out,"res_nash_vs_not_genescore.csv"),sep=";")

#IV) DEGs analysis
#reformat expression matrix
res_de<-fread("../ewas_de_data/PreciNASH_DE_NASH_gene.csv.gz",select = c(12,2,3,6,7),
            col.names = c("gene","baseMean","log2FoldChange","pvalue","padj"))

res_de[padj<0.001]
fwrite(res_de,fp(out,"res_degs_nash_vs_not.csv"),sep=";")


#V) Validation GeneScore
res_meth<-fread(fp(out,"res_nash_vs_not_genescore.csv"))
res_merge<-merge(res_meth,res_de,by="gene")


fwrite(res_merge[order(gene)],
       fp(out,"res_nash_vs_not_meth_and_degs_merge.csv"),
       sep=";")
res_merge[is.na(padj),padj:=1]
ggplot(unique(res_merge[in_eQTR==F],by="gene"))+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score)))+scale_y_log10()

ggplot(res_merge[padj<0.05])+geom_point(aes(x=abs(log2FoldChange),y=abs(gene_score)))

