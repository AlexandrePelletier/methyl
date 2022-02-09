#make CpG matrixes Placenta, Pancreas (and put liver in same format if necessary)
library(data.table)
library(stringr)
library(limma)
library(ggplot2)
source("scripts/utils/new_utils.R")
#PLACENTA : 
out<-"analyses/placenta"
dir.create(out)
#methylation data
#I) Get DMC 
#1)preprocess
#methylation data
meth<-fread("datasets/data_placenta/run44_NICHD_ctrl-bkg.txt")
dim(meth)
head(meth[,1:10])

cols<-c("TargetID",colnames(meth)[str_detect(colnames(meth),"AVG_Beta")])

meth<-meth[,.SD,.SDcols=cols]

meth2<-fread("datasets/data_placenta/run47_NICHD_ctrl_bkg.txt")
cols2<-c("TargetID",colnames(meth2)[str_detect(colnames(meth2),"AVG_Beta")])
meth2<-meth2[,.SD,.SDcols=cols2]
samples1<-str_remove(cols[-1],"\\.AVG_Beta")
samples2<-str_remove(cols2[-1],"\\.AVG_Beta")

meth<-merge(meth,meth2,by="TargetID")
meth[,1:10]
meth<-meth[,cpg_id:=TargetID][,-"TargetID"]
meth<-meth[,.SD,.SDcols=c(ncol(meth),1:(ncol(meth)-1))]
colnames(meth)[-1]<-paste0("p",str_remove_all(colnames(meth)[-1],"\\.AVG_Beta|TB"))
meth
ncol(meth)
fwrite(meth,"datasets/data_placenta/312_samples_beta_value_matrix.csv",sep=";")
meth<-fread("datasets/data_placenta/312_samples_beta_value_matrix.csv")
meth<-as.matrix(meth,rownames = "cpg_id")

#metadata samples
mtd<-fread("datasets/data_placenta/batch_outlier_methyl_cohort1_2_040716.txt")
mtd[,batch:=c.rep.1..174...rep.2..129..]
mtd[,sample:=paste0("p",str_remove(Tycko.lab.ID,"TB"))]
mtd[,sex:=as.factor(ifelse(GENDER==1,"M","F"))]
fwrite(mtd,"datasets/data_placenta/batch_outlier_methyl_cohort1_2_040716.txt")
#infos cpgs
anno<-fread("ref/all_anno_450K_hg38.csv")
qual<-fread("ref/HM450.hg38.manifest.tsv.gz")
anno[seqnames=="chrY"]
meth["cg13654344","p8798"] #there is sexual chr to rm

#cpg filtering
cpgs_sex<a
meth_sex<-meth[rownames(meth)%in%anno[!seqnames%in%c("chrY","chrX")]$locisID,]
dim(meth_sex) #474053    312


meth_f<-meth_sex[!rownames(meth_sex)%in%qual[MASK_general==T]$probeID,]
dim(meth_f) #411465    312

#[to do if necessary] preprocess / normalize with ENmix (https://www.bioconductor.org/packages/release/bioc/vignettes/ENmix/inst/doc/ENmix.pdf)
# M<-matrix_for_methylated_intensity
# U<-matrix_for_unmethylated_intensity
# pheno<-as.data.frame(cbind(colnames(M), colnames(M)))
# names(pheno)<-c("Basename","filenames")
# rownames(pheno)<-pheno$Basename
# pheno<-AnnotatedDataFrame(data=pheno)
# anno<-c("IlluminaHumanMethylation450k", "ilmn12.hg19")
# names(anno)<-c("array", "annotation")
# mdat<-MethylSet(Meth = M, Unmeth = U, annotation=anno,phenoData=pheno)

#2) limma compa by sex

design<-model.matrix(~0+sex,data = data.frame(mtd,row.names = "sample"))
head(design)
fit <- lmFit(meth_f[,rownames(design)], design)

colnames(design)
cont.matrix <- makeContrasts(sex= "sexF-sexM",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
res<-topTable(fit2,coef = "sex",n = Inf)
res<-data.table(res,keep.rownames = "cpg_id")
res<-res[,pval:=P.Value][,pval_adj:=adj.P.Val][,meth.change:=logFC][,avg_pct.meth:=AveExpr][,.(cpg_id,meth.change,pval,pval_adj,avg_pct.meth)]
res[pval_adj<0.05] #17117 / 411465
#save
ggplot(res[!is.na(pval)])+geom_point(aes(x=meth.change,y=-log10(pval),col=abs(meth.change)>0.10&pval_adj<0.05))
ggsave("analyses/placenta/volcano_FvsM.png")
fwrite(res,"analyses/placenta/res_FvsM.csv",sep=";")


#[bonus] by iugr like vs not
mtd[,iugr_like:=ifelse(IUGR.status..Y.N.=="Yes",T,F)]
design<-model.matrix(~0+iugr_like+sex,data = data.frame(mtd,row.names = "sample"))
head(design)
fit <- lmFit(meth_f[,rownames(design)], design)

colnames(design)
cont.matrix <- makeContrasts(iugr_status= "iugr_likeTRUE-iugr_likeFALSE",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
res<-topTable(fit2,coef = "iugr_status",n = Inf)
res<-data.table(res,keep.rownames = "cpg_id")
res<-res[,pval:=P.Value][,pval_adj:=adj.P.Val][,meth.change:=logFC][,avg_pct.meth:=AveExpr][,.(cpg_id,meth.change,pval,pval_adj,avg_pct.meth)]
res[pval_adj<0.5] #0 / 411465
#save
ggplot(res[!is.na(pval)])+geom_point(aes(x=meth.change,y=-log10(pval),col=abs(meth.change)>0.02&pval<0.001))
ggsave("analyses/placenta/volcano_iugr_like_vs_not.png")
fwrite(res,"analyses/placenta/res_iugr_like_vs_not.csv",sep=";")
#[end bonus] by iugr like vs not

#II) functional CpGs annotation
anno[,cpg_id:=locisID][,chr:=seqnames][,pos:=start]
fwrite(anno,"ref/all_anno_450K_hg38.csv",sep=";")
cpgs<-anno[,.(chr,pos,cpg_id)][cpg_id%in%c(res$cpg_id)]
rm(res)

#1)annot for cpgs_genes_links
# a)by distance to tss
tss<-fread("ref/hg38/hg38_refseq_curated.txt.gz")
tss<-tss[,chr:=chrom][strand=="+",start:=txStart][strand=="-",start:=txEnd][,end:=start+1][,gene:=name2][,.(chr,start,end,gene,strand)]
fwrite(tss[order(chr,start)],"ref/hg38/tss_genes.bed",sep="\t")

#first genes within +/-200kb

genes_in_200kb<-bed_inter(cpgs[,start:=pos-200e3][start<0,start:=0][,end:=pos+200e3][,.(chr,start,end,cpg_id)][order(chr,start)],
          "ref/hg38/tss_genes.bed",
          select = c(5,6,8,9,4),col.names = c("chr","tss_pos","gene","strand","cpg_id"))

genes_in_200kb
cpgs_genes_tss<-merge(unique(cpgs),unique(genes_in_200kb),all.x=T,by=c("chr","cpg_id"))[,.(chr,pos,cpg_id,gene,tss_pos,strand)]
cpgs_genes_tss[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]

cpgs_genes_1<-cpgs_genes_tss[,closest.gene:=abs(tss_dist)==min(abs(tss_dist)),by=.(cpg_id)][closest.gene==T][,-"closest.gene"]
summary(abs(cpgs_genes_1$tss_dist))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#       0     269    1959   17162   17985  199970


fwrite(cpgs_genes_1,fp(out,"cpgs_closest_gene_tss_within_200kb_around.csv"),sep=";")

# b)by presence in eQTL region
eqtl<-fread("datasets/data_placenta/eqtl_perms_full_090117.txt",select = c(1,7,8,9,10,16,19),col.names = c("refseq_mrna","tss_distance","snp_id","chr","pos.hg19","pval","qval"))

#trans in hg38
fwrite(unique(eqtl[,.(chr,pos.hg19,pos.hg19,snp_id)][order(chr,pos.hg19)][!is.na(chr)]),"temp_placenta_eqtl_hg19.bed",col.names = F,sep="\t")
system("CrossMap.py bed ref/hg19ToHg38.over.chain.gz temp_placenta_eqtl_hg19.bed temp_placenta_eqtl_hg38.bed")
trans<-fread("temp_placenta_eqtl_hg38.bed",select=c(1,2,4),col.names = c("chr","pos","snp_id"))
eqtl<-merge(eqtl,trans[,-"chr"])

#trans nm in gene symbol
refseq_trans<-fread("ref/refseq_mrna_to_hgnc_symbol_translator.csv")
refseq_trans<-refseq_trans[,gene:=hgnc_symbol][,-"hgnc_symbol"]
eqtl<-merge(eqtl,refseq_trans,by="refseq_mrna")
fwrite(eqtl,"ref/eQTL/placenta_eQTL_fabien.csv",sep=";")
eqtl[,start.eQTR:=pos-500]
eqtl[,end.eQTR:=pos+500]
cpgs_eQTR<-bed_inter(a=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b=eqtl[,chr:=paste0("chr",chr)][,.(chr,start.eQTR,end.eQTR,gene)][order(chr,start.eQTR)],
          select = c(4,8),col.names = c("cpg_id","gene"))

cpgs_eQTR #23k  cpgs - snp+/-500pb match
unique(cpgs_eQTR) #12k cpg-gene match
unique(cpgs_eQTR,by="cpg_id")#7k cpgs linked to a gene
unique(cpgs_eQTR,by="gene")#5k genes linked to a cpg
cpgs_eQTR[,in_eQTR:=T]

#merge the 2 links and calculate linksWeight
cpgs_genes_1[,in_eQTR:=F]
cpgs_genes<-merge(cpgs_genes_1[,.(cpg_id,gene,in_eQTR,tss_dist)],cpgs_eQTR[,.SD,.SDcols=c("cpg_id","in_eQTR","gene")],all=T)
cpgs_genes[in_eQTR==T]

cpgs_genes<-merge(cpgs[,.(cpg_id,chr,pos)],cpgs_genes,all.x = T,by="cpg_id")#merge with coord


cpgs_genes<-unique(cpgs_genes[!is.na(gene)&gene!=""])
unique(cpgs_genes,by="cpg_id") #400k/~760k cpgs link to a gene

unique(cpgs_genes[in_eQTR==T],by="cpg_id") #dont 7k linked thx to eQTR
cpgs_genes[,double.linked:=any(in_eQTR==T)&any(in_eQTR==F),by=c("cpg_id","gene")]
unique(cpgs_genes[double.linked==T],by="cpg_id")#dont 1500 also gene linked based on tss

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
          b="ref/ensembl_regulatory/ensembl_regulatory_hg38.bed",
          select = c(4,8),col.names = c("cpg_id","ensembl_regulatory_domain"))



cpgs_ensembl[,ensembl_regulatory_domain:=paste(ensembl_regulatory_domain,collapse = "/"),by="cpg_id"] 
cpgs_ensembl<-unique(cpgs_ensembl)
cpgs_ensembl

# b)chromatin regulatory feature  matching
chrine_feat<-fread("datasets/data_placenta/placenta_7_dense.bed",skip = 1,select = c(1,2,3,4),col.names = c("chr","start","end","chromatin_state"))
chrine_feat[,chrine_id:=1:.N]
#trans in hg38
fwrite(unique(chrine_feat[,.(chr,start,end,chrine_id)][order(chr,start)][!is.na(chr)]),"temp_placenta_chrine_feat_hg19.bed",col.names = F,sep="\t")
system("CrossMap.py bed ref/hg19ToHg38.over.chain.gz temp_placenta_chrine_feat_hg19.bed temp_placenta_chrine_feat_hg38.bed")

trans<-fread("temp_placenta_chrine_feat_hg38.bed",select=c(1,2,3,4),col.names = c("chr","start","end","chrine_id"))
chrine_feat<-merge(chrine_feat[,-c("chr","start","end")],trans)

states<-data.table(chromatin_state=1:7,chromatin_feature=c("Promoter","Active Enhancer",rep("Inactive Enhancer",2),"state5","state6","Gene Body"))
chrine_feat<-merge(chrine_feat,states,by='chromatin_state')
dir.create("ref/Chromatin_Annot/Placenta/")
fwrite(chrine_feat[,.(chr,start,end,chromatin_state,chromatin_feature,chrine_id)],"ref/Chromatin_Annot/Placenta/chromatin_states_hg38.bed",sep="\t")
cpgs_chromatin<-bed_inter(a=unique(cpgs,by="cpg_id")[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          b=chrine_feat[,.(chr,start,end,chromatin_feature)][order(chr,start)],
          select = c(4,8),col.names = c("cpg_id","chromatin_feature"))

cpgs_chromatin
cpgs_reg<-merge(cpgs_chromatin,cpgs_ensembl,all.x=T,by="cpg_id")
cpgs_reg<-merge(cpgs[,.(cpg_id,chr,pos)],cpgs_reg,all.x = T,by="cpg_id")#merge with coord


#c) regulatory_weight
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
source("scripts/utils/new_utils.R")
out<-"analyses/placenta/"
res<-fread("analyses/placenta/res_FvsM.csv")
cpgs_anno<-fread(fp(out,"cpgs_genes_annot_and_weight_reference.tsv"))

res<-merge(res[,.(cpg_id,pval,meth.change)],unique(cpgs_anno[order(cpg_id,gene,-links_score)],by=c("cpg_id","gene")),by=c("cpg_id"))
res[,cpg_score:=(-log10(pval)/4*(meth.change*100))*regul_weight*links_weight]

#GeneScore
res[,n_cpg.gene:=.N,by=.(gene)]
res[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
res[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

plot(density(unique(res,by="gene")$n_cpg.gene))
summary(unique(res,by="gene")$n_cpg.gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.00    4.00   10.00   13.77   18.00  461.00 

plot(density(unique(res,by="gene")$gene_score))

ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
fwrite(res,fp(out,"res_FvsM_genescore.csv"),sep=";")

#IV) DEGs analysis
#reformat expression matrix
expr<-fread("datasets/data_placenta/expression_84samples_Placenta_031715.txt")
expr<-expr[,-c("chr","start","stop")]
colnames(expr)<-c("gene",paste0("p",str_remove(colnames(expr)[-1],"BM")))
expr
ncol(expr)
fwrite(expr,"datasets/data_placenta/84_samples_expression_matrix.csv",sep=";")

#DE analyis
library(DESeq2)
expr<-fread("datasets/data_placenta/84_samples_expression_matrix.csv")
expr[,expr_mean:=rowMeans(as.matrix(expr[,-"gene"]))]
expr<-unique(expr[order(gene,-expr_mean)][,-"expr_mean"],by="gene")
counts<-as.matrix(data.frame(expr,row.names = 1))
head(counts)
dim(counts) #25548    83
samples<-colnames(counts)
mtd<-fread("datasets/data_placenta/batch_outlier_methyl_cohort1_2_040716.txt")
mtd[,sex:=factor(sex,levels = c("M","F"))]
table(mtd$sample%in%samples) #74/83 in common
samples<-samples[samples%in%mtd$sample]
length(samples)
dds1 <- DESeqDataSetFromMatrix(round(counts)[,samples], 
                               colData = data.frame(mtd[,.(sample,sex)],row.names ="sample",drop=F)[samples,], 
                               design = ~ sex )
dds1 <- DESeq(dds1)
resultsNames(dds1)
# resultsNames(dds)
res1 <- results(dds1, 
                contrast = c("sex","F","M"),
                alpha = 0.05)

res1 <- lfcShrink(dds1, 
                  coef = "sex_F_vs_M" ,
                  res=res1)
res1<-data.table(data.frame(res1),keep.rownames = "gene")

res1[padj<0.05]
fwrite(res1,"analyses/placenta/res_degs_FvsM.csv",sep=";")


#V) Validation GeneScore
out<-"analyses/placenta"
res_de<-fread("analyses/placenta/res_degs_FvsM.csv")

res_meth<-fread(fp(out,"res_FvsM_genescore.csv"))
res_merge<-merge(res_meth,res_de,by="gene")


fwrite(res_merge[order(gene)],
       fp(out,"res_FvsM_meth_and_degs_merge.csv"),
       sep=";")
res_merge[is.na(padj),padj:=1]
ggplot(unique(res_merge,by="gene"))+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score)))+scale_y_log10()

ggplot(res_merge[padj<0.05])+geom_point(aes(x=abs(log2FoldChange),y=abs(gene_score)))
