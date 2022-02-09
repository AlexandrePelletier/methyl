

#validation correl FeatureRegion meth.change with Fold change of genes +/-5kb.
#res attendu : promoter (ctl +), enhancer (?), in other region (ctl -), in cbps, liver, placenta
source("scripts/utils/new_utils.R")
out<-"analyses/genescore_by_feature_region/"
dir.create(out)
# FOR CBPs
# - summarize Meth change by chromatin feature region :
regions<-fread("ref/Chromatin_Annot/CBP/CD34_all_chromatin_feature.csv",select=c(1:3,5,4),
               col.names=c("chr","start",
                          "end","length",
                          "chromatin_state"
                          ))

regions[,region_id:=1:.N]
chrine_feat<-data.table(chromatin_state=0:6,chromatin_feature=c("Heterochromatin",rep("Gene Body",3),
                                                                "Active Enhancer","Poised Enhancer",
                                                                "Promoter"))
regions<-merge(regions,chrine_feat,by="chromatin_state")


#overlap cpgs  and regions
res<-fread("analyses/2020-12-08_cbps_new_genescore/res_gene_score_lgaF_vs_ctrlF.csv")
res<-res[,cpg_id:=cpgID][,-"cpgID"]
cpgs<-unique(res[,.(chr,pos,cpg_id)])
inter_reg_cpgs<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","cpg_id"))

regions_cpgs<-merge(regions,inter_reg_cpgs,all.x=T,by="region_id")

# add meth change info
regions_cpgs_meth<-merge(regions_cpgs,unique(res[,.(cpg_id,meth.change,pval)]),by="cpg_id")
regions_cpgs_meth[pval<0.001]
regions_cpgs_meth[,avg.pval.region:=mean(pval),by="region_id"]
regions_cpgs_meth[,avg.mlog10.pval.region:=mean(-log10(pval)),by="region_id"]
regions_cpgs_meth[,avg.meth.change.region:=mean(meth.change),by="region_id"]
regions_cpgs_meth[,dmc_score:=-log10(pval)/4*meth.change]
regions_cpgs_meth[,avg.dmc_score:=mean(dmc_score),by="region_id"]
regions_cpgs_meth[,max.dmc_score:=dmc_score[which.max(abs(dmc_score))],by="region_id"]
regions_cpgs_meth[,med.dmc_score:=median(dmc_score),by="region_id"]

#gene_score based
regions_cpgs_meth[,sum.dmc_score:=sum(dmc_score),by="region_id"]

regions_cpgs_meth[,n_cpg.reg:=.N,by=.(region_id)]
summary(regions_cpgs_meth$n_cpg.reg)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 1.00    2.00    4.00   10.39   12.00  420.00 


regions_cpgs_meth[,n_cpg_weight:=(1/sum(1/(abs(dmc_score)+1)))^0.25,by="region_id"]

regions_cpgs_meth[,gs_based.dmc_score:=sum(dmc_score)*n_cpg_weight,by="region_id"]
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=sum.dmc_score))
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=gs_based.dmc_score))
fwrite(regions_cpgs_meth,fp(out,"cbps_regions_meth.change.csv.gz"),sep=";")

#++ meth change on enhancer and promoter :
ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=avg.pval.region,col=chromatin_feature))

ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=avg.mlog10.pval.region,col=chromatin_feature))

ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change.region,col=chromatin_feature))

ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=avg.dmc_score,col=chromatin_feature))

#add Fold change of genes +/-5kb
tss_genes<-fread("ref/hg19/hg19_refseq_curated.txt.gz",select = c(2,3,4,5,6,13),
                 col.names = c("gene_id","chr","strand","start","end","gene"))
tss_genes[strand=="+",tss_pos:=start]
tss_genes[strand=="-",tss_pos:=end]
tss_genes<-tss_genes[chr%in%paste0("chr",c(1:21,"X","Y","MT"))]
fwrite(tss_genes[,.(gene_id,chr,tss_pos,strand,gene)],"ref/hg19/tss_pos_hg19_refseq_curated_250121.txt.gz")
tss_genes[,start:=tss_pos-5000][start<0,start:=0][,end:=tss_pos+5000]
tss_genes[order(chr,start)]
inter_reg_genes<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=tss_genes[,.(chr,start,end,gene_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","gene_id"))
regions_cpgs_meth_genes<-merge(regions_cpgs_meth,inter_reg_genes,by="region_id")

expr_change<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-08-19_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
expr_change<-merge(expr_change,unique(tss_genes[,.(gene_id,gene,tss_pos)]))[,.(gene_id,gene,tss_pos,baseMean,log2FoldChange,pvalue,padj)]
expr_change[is.na(padj),padj:=1]
regions_cpgs_meth_genes_expr<-merge(regions_cpgs_meth_genes,expr_change,by="gene_id")

fwrite(regions_cpgs_meth_genes_expr,fp(out,"cpbs_regions_meth.change_degs.csv.gz"),sep=";")


#correl FeatureRegion meth.change with Fold change
#boxplot meth.change for DEGs
ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.pval.region),fill=chromatin_feature))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.mlog10.pval.region),fill=chromatin_feature))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.meth.change.region),fill=chromatin_feature))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.dmc_score),fill=chromatin_feature))



#boxplot degs vs not 
ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.meth.change.region),col=padj<0.05))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.pval.region),col=padj<0.05))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.dmc_score),col=padj<0.05))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(max.dmc_score),col=padj<0.05))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(med.dmc_score),col=padj<0.05))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(gs_based.dmc_score),col=padj<0.2))

length(unique(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05]$gene))#102

#geom_point
ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_point(aes(x=log2FoldChange,y=gs_based.dmc_score,col=chromatin_feature))+facet_wrap("chromatin_feature")


ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_point(aes(x=log2FoldChange,y=max.dmc_score,col=chromatin_feature))

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_point(aes(x=log2FoldChange,y=max.dmc_score,col=chromatin_feature))+facet_wrap("chromatin_feature")

ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id"))[padj<0.05])+
  geom_point(aes(x=chromatin_feature,y=abs(max.dmc_score),col=padj<0.1))



#SAME FOR LIVER PLACENTA, NASH

#LIVER
# - summurarize Meth change by chromatin feature region :
regions<-fread("ref/Chromatin_Annot/Liver/chromatin_features_reg.bed",
               col.names=c("chr","start",
                          "end","chromatin_state",
                          "chromatin_feature"))

regions<-regions[order(chr,start)][,region_id:=1:.N]
#overlap cpgs  and regions
res<-fread("analyses/liver/res_gene_score_HCCvsCtrl_liver.csv")
res<-res[,cpg_id:=cpgID][,-"cpgID"]
cpgs<-unique(res[,.(chr,pos,cpg_id)])
inter_reg_cpgs<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","cpg_id"))

regions_cpgs<-merge(regions,inter_reg_cpgs,all.x=T,by="region_id")

# add meth change info
regions_cpgs_meth<-merge(regions_cpgs,unique(res[,.(cpg_id,meth.change,pval)]),by="cpg_id")
regions_cpgs_meth[pval<0.001]
regions_cpgs_meth[,avg.pval.region:=mean(pval),by="region_id"]
regions_cpgs_meth[,avg.mlog10.pval.region:=mean(-log10(pval)),by="region_id"]
regions_cpgs_meth[,avg.meth.change.region:=mean(meth.change),by="region_id"]
regions_cpgs_meth[,dmc_score:=-log10(pval)/4*meth.change*100]
regions_cpgs_meth[,avg.dmc_score:=mean(dmc_score),by="region_id"]
regions_cpgs_meth[,max.dmc_score:=dmc_score[which.max(abs(dmc_score))],by="region_id"]
regions_cpgs_meth[,med.dmc_score:=median(dmc_score),by="region_id"]

#gene_score based
regions_cpgs_meth[,n_cpg.reg:=.N,by=.(region_id)]
summary(regions_cpgs_meth$n_cpg.reg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1.00    8.00   18.00   34.78   36.00 1944.00
regions_cpgs_meth[,n_cpg_weight:=(1/sum(1/(abs(dmc_score)+1)))^0.25,by="region_id"]
regions_cpgs_meth[,gs_based.dmc_score:=sum(dmc_score)*n_cpg_weight,by="region_id"]
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=gs_based.dmc_score))
fwrite(regions_cpgs_meth,fp(out,"liver_regions_meth.change.csv.gz"),sep=";")


#add Fold change of genes +/-5kb
tss_genes<-fread("ref/hg19/tss_pos_hg19_refseq_curated_250121.txt.gz")
tss_genes[,start:=tss_pos-5000][start<0,start:=0][,end:=tss_pos+5000]
tss_genes[order(chr,start)]
inter_reg_genes<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=tss_genes[,.(chr,start,end,gene_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","gene_id"))
regions_cpgs_meth_genes<-merge(regions_cpgs_meth,inter_reg_genes,by="region_id",allow.cartesian=TRUE)

expr_change<-unique(fread("analyses/liver/2021-01-11_res_1gene_by_cpg_bonne_correl_with_degs.csv"),by="gene")
expr_change<-merge(expr_change,unique(tss_genes[,.(gene_id,gene,tss_pos)]))[,.(gene_id,gene,tss_pos,log2FoldChange,pvalue,padj)]
expr_change[is.na(padj),padj:=1]
regions_cpgs_meth_genes_expr<-merge(regions_cpgs_meth_genes,expr_change,by="gene_id")

ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.dmc_score),col=chromatin_feature))+scale_y_log10()
ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(gs_based.dmc_score),col=padj<0.2))+scale_y_log10()

fwrite(regions_cpgs_meth_genes_expr,fp(out,"liver_regions_meth.change_degs.csv.gz"),sep=";")


#PLACENTA
# - summurarize Meth change by chromatin feature region :
regions<-fread("ref/Chromatin_Annot/Placenta/chromatin_states_hg38.bed")

regions<-regions[order(chr,start)][,region_id:=1:.N][,-"chrine_id"]
#overlap cpgs  and regions
res<-fread("analyses/placenta/res_FvsM_genescore.csv")

cpgs<-unique(res[,.(chr,pos,cpg_id)])
inter_reg_cpgs<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","cpg_id"))

regions_cpgs<-merge(regions,inter_reg_cpgs,all.x=T,by="region_id")

# add meth change info
regions_cpgs_meth<-merge(regions_cpgs,unique(res[,.(cpg_id,meth.change,pval)]),by="cpg_id")
regions_cpgs_meth[pval<0.001]
regions_cpgs_meth[,avg.pval.region:=mean(pval),by="region_id"]
regions_cpgs_meth[,avg.mlog10.pval.region:=mean(-log10(pval)),by="region_id"]
regions_cpgs_meth[,avg.meth.change.region:=mean(meth.change),by="region_id"]
regions_cpgs_meth[,dmc_score:=-log10(pval)/4*meth.change*100]
regions_cpgs_meth[,avg.dmc_score:=mean(dmc_score),by="region_id"]
regions_cpgs_meth[,max.dmc_score:=dmc_score[which.max(abs(dmc_score))],by="region_id"]
regions_cpgs_meth[,med.dmc_score:=median(dmc_score),by="region_id"]

#gene_score based
regions_cpgs_meth[,n_cpg.reg:=.N,by=.(region_id)]
summary(regions_cpgs_meth$n_cpg.reg)
 #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #    1.0     3.0     7.0    15.1    14.0   305.0 

regions_cpgs_meth[,n_cpg_weight:=(1/sum(1/(abs(dmc_score)+1)))^0.2,by="region_id"]
regions_cpgs_meth[,sum.dmc_score:=sum(dmc_score),by="region_id"]
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=sum.dmc_score))

regions_cpgs_meth[,length:=end-start]
regions_cpgs_meth[n_cpg.reg>200&sum.dmc_score<(-250)][order(pval)]
regions_cpgs_meth[,gs_based.dmc_score:=sum(dmc_score)*n_cpg_weight,by="region_id"]
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=gs_based.dmc_score))
fwrite(regions_cpgs_meth,fp(out,"placenta_regions_meth.change.csv.gz"),sep=";")


#add Fold change of genes +/-5kb
tss_genes<-fread("ref/hg38/hg38_refseq_curated.txt.gz",select = c(2,3,4,5,6,13),
                 col.names = c("gene_id","chr","strand","start","end","gene"))
tss_genes[strand=="+",tss_pos:=start]
tss_genes[strand=="-",tss_pos:=end]
tss_genes<-tss_genes[chr%in%paste0("chr",c(1:21,"X","Y","MT"))]
fwrite(tss_genes[,.(gene_id,chr,tss_pos,strand,gene)],"ref/hg38/tss_pos_hg38_refseq_curated_250121.txt.gz")
tss_genes[,start:=tss_pos-5000][start<0,start:=0][,end:=tss_pos+5000]
tss_genes[order(chr,start)]
inter_reg_genes<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=tss_genes[,.(chr,start,end,gene_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","gene_id"))
regions_cpgs_meth_genes<-merge(regions_cpgs_meth,inter_reg_genes,by="region_id",allow.cartesian=TRUE)

expr_change<-unique(fread("analyses/placenta/res_degs_FvsM.csv"),by="gene")
expr_change<-merge(expr_change,unique(tss_genes[,.(gene_id,gene,tss_pos)]))[,.(gene_id,gene,tss_pos,baseMean,log2FoldChange,pvalue,padj)]
expr_change[is.na(padj),padj:=1]
regions_cpgs_meth_genes_expr<-merge(regions_cpgs_meth_genes,expr_change,by="gene_id")

ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.meth.change.region),col=chromatin_feature))+scale_y_log10()
ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(sum.dmc_score),col=padj<0.2))+scale_y_log10()

fwrite(regions_cpgs_meth_genes_expr,fp(out,"placenta_regions_meth.change_degs.csv.gz"),sep=";")



#NASH
# - summurarize Meth change by chromatin feature region :
regions<-fread("ref/Chromatin_Annot/Liver/chromatin_features_reg.bed",
               col.names=c("chr","start",
                          "end","chromatin_state",
                          "chromatin_feature"))

regions<-regions[order(chr,start)][,region_id:=1:.N]

#overlap cpgs  and regions
res<-fread("analyses/nash/res_meth_nash_vs_not.csv.gz")

cpgs<-unique(res[,.(chr,pos,cpg_id)])[!(is.na(chr)|chr=="")]
inter_reg_cpgs<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","cpg_id"))

regions_cpgs<-merge(regions,inter_reg_cpgs,all.x=T,by="region_id")

# add meth change info
regions_cpgs_meth<-merge(regions_cpgs,unique(res[,.(cpg_id,meth.change,pval)]),by="cpg_id")
regions_cpgs_meth[pval<0.001]
regions_cpgs_meth[,avg.pval.region:=mean(pval),by="region_id"]
regions_cpgs_meth[,avg.mlog10.pval.region:=mean(-log10(pval)),by="region_id"]
regions_cpgs_meth[,avg.meth.change.region:=mean(meth.change),by="region_id"]
regions_cpgs_meth[,dmc_score:=-log10(pval)/4*meth.change*100]
regions_cpgs_meth[,avg.dmc_score:=mean(dmc_score),by="region_id"]
regions_cpgs_meth[,max.dmc_score:=dmc_score[which.max(abs(dmc_score))],by="region_id"]
regions_cpgs_meth[,med.dmc_score:=median(dmc_score),by="region_id"]

#gene_score based
regions_cpgs_meth[,n_cpg.reg:=.N,by=.(region_id)]
summary(regions_cpgs_meth$n_cpg.reg)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #    1.0     4.0     9.0    12.4    16.0   194.0 

regions_cpgs_meth[,n_cpg_weight:=(1/sum(1/(abs(dmc_score)+1)))^0.2,by="region_id"]
regions_cpgs_meth[,sum.dmc_score:=sum(dmc_score),by="region_id"]
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=sum.dmc_score))

regions_cpgs_meth[,length:=end-start]
regions_cpgs_meth[n_cpg.reg>200&sum.dmc_score<(-250)][order(pval)]
regions_cpgs_meth[,gs_based.dmc_score:=sum(dmc_score)*n_cpg_weight,by="region_id"]
ggplot(unique(regions_cpgs_meth,by="region_id"))+geom_point(aes(x=n_cpg.reg,y=gs_based.dmc_score))
fwrite(regions_cpgs_meth,fp(out,"nash_regions_meth.change.csv.gz"),sep=";")

#add Fold change of genes +/-5kb
tss_genes<-fread("ref/hg19/tss_pos_hg19_refseq_curated_250121.txt.gz")
tss_genes[,start:=tss_pos-5000][start<0,start:=0][,end:=tss_pos+5000]
tss_genes[order(chr,start)]
inter_reg_genes<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
          b=tss_genes[,.(chr,start,end,gene_id)][order(chr,start)],
          select = c(4,8),col.names = c("region_id","gene_id"))
regions_cpgs_meth_genes<-merge(regions_cpgs_meth,inter_reg_genes,by="region_id",allow.cartesian=TRUE)

expr_change<-unique(fread("analyses/nash/res_degs_nash_vs_not.csv"),by="gene")
expr_change<-merge(expr_change,unique(tss_genes[,.(gene_id,gene,tss_pos)]))[,.(gene_id,gene,tss_pos,baseMean,log2FoldChange,pvalue,padj)]
expr_change[is.na(padj),padj:=1]
regions_cpgs_meth_genes_expr<-merge(regions_cpgs_meth_genes,expr_change,by="gene_id")

ggplot(unique(regions_cpgs_meth,by="region_id"))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(avg.meth.change.region),col=chromatin_feature))+scale_y_log10()
ggplot(unique(regions_cpgs_meth_genes_expr,by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=abs(gs_based.dmc_score),col=padj<0.05))+scale_y_log10()

fwrite(regions_cpgs_meth_genes_expr,fp(out,"nash_regions_meth.change_degs.csv.gz"),sep=";")

#ccl : focus on cbps car others doesn't have correlation in promoter methylation / gene expression

#test correl enhancer FeatureRegion with Fold change of genes +/-100 kb. res attendu :  promoter (correl baisse), enhancer (correl monte), in other region (no correl), in cbps, liver, placenta
#run GSR-1



#test correl enhancer FeatureRegion with Fold change of genes in eQTR. res attendu : enhancer (correl monte++), in other region (no correl), in cbps, liver, placenta

#find most signif eQTL locally
eqtls<-list()
#liver
eqtls[["liver"]]<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz")
nrow(eqtls[["liver"]]) #629559 snp-gene

  # extract chr and pos hg38 of snps
eqtls[["liver"]][,chr.hg38:=str_extract(variant_id,"^chr[0-9XYM]{1,2}") ]
eqtls[["liver"]]<-eqtls[["liver"]][!is.na(chr.hg38)][,pval_link:=pval_nominal][,pval_type:="nominal"]

eqtls[["liver"]][,pos.hg38:=as.numeric(str_sub(str_extract(variant_id,"_[0-9]+"),2)) ]
eqtls[["liver"]]

  # calculate dist snps asso to a same gene")
GetDistToNext<-function(positions){
  return(c(sapply(1:(length(positions)-1), function(i)abs(positions[i]-positions[i+1])),NA))
}

eqtls[["liver"]][,dist_to_next:=GetDistToNext(pos.hg38),by="gene_id"]

summary(eqtls[["liver"]]$dist_to_next)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  #     0     114     414    2715    1283 1862634    5734 


  # search most signif eQTL locally (at +/-10kb)
eqtls[["liver"]][,minLocal.cand:=pval_link<quantile(pval_link,0.25),by=c("gene_id")]
eqtls[["liver"]][minLocal.cand==T]#128k/629k

eqtls[["liver"]][,n_eqtl.gene:=.N,by="gene_id"]
summary(unique(eqtls[["liver"]],by="gene_id")$n_eqtl.gene)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #   1.0     8.0    38.0   109.8   115.0  4120.0 

is.minLocal<-function(dists,pvals){
  isMins<-sapply(1:length(dists), function(i){
    locisA5kb<-which(dists>(dists[i]-5000)&dists<(dists[i]+5000))
    
    if(all(pvals[locisA5kb]>=pvals[i])){
      return(T)
    }else{
      return(F)
    }
  })
  
  return(isMins)
}

eqtls[["liver"]][minLocal.cand==T,minLocal:=is.minLocal(tss_distance,pval_link),by=c("gene_id")]
eqtls[["liver"]][minLocal==T] #53901/ 629559

  # instead of 
egenes_liver<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz")
egenes_liver#22k cpgs




#whole blood
tissue<-"whole_blood"

eqtls[[tissue]]<-fread(paste0("ref/eQTL/GTEx_Analysis_v8_eQTL/",str_replace_all(str_to_title(str_replace_all(tissue,"_"," "))," ","_"),".v8.signif_variant_gene_pairs.txt.gz"))
nrow(eqtls[[tissue]]) #2 414 653 snp-gene

  # extract chr and pos hg38 of snps
eqtls[[tissue]][,chr.hg38:=str_extract(variant_id,"^chr[0-9XYM]{1,2}") ]
eqtls[[tissue]]<-eqtls[[tissue]][!is.na(chr.hg38),][,pval_link:=pval_nominal][,pval_type:="nominal"]

eqtls[[tissue]][,pos.hg38:=as.numeric(str_sub(str_extract(variant_id,"_[0-9]+"),2)) ]
eqtls[[tissue]]

  # calculate dist snps asso to a same gene")
eqtls[[tissue]][,dist_to_next:=GetDistToNext(pos.hg38),by="gene_id"]
summary(eqtls[[tissue]]$dist_to_next)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   #    0     113     394    2679    1153 1852417   12360 

# search most signif eQTL locally (at +/-10kb)
eqtls[[tissue]][,minLocal.cand:=pval_link<quantile(pval_link,0.25),by=c("gene_id")]
eqtls[[tissue]][minLocal.cand==T]#570k/2.4M

eqtls[[tissue]][,n_eqtl.gene:=.N,by="gene_id"]
summary(unique(eqtls[[tissue]],by="gene_id")$n_eqtl.gene)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #    1.0    22.0    87.0   195.4   230.0  8571.0 

eqtls[[tissue]][minLocal.cand==T,minLocal:=is.minLocal(tss_distance,pval_link),by=c("gene_id")]
eqtls[[tissue]][minLocal==T] #160k/2.4M


#meta
tissue<-"tissue_wide"
meta<-fread(paste0("ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz"))[PVALUE_Q>0.1&n.tissue_sig>30&PVALUE_FE<10e-30]
meta #552921 snp-gene 
meta<-meta[,pos.hg38:=pos][,chr.hg38:=chr][,.(variant_id,gene_id,chr.hg38,pos.hg38,PVALUE_Q,n.tissue_sig,PVALUE_FE)]
eqtls[[tissue]]<-meta
eqtls[[tissue]][,pval_link:=PVALUE_FE][,pval_type:="PVALUE_FE"]
  # calculate dist snps asso to a same gene")
eqtls[[tissue]][,dist_to_next:=GetDistToNext(pos.hg38),by="gene_id"]
summary(eqtls[[tissue]]$dist_to_next)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
 #      0     184     751    5789    2655 1940632    7246 

# search most signif eQTL locally (at +/-10kb)
eqtls[[tissue]][,minLocal.cand:=pval_link<quantile(pval_link,0.25),by=c("gene_id")]
eqtls[[tissue]][minLocal.cand==T]#130k/550k

eqtls[[tissue]][,n_eqtl.gene:=.N,by="gene_id"]
summary(unique(eqtls[[tissue]],by="gene_id")$n_eqtl.gene)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #   1.00    4.00   16.00   76.31   63.00 5137.00 

eqtls[[tissue]][minLocal.cand==T,minLocal:=is.minLocal(pos.hg38,pval_link),by=c("gene_id")]
eqtls[[tissue]][minLocal==T] #58k/550k

#placenta
tissue<-"placenta"
eqtls[[tissue]]<-fread("ref/eQTL/placenta_eQTL_fabien.csv")
nrow(eqtls[[tissue]]) #26k snp-gene
eqtls[[tissue]][,pos.hg38:=pos][,chr.hg38:=paste0("chr",chr)][,pval_link:=pval][,pval_type:="nominal"]
trans<-TransNMtoEnsemblVers(eqtls$placenta$refseq_mrna)
trans<-unique(trans,by="refseq_mrna")
eqtls[[tissue]]<-merge(eqtls[[tissue]],trans[,gene_id:=ensembl_gene_id_version][,-"ensembl_gene_id_version"])
eqtls[[tissue]][,variant_id:=snp_id]
  # calculate dist snps asso to a same gene")
eqtls[[tissue]][,dist_to_next:=GetDistToNext(pos.hg38),by="gene_id"]
summary(eqtls[[tissue]]$dist_to_next)
    # Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
    #     0         0         0     85620         0 179625114     16845 
  # => very few, keep all snps-gene
eqtls[[tissue]][,minLocal:=T]

#merge all eQTL
eqtls_dt<-Reduce(rbind,lapply(names(eqtls),function(tiss)eqtls[[tiss]][,tissue:=tiss][,.(variant_id,gene_id,chr.hg38,pos.hg38,tissue,pval_link,pval_type,minLocal)]))

#trans in hg19
fwrite(unique(eqtls_dt[,.(chr.hg38,pos.hg38,pos.hg38,variant_id)][order(chr.hg38,pos.hg38)][!is.na(chr.hg38)]),"temp_eqtls_hg38.bed",col.names = F,sep="\t")
system("CrossMap.py bed ref/hg19ToHg38.over.chain.gz temp_eqtls_hg38.bed temp_eqtls_hg19.bed")
trans<-fread("temp_eqtls_hg19.bed",select=c(1,2,4),col.names = c("chr.hg19","pos.hg19","variant_id"))
eqtls_dt<-merge(eqtls_dt,trans,all.x=T)
rm(trans)

eqtls_dt


#trans ens in symbol
trans_sy<-TransEnsemblVerstoSymbol(unique(eqtls_dt$gene_id))
trans_sy<-trans_sy[,gene:=hgnc_symbol][,gene_id:=ensembl_gene_id_version][,.(gene,gene_id)]
eqtls_dt<-merge(eqtls_dt,trans_sy,all.x=T,by="gene_id")

#add ensemble version for tissue wide gene
mart<-GetMartGenes()
attr<-GetBiomartAttrs(mart)
attr[str_detect(name,"start"),]


trans_vers<-data.table(getBM(attributes = c('ensembl_gene_id',"ensembl_gene_id_version"),
      filters = 'ensembl_gene_id', 
      values = unique(eqtls_dt$gene_id), 
      mart = mart))

eqtls_dt[tissue=="tissue_wide",ensembl_gene_id:=gene_id]
eqtls_dt<-merge(eqtls_dt,trans_vers[,.(ensembl_gene_id,ensembl_gene_id_version)],by="ensembl_gene_id",all.x=T)
eqtls_dt<-eqtls_dt[tissue=="tissue_wide",gene_id:=ensembl_gene_id_version][,-c("ensembl_gene_id","ensembl_gene_id_version")]
eqtls_dt<-eqtls_dt[!is.na(gene_id)]

#add tss_dist
tss_poss<-data.table(getBM(attributes = c('ensembl_gene_id_version', 'transcription_start_site',"strand"),
      filters = 'ensembl_gene_id_version', 
      values = unique(eqtls_dt$gene_id), 
      mart = mart))

tss_poss[,gene_id:=ensembl_gene_id_version][,tss_pos:=transcription_start_site][,strand:=ifelse(strand==1,"+","-")]
tss_poss<-unique(tss_poss,by="gene_id") #10k/30k gene_id n'ont pas de tss pos, ont les rm on mergeant
eqtls_dt<-merge(eqtls_dt,tss_poss[,.(gene_id,tss_pos,strand)])
eqtls_dt[is.na(strand)] #ok 
fwrite(eqtls_dt,"ref/eQTL/all_eQTLs_locally_signif.csv.gz",sep=";")
saveRDS(eqtls,"ref/eQTL/all_eQTLs_list.rds")
rm(eqtls,eqtls_dt)

#uniformize meth change by region
regions_meths<-Reduce(function(x,y)rbind(x,y,fill=T),lapply(c("liver","cbps","nash","placenta"),function(tissue)fread(paste0(out,tissue,"_regions_meth.change.csv.gz"))[,exp:=tissue][,genome_ref:=ifelse(exp=="placenta","hg38","hg19")]))
regions_meths[exp%in%c("liver","placenta"),tissue_ref:=exp]
regions_meths[exp=="cbps",tissue_ref:="whole_blood"]
regions_meths[exp=="nash",tissue_ref:="liver"]
regions_meths
regions_meths[,length.reg:=end-start]
unique(regions_meths$exp)
fwrite(regions_meths,"analyses/genescore_by_feature_region/all_tissue_regions_meth.change.csv.gz",sep=";")

#uniformize expr change 
expr_change<-Reduce(rbind,lapply(c("nash","liver","cbps","placenta"),function(expe)unique(fread(paste0("analyses/genescore_by_feature_region/",expe,"_regions_meth.change_degs.csv.gz")),by="gene")[,.(gene,log2FoldChange,pvalue,padj)][,exp:=expe]))
  # add gene_id
trans_sy<-TransSymboltoEnsemblVers(unique(expr_change$gene))
trans_sy[,gene_id:=ensembl_gene_id_version][,gene:=hgnc_symbol]
expr_change<-merge(expr_change,trans_sy[,.(gene,gene_id)],all.x=T)
  # filter
expr_change<-expr_change[!is.na(pvalue)]

fwrite(expr_change,"analyses/genescore_by_feature_region/all_tissue_regions_expr.change.csv.gz")

#Correl enhancer FeatureRegion with Fold change of genes eQTL linked ?
  #gene link to region if eQTL in the region
source("scripts/utils/new_utils.R")
out<-"analyses/genescore_by_feature_region/correl_meth.change_genes_expr_eqtl_links"
dir.create(out)

regions_meths<-fread("analyses/genescore_by_feature_region/all_tissue_regions_meth.change.csv.gz")
eQTLs<-fread("ref/eQTL/all_eQTLs_locally_signif.csv.gz")
unique(eQTLs[minLocal==T]$tissue)

i<-0
for(expe in unique(regions_meths$exp)){
  message(expe)
  regions_meths_tissue<-regions_meths[exp==expe][,.(cpg_id,region_id,chr,start,end,length.reg,chromatin_state,
                                                       chromatin_feature,meth.change,pval,n_cpg.reg,
                                                       tissue_ref,genome_ref,exp)]
  regions<-unique(regions_meths_tissue,by="region_id")
  
  #add Fold change of genes with eQTL in the region [TO FINISH]
  
  eQTL_tissue<-eQTLs[minLocal==T&tissue%in%c(unique(regions$tissue_ref),"tissue_wide")]
  
  
  if(unique(regions_meths_tissue$genome_ref)=="hg19"){
    eQTL_tissue[,start:=pos.hg19][,end:=pos.hg19+1][,chr:=chr.hg19]
    
  }else{
    eQTL_tissue[,start:=pos.hg38][,end:=pos.hg38+1][,chr:=chr.hg38]
    
  }
  eQTL_tissue<-eQTL_tissue[!is.na(start)&!is.na(chr)]

  message(nrow(eQTL_tissue),"eQTLs matching with chromatin regions")
    inter_reg_eqtl<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
            b=unique(eQTL_tissue[,.(chr,start,end,variant_id)])[order(chr,start)],
            select = c(4,8),col.names = c("region_id","variant_id"))
  
  message(nrow(inter_reg_eqtl),"match")
  
  
  regions_cpgs_meth_eqtls<-merge(regions_meths_tissue,inter_reg_eqtl,by="region_id",allow.cartesian=TRUE)
  
  message("number of eQTLs by chromatin_feature:")
  print(table(unique(regions_cpgs_meth_eqtls,by=c("region_id","variant_id"))$chromatin_feature))
  
  regions_cpgs_meth_eqtls[,n.eqtl.reg:=length(unique(variant_id)),by="region_id"]
  regions_cpgs_meth_eqtls[,avg.n.eqtl.reg:=mean(n.eqtl.reg),by="chromatin_feature"]
  regions_cpgs_meth_eqtls[,avg.n.eqtl.reg_norm:=avg.n.eqtl.reg/mean(length.reg),by="chromatin_feature"]                 
                     
  p<-ggplot(unique(regions_cpgs_meth_eqtls,by=c("region_id","variant_id")))+
    geom_col(aes(x=chromatin_feature,y=avg.n.eqtl.reg_norm))
  ggsave(fp(out,paste0(expe,"normalized_number_of_eqtl_by_chromatin_feature.png")),plot = p)
  
  regions_cpgs_meth_eqtls_genes<-merge(regions_cpgs_meth_eqtls,eQTL_tissue[,.(variant_id,gene_id,gene,tss_pos,strand)],by="variant_id",allow.cartesian=TRUE)
  
  regions_cpgs_meth_eqtls_genes[strand=="+",tss_dist:=round(mean(c(start[1],end[1])))-tss_pos,by=c("region_id","gene_id")]
  regions_cpgs_meth_eqtls_genes[strand=="-",tss_dist:=tss_pos-round(mean(c(start[1],end[1]))),by=c("region_id","gene_id")]
  i<-i+1
  if(i==1){
    regions_cpgs_meth_eqtls_genes_all<-regions_cpgs_meth_eqtls_genes
  }else{
    regions_cpgs_meth_eqtls_genes_all<-rbind(regions_cpgs_meth_eqtls_genes_all,regions_cpgs_meth_eqtls_genes)
  }
  
}

#harmonize a little
regions_cpgs_meth_eqtls_genes_all[,chromatin_feature:=ifelse(chromatin_feature[1]=="Inactive Enhancer","Poised Enhancer",chromatin_feature[1]),by="region_id"]
  #save
fwrite(regions_cpgs_meth_eqtls_genes_all,fp(out,"all_tissue_regions_meth.change.eqtls_genes_linked.csv.gz"),sep=";")

    #stats
      #number of eQTLs in promoter, enhancer..
table(unique(regions_cpgs_meth_eqtls_genes_all,by=c("region_id","variant_id"))[,.(chromatin_feature,exp)])
# chromatin_feature    cbps liver  nash placenta
#   Active Enhancer    5021  1413   184      215
#   Gene Body         10529  6117   112      332
#   Heterochromatin   14293 52163  1137        0
#   Inactive Enhancer     0  2371   349     1957
#   Poised Enhancer    1995     0     0        0
#   Promoter           3437  1525    26       45
#   state5                0     0     0    49297
#   state6                0     0     0       43

      #number of genes
table(unique(regions_cpgs_meth_eqtls_genes_all,by=c("region_id","gene_id"))[,.(chromatin_feature,exp)])
#                    exp
# chromatin_feature    cbps liver  nash placenta
#   Active Enhancer    4730  1318   183      224
#   Gene Body          9965  4169    88      271
#   Heterochromatin   10736 26468   602        0
#   Inactive Enhancer     0  2269   341     1971
#   Poised Enhancer    1965     0     0        0
#   Promoter           3354  1511    26       46
#   state5                0     0     0    30559
#   state6                0     0     0       48

      #number of eQTL by kb in promoter, enhancer..
regions_cpgs_meth_eqtls_genes_all[,n_eQTL:=length(unique(variant_id)),by=c("region_id","exp")]

regions_cpgs_meth_eqtls_genes_all[,n_eQTL.kb:=n_eQTL/length.reg*1000]


ggplot(unique(regions_cpgs_meth_eqtls_genes_all,by=c("region_id","exp")))+geom_boxplot(aes(x=chromatin_feature,y=n_eQTL.kb))+facet_wrap("exp")+scale_y_log10()


    #merge with degs
expr_change<-fread("analyses/genescore_by_feature_region/all_tissue_regions_expr.change.csv.gz")
regions_cpgs_meth_eqtls_genes_expr_all<-merge(regions_cpgs_meth_eqtls_genes_all[,-"gene"],expr_change,by=c("gene_id","exp"))

    #interpret
      #number of DMC in DEGs associated region in promoter, enhancer..
regions_cpgs_meth_eqtls_genes_expr_all[,degs_asso_reg:=any(padj<0.2),by=c("region_id","exp")]
regions_cpgs_meth_eqtls_genes_expr_all[,ratio.dmc:=sum(pval<0.01)/length(pval),by=c("region_id","exp")]


ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="cbps"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=ratio.dmc,fill=degs_asso_reg),outlier.shape = NA)+coord_cartesian(ylim=c(0,0.25))


ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="nash"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=ratio.dmc,fill=degs_asso_reg),outlier.shape = NA)+coord_cartesian(ylim=c(0,0.1))


ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="placenta"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=ratio.dmc,fill=degs_asso_reg),outlier.shape = NA)+coord_cartesian(ylim=c(0,0.1))


      #meth change sig avg correl with expr change ?
regions_cpgs_meth_eqtls_genes_expr_all[padj<0.2,deg:=ifelse(log2FoldChange>0,"upreg","downreg")]
regions_cpgs_meth_eqtls_genes_expr_all[padj>=0.2,deg:="no_change"]

regions_cpgs_meth_eqtls_genes_expr_all[pval<0.01,avg.meth.change_sig:=mean(meth.change),by=c("region_id","exp")]

regions_cpgs_meth_eqtls_genes_expr_all[,deg:=factor(deg,levels = c("no_change","downreg","upreg"))]

ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[!is.na(padj)&pval<0.01&exp=="cbps"],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change_sig,col=deg),outlier.shape = NA)
  
ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[!is.na(padj)&pval<0.01&exp=="nash"],by=c("region_id","gene_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change_sig,col=deg),outlier.shape = NA)

      # meth.change sd in promoter, enhancer.. [plus de sd in enhancer region ]
regions_cpgs_meth_eqtls_genes_expr_all[,meth.change.sd:=sd(meth.change),by=c("region_id","exp")]
    
ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[degs_asso_reg==T&exp=="cbps"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=meth.change.sd),outlier.shape = NA)
  
  # meth.change sd in DEGs associated region compared to non DEGs asso regions in promoter, enhancer.. [plus de sd diff btw degs/non-degs regions of enhaner feature than other region]
ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="cbps"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=meth.change.sd,fill=degs_asso_reg),outlier.shape = NA)
  

  # meth.change avg in DEGs associated region compared to non DEGs asso regions in promoter, enhancer.. [plus de meth.change in enhancer than non regulatory region ?]
regions_cpgs_meth_eqtls_genes_expr_all[,meth.change.avg:=mean(meth.change),by=c("region_id","exp")]

ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="cbps"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=meth.change.avg,fill=degs_asso_reg),outlier.shape = NA)
  
regions_cpgs_meth_eqtls_genes_expr_all[,mlog10pval.avg:=mean(-log10(pval)),by=c("region_id","exp")]
ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="cbps"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=mlog10pval.avg,fill=degs_asso_reg),outlier.shape = NA)
  


ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="liver"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=meth.change.avg,fill=degs_asso_reg),outlier.shape = NA)
  
ggplot(unique(regions_cpgs_meth_eqtls_genes_expr_all[exp=="liver"],by=c("region_id")))+
  geom_boxplot(aes(x=chromatin_feature,y=mlog10pval.avg,fill=degs_asso_reg),outlier.shape = NA)
  
fwrite(regions_cpgs_meth_eqtls_genes_expr_all,fp(out,"all_tissue_regions_meth.change.eqtls_genes_linked_expr.change.csv.gz"))

#test increase window of eQTL in cbps

#validate most signif closest to enhancer / regulatory element than less signif

