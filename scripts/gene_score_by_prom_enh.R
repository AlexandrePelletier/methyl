#genescore by prom/ enh, best correl with Foldchange of DEGs ?
source("scripts/utils/new_utils.R")
out<-"outputs/gene_score_by_prom_enh"
dir.create(out)
res_cl<-fread("outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")


#gene_score calculation
res_cl[,region_type:=ifelse(abs(tss_dist)<=2000,"promoter","other"),by="gene"]
res_cl[,n_cpg_weight:=(1/sum(1/(abs(CpGScore)+1)))^(1/4),by=c('region_type',"gene")]
res_cl[,gene_score:=sum(CpGScore)*n_cpg_weight,by=c('region_type',"gene")]

#classical approach calculation
res_cl[,avg.meth.change:=mean(meth.change),by=c('region_type',"gene")]
res_cl[,avg.m.log10.pval:=mean(-log10(pval)),by=c('region_type',"gene")]
res_cl[,avg.dmc_score:=mean(-log10(pval)*meth.change),by=c('region_type',"gene")]


#find DMR
dmr_res<-fread("outputs/model13_all_compas/DMR_with_comb_p/C.L.regions-p.bed.gz")
#merge with gene
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("org.Hs.eg.db")
dmr_res<-dmr_res[,chr:=`#chrom`][,-1]
dmr_res<-dmr_res[,.(chr,start,end,min_p,n_probes,z_p,z_sidak_p)]

dmr<-makeGRangesFromDataFrame(dmr_res[,1:3])

anno_dmr <- annotatePeak(dmr, tssRegion=c(-2000, 2000),
                       TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                       annoDb="org.Hs.eg.db")
anno_dmr<-as.data.frame(anno_dmr)
head(anno_dmr,10)
fwrite(anno_dmr,fp(out,"chipseaker_anno_dmr.csv.gz"),sep=";")
anno_dmr<-fread(fp(out,"chipseaker_anno_dmr.csv.gz"),sep=";")

dmr_tss<-data.table(anno_dmr)[,chr:=seqnames][,tss_dist:=distanceToTSS][,gene:=SYMBOL][,.(chr,start,end,tss_dist,gene)]

dmr_res<-merge(dmr_res,dmr_tss,all.x=T)
dmr_res[z_sidak_p<0.05&abs(tss_dist)<2000]$gene
dmr_res[gene=="SOCS3"]
res_cl[gene=="SOCS3"&abs(tss_dist)<2000]

#r2 avec DEGs (top1000)

res_de_cl<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")
res_de_cl[,top1000:=pvalue<=sort(pvalue)[1000]]

res_cl_1000<-merge(unique(res_cl[order(pval)],by=c("gene","region_type")),res_de_cl[top1000==T],by=c("gene"))
res_cl_1000[,r2:=cor(abs(gene_score),abs(log2FoldChange))^2,by=c("region_type")]

res_cl_1000



