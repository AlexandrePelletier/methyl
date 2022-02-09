source("scripts/utils/new_utils.R")
out<-"analyses/genescore_by_feature_region/correl_meth.change_genes_expr"
regions_meths<-Reduce(function(x,y)rbind(x,y,fill=T),lapply(c("liver","cbps","nash","placenta"),function(tissue)fread(paste0(out,tissue,"_regions_meth.change.csv.gz"))[,tissue:=tissue][,genome_ref:=ifelse(tissue=="placenta","hg38","hg19")]))

for(tiss in unique(regions_meths$tissue)){
  regions_meths_tissue<-regions_meths[tissue==tiss][,.(cpg_id,region_id,chr,start,end,chromatin_state,
                                                       chromatin_feature,meth.change,pval,n_cpg.reg,
                                                       tissue,genome_ref)]
  regions<-unique(regions_meths_tissue,by="region_id")
  #add Fold change of genes +/-100kb
  if(unique(regions_meths_tissue$genome_ref)=="hg19"){
    tss_genes<-fread("ref/hg19/tss_pos_hg19_refseq_curated_250121.txt.gz")
  }else{
    tss_genes<-fread("ref/hg38/tss_pos_hg38_refseq_curated_250121.txt.gz")
  }

  tss_genes[,start:=tss_pos-100000][start<0,start:=0][,end:=tss_pos+100000]
  tss_genes[order(chr,start)]
  inter_reg_genes<-bed_inter(a=regions[,.(chr,start,end,region_id)][order(chr,start)],
            b=tss_genes[,.(chr,start,end,gene_id)][order(chr,start)],
            select = c(4,8),col.names = c("region_id","gene_id"))
  regions_cpgs_meth_genes<-merge(regions_meths_tissue,inter_reg_genes,by="region_id",allow.cartesian=TRUE)
  
  expr_change<-unique(fread(paste0("analyses/genescore_by_feature_region/",tiss,"_regions_meth.change_degs.csv.gz")),by="gene")
  expr_change<-merge(expr_change,unique(tss_genes[,.(gene_id,gene,tss_pos,strand)]))[,.(gene_id,gene,tss_pos,strand,log2FoldChange,pvalue,padj)]
  expr_change[is.na(padj),padj:=1]
  regions_cpgs_meth_genes_expr<-merge(regions_cpgs_meth_genes,expr_change,by="gene_id",all.x=T)
  regions_cpgs_meth_genes_expr[strand=="+",tss_dist:=round(mean(c(start[1],end[1])))-tss_pos,by=c("region_id","gene_id")]
  regions_cpgs_meth_genes_expr[strand=="-",tss_dist:=tss_pos-round(mean(c(start[1],end[1]))),by=c("region_id","gene_id")]
  
  regions_cpgs_meth_genes_expr[pval<0.001,avg.meth.change.cpg_sig:=mean(meth.change),by="region_id"]
  regions_cpgs_meth_genes_expr[pval<0.001,med.meth.change.cpg_sig:=median(meth.change),by="region_id"]
  
  regions_cpgs_meth_genes_expr[padj<0.1,deg:=ifelse(log2FoldChange>0,"upreg","downreg")]
  regions_cpgs_meth_genes_expr[padj>=0.1,deg:="no_change"]
  regions_cpgs_meth_genes_expr[,deg:=factor(deg,levels = c("no_change","downreg","upreg"))]
  
  #all regions_genes link :
  p1<-ggplot(unique(regions_cpgs_meth_genes_expr[!is.na(padj)&pval<0.001],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change.cpg_sig,col=deg),outlier.shape = NA)
  
  p2<-ggplot(unique(regions_cpgs_meth_genes_expr[!is.na(padj)&pval<0.001],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=med.meth.change.cpg_sig,col=deg),outlier.shape = NA)
  
  p_all<-p1+p2
  ggsave(fp(out,paste0(tiss,"_all_degs_100kb_correl_with_feature_region_meth.change.png")),plot = p_all,width = 12)
  
  #for each regions, best genes in term of degs :
  regions_cpgs_meth_genes_expr[!is.na(padj),best_deg.reg:=gene==gene[which.min(padj)],by="region_id"]
  
  p1<-ggplot(unique(regions_cpgs_meth_genes_expr[!is.na(padj)&pval<0.001&best_deg.reg==T],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=tss_dist,fill=tss_dist>0),outlier.shape = NA)+scale_y_log10()
  
  p2<-ggplot(unique(regions_cpgs_meth_genes_expr[!is.na(padj)&pval<0.001&best_deg.reg==T],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change.cpg_sig,fill=deg),outlier.shape = NA)
  
  p3<-ggplot(unique(regions_cpgs_meth_genes_expr[!is.na(padj)&pval<0.001&best_deg.reg==T],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=med.meth.change.cpg_sig,fill=deg),outlier.shape = NA)
  
  p_all<-p1+p2+p3
  ggsave(fp(out,paste0(tiss,"_best_deg_100kb_correl_with_feature_region_meth.change.png")),plot = p_all,width = 12)
  
  
  fwrite(regions_cpgs_meth_genes_expr,fp(out,paste0(tiss,"_all_regions_meth.change_and_degs_100kb.csv.gz")))
  
}
