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
    egions_cpgs_meth_eqtls_genes_all<-regions_cpgs_meth_eqtls_genes
  }else{
    regions_cpgs_meth_eqtls_genes_all<-rbind(regions_cpgs_meth_eqtls_genes_all,regions_cpgs_meth_eqtls_genes)
  }
  
  expr_change<-fread("analyses/genescore_by_feature_region/all_tissue_regions_expr.change.csv.gz")[exp==expe]
  expr_change<-expr_change[,.(gene,log2FoldChange,pvalue,padj)]
  expr_change[is.na(padj),padj:=1]
  regions_cpgs_meth_eqtls_genes_expr<-merge(regions_cpgs_meth_eqtls_genes,expr_change,by="gene",all.x=T)
  
  regions_cpgs_meth_eqtls_genes_expr[pval<0.001,avg.meth.change.cpg_sig:=mean(meth.change),by="region_id"]
  regions_cpgs_meth_eqtls_genes_expr[pval<0.001,med.meth.change.cpg_sig:=median(meth.change),by="region_id"]
  
  regions_cpgs_meth_eqtls_genes_expr[padj<0.1,deg:=ifelse(log2FoldChange>0,"upreg","downreg")]
  regions_cpgs_meth_eqtls_genes_expr[padj>=0.1,deg:="no_change"]
  regions_cpgs_meth_eqtls_genes_expr[,deg:=factor(deg,levels = c("no_change","downreg","upreg"))]
  
  #all regions_genes link :
  p1<-ggplot(unique(regions_cpgs_meth_eqtls_genes_expr[!is.na(padj)&pval<0.001],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change.cpg_sig,col=deg),outlier.shape = NA)
  
  p2<-ggplot(unique(regions_cpgs_meth_eqtls_genes_expr[!is.na(padj)&pval<0.001],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=med.meth.change.cpg_sig,col=deg),outlier.shape = NA)
  
  p_all<-p1+p2
  ggsave(fp(out,paste0(expe,"_all_degs_eQTL_correl_with_feature_region_meth.change.png")),plot = p_all,width = 12)
  
  #for each regions, best genes in term of degs :
  regions_cpgs_meth_eqtls_genes_expr[!is.na(padj),best_deg.reg:=gene==gene[which.min(padj)],by="region_id"]
  
  p1<-ggplot(unique(regions_cpgs_meth_eqtls_genes_expr[!is.na(padj)&pval<0.001&best_deg.reg==T],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=tss_dist,fill=tss_dist>0),outlier.shape = NA)+scale_y_log10()
  
  p2<-ggplot(unique(regions_cpgs_meth_eqtls_genes_expr[!is.na(padj)&pval<0.001&best_deg.reg==T],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=avg.meth.change.cpg_sig,fill=deg),outlier.shape = NA)
  
  p3<-ggplot(unique(regions_cpgs_meth_eqtls_genes_expr[!is.na(padj)&pval<0.001&best_deg.reg==T],by=c("region_id","gene_id")))+
    geom_boxplot(aes(x=chromatin_feature,y=med.meth.change.cpg_sig,fill=deg),outlier.shape = NA)
  
  p_all<-p1+p2+p3
  ggsave(fp(out,paste0(expe,"_best_deg_eQTL_correl_with_feature_region_meth.change.png")),plot = p_all,width = 12)
  
  
  fwrite(regions_cpgs_meth_eqtls_genes_expr,fp(out,paste0(expe,"_all_regions_meth.change_and_degs_eQTL_link.csv.gz")))
  
}
