### Project Setup ==================================================================================
library(here)
out <- here("outputs", "03-pathway_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  source("scripts/utils/new_utils.R")
  library(clusterProfiler)
  library(parallel)
})


### Tables and Figures Theme =======================================================================
# theme_set(theme_light())


### Functions ======================================================================================


### Analysis =======================================================================================

gwas_genes_rep10<-fread("ref/gwas/reported_gene_traits_GWAS_10genes.csv")
res_perm<-fread("outputs/02-gene_score_calculation_and_validation/res_1000perm_genescore_add.csv.gz")

res_gw<-fread("outputs/03-pathway_analysis/res_gsea_gwas.csv")
traits_sig<-res_gw[p.adjust<0.001]$ID #already filter but for understanding
gwas_genes_sig<-gwas_genes_rep10[disease_trait%in% traits_sig]

res_gw_perm<-Reduce(rbind,mclapply(1:1000,function(i){
  resg<-res_perm[perm==i]
  gene_scores<-resg$gene_score_add
  names(gene_scores)<-resg$gene
  gene_scores<-sort(gene_scores,decreasing = T)

  res_gsea_gwas<- data.table(as.data.frame(GSEA(geneList = rank(gene_scores),
                                              TERM2GENE = gwas_genes_sig[,.(disease_trait,reported_gene)],
                                              maxGSSize    = 500,
                                              eps = 0,
                                              pvalueCutoff = 1
                                              )))
  
  return(res_gsea_gwas[,perm:=i][,.(ID,p.adjust,perm)])},mc.cores = 20,mc.silent=T))


res_gwp<-merge(res_gw,res_gw_perm,all=T)
res_gwp[,p.perm:=sum(p.adjust[is.na(perm)]>=p.adjust[!is.na(perm)])/1000,by="ID"]

message(nrow(res_gwp[is.na(perm)&p.perm<0.01]), " traits are signif !") 

fwrite(res_gwp[is.na(perm)][,-"perm"],"outputs/03-pathway_analysis/res_gsea_gwas_perm.csv")

### Complete =======================================================================================
message("Success!", appendLF = TRUE)
