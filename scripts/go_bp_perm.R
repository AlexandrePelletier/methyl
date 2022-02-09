### Project Setup ==================================================================================
library(here)
out <- here("outputs", "03-pathway_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  source("scripts/utils/new_utils.R")
  library(clusterProfiler)
  library(parallel)
  library(org.Hs.eg.db)

})


### Tables and Figures Theme =======================================================================
# theme_set(theme_light())


### Functions ======================================================================================


### Analysis =======================================================================================

res_g<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")
go_sig<-res_g[p.adjust<0.001]$ID #already filter but for understanding
res_perm<-fread("outputs/02-gene_score_calculation_and_validation/res_1000perm_genescore_add.csv.gz")

genes.df<-bitr(res_perm[perm==1]$gene,
                 fromType = 'SYMBOL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)

res_g_perm<-Reduce(function(x,y)rbind(x,y,fill=T),mclapply(1:1000,function(i){
  print(i)
  resg<-res_perm[perm==i]
  gene_scores<-resg$gene_score_add
  names(gene_scores)<-resg$gene
  
  gene_scores<-gene_scores[genes.df$SYMBOL]
  names(gene_scores)<-genes.df$ENTREZID
  gene_scores<-sort(gene_scores,decreasing = T)
  
  resgo<-gseGO(geneList     = rank(gene_scores), 
                                               ont = "BP",
                                                minGSSize    = 50,
                                                pvalueCutoff = 1,
                                                OrgDb = org.Hs.eg.db,
                                               verbose=F)
  if("gseaResult" %in% class(resgo)){
      res_gsea_go<- data.table(as.data.frame(resgo))
      return(res_gsea_go[,perm:=i][ID%in%go_sig][,.(ID,p.adjust,perm)])
  }else{
    return(data.table())
        }
  
  },mc.cores=10))
  

res_g[,perm:=0]
res_gp<-merge(res_g,res_g_perm,all=T)
res_gp[,p.perm:=sum(p.adjust[perm==0]>=p.adjust[perm>0],na.rm=T)/(sum(perm>0,na.rm=T),by="ID"]

message(nrow(res_gp[is.na(perm)&p.perm<0.01]), " traits are signif !") 

fwrite(res_gp[is.na(perm)][,-"perm"],"outputs/03-pathway_analysis/res_gsea_go_perm.csv")


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
