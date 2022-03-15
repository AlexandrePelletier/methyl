out<-"outputs/19-gene_expression_analysis"
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

#top express genes - gsea
cbpsh<-subset(readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds"),hto==T)


genes_lin<-AverageExpression(cbpsh,assays = "SCT",group.by = "lineage_hmap",slot = "data")$SCT

genes_lin[order(genes_lin[,"HSC"],decreasing = T),"HSC"]

trad_genes<-data.table(bitr(rownames(genes_lin),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F))

res_gsea<-Reduce(rbind,lapply(colnames(genes_lin), function(lin){
  gl<-genes_lin[order(genes_lin[,lin],decreasing = T),lin]
  res_gs<-gseGO(gl,
                ont = "BP",
                org.Hs.eg.db,
                keyType = "SYMBOL",
                exponent = 1,
                minGSSize = 10,
                maxGSSize = 500,
                eps = 1e-50,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose = TRUE,
                seed = FALSE,
                by = "fgsea"
              )
  
  res_gs<-data.table(as.data.frame(res_gs))
  return(res_gs[,lineage:=lin])
}))
fwrite(res_gsea,fp(out,"res_gsea_top_express_genes_by_lineage_go_bp_padj0.05.csv.gz"))

res_gsea[,top30goterm:=rank(p.adjust)<=30,by="lineage"]
split(res_gsea[top30goterm==T,Description,by="lineage"],by="lineage") #pas de réponse au stress dans les tops go enriched

#regulons
cbpsh<-subset(readRDS("outputs/10-SCENIC/cbps_with_regulons_activity.rds"),hto==T)

regulon_lin<-AverageExpression(cbpsh,assays = "TF_AUC",group.by = "lineage_hmap",slot = "data")$TF_AUC
regulon_lin_dt<-melt(data.table(regulon_lin,keep.rownames = "regulon"),variable.name = "lineage",value.name = "activity",id.vars = "regulon")
regulon_lin_dt<-regulon_lin_dt[order(lineage,-activity)]
regulon_lin_dtf<-regulon_lin_dt[!str_detect(regulon,"e$")]
regulon_lin_dtf[,top10regulon:=rank(-activity)<=10,by="lineage"]

regulon_lin_dtf[,top30regulon:=rank(-activity)<=30,by="lineage"]
split(regulon_lin_dtf[top30regulon==T],by="lineage") #pas de réponse au stress dans les tops go enriched

fwrite(regulon_lin_dtf[top30regulon==T],fp(out,"top30_regulons_activity_by_lineage.csv"))

lin_of_int<-c("LT-HSC","HSC","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas")
top_regulon_lin<-regulon_lin[unique(regulon_lin_dtf[lineage%in%lin_of_int&top30regulon==T]$regulon),]

library(ComplexHeatmap)
pdf(fp(out,"heatmap_top30_regulon_by_lineage_hto_only.pdf"),width = 10,height = 6)
Heatmap(t(top_regulon_lin[,lin_of_int]),col = c("white","orange",'red',"darkred"),column_names_gp = gpar(fontsize = 8))
dev.off()
