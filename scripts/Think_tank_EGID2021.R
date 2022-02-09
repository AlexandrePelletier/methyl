source("scripts/utils/new_utils.R")

#genes of stemness / growth 
res_k_meth<-fread("outputs/03-pathway_analysis/res_gsea_kegg.csv")

stem_id<-res_k_meth[str_detect(Description,"Wnt|Hippo|Fox|pluri")]$ID

growth_id<-res_k_meth[str_detect(Description,"PI3K|Insulin sig|Cell cycle|Growth horm")]$ID

getGenesKEGGPathw<-function(pathID){
  require(KEGGREST)
  g<-keggGet(pathID)[[1]]$GENE
  g<-g[1:length(g)%%2==0]
  return(as.vector(sapply(g,function(x)strsplit(x,";")[[1]][1])))
}

stem_g<-getGenesKEGGPathw("")
growth_g<-getGenesKEGGPathw("")

res_h<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_all_cbps.csv.gz")
ggplot(res_h)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.1&abs(log2FoldChange)>0.5))+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = "bottom")

res_h[padj<0.1&abs(log2FoldChange)>0.5]
res_h[padj<0.1&log2FoldChange<(-0.5)] #122
res_h[padj<0.1&log2FoldChange>0.5]

#pathway
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
res_kegg<-enrichKEGG(bitr(res_h[padj<0.1&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1)
res_kegg_dt<-data.table(as.data.frame(res_kegg))
res_kegg_dt[p.adjust<0.05]

emapplot(pairwise_termsim(res_kegg,showCategory = 15),showCategory = 15)

#proportion HSC/MPP  activated
mtd<-data.table(readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/all_cbps_with_regulons_activity.rds")@meta.data)
mtd$lineage[mtd$lineage=="unknown-1"]<-"MPP"
mtd[,n.sample:=.N,.(sample,orig.ident)]
mtd[,pct.lin:=.N/n.sample,.(sample,orig.ident,lineage)]

mtdsl<-unique(mtd[ambigous==F&group!="iugr"],by=c("sample","orig.ident","lineage"))
ggplot(mtdsl[lineage%in%c('HSC',"MPP")&hto==T])+
  geom_boxplot(aes(x=group,y=pct.lin,fill=group))+facet_wrap("lineage",scales = "free")+
  theme_bw()



####supp####
res_b<-fread("outputs/07-LGA_vs_Ctrl_Basal/res_pseudobulkDESeq2_all_cbps.csv.gz")[,stimulation:=F]
res_h<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_all_cbps.csv.gz")[,stimulation:=T]


res<-rbind(res_b,res_h)

ggplot(res)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.1&abs(log2FoldChange)>0.5))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()
                           
paste(res[padj<0.1&log2FoldChange>0.5]$gene,collapse = ", ")

res_b_hsc<-fread("outputs/07-LGA_vs_Ctrl_Basal/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,stimulation:=F]

res_h_hsc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,stimulation:=T]

res_hsc<-rbind(res_b_hsc,res_h_hsc)

ggplot(res_hsc)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()
                           

ggplot(res_hsc)+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.6))+
  facet_wrap("stimulation")+scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = 'bottom')
                           

#TF Acutivity GATA1 SPI1 CEBPA
library(Seurat)
cbps<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/all_cbps_with_regulons_activity.rds")

DefaultAssay(cbps)<-"TF_AUC"
Idents(cbps)<-"lineage"
Idents(cbps)[Idents(cbps)=="unknown-1"]<-"MPP"
FeaturePlot(cbps,c("GATA1"),label=T,cols = c("grey","red"),max.cutoff = 0.3,min.cutoff = 0.1)
