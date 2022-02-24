source("scripts/utils/new_utils.R")
out<-"outputs/cst_2e_annee"
dir.create(out)


#hto specific gene expression
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps6<-cbps[,str_detect(cbps$batch,"[0-6]")]

#major gene expr differences
library(ComplexHeatmap)
genes_diff<-data.table(FindMarkers(group.by = "hto",ident.1 = "TRUE",ident.2 ="FALSE" ),keep.rownames = "gene")
genes_diff #4078 DEGs
avg_expr<-AverageExpression(cbps6,group.by = "sample",assays = "SCT",slot = "data")$SCT

Heatmap(t(scale(t(avg_expr[genes_diff$gene,]))),
                           top_annotation = HeatmapAnnotation(df=data.frame(unique(data.table(cbps6@meta.data),by="sample"),row.names = "sample")[colnames(avg_expr),c('hto',"batch",'group','sex')] ),
                           show_row_names = F,name = "Expression")
                   

avg_expr_diff<-merge(data.table(avg_expr$SCT,keep.rownames = "gene"),data.table(genes_diff,keep.rownames = "gene"),by="gene")
ggplot(avg_expr_diff)+geom_point(aes(x=log1p(HTO),y=log1p(non_HTO),col=p_val_adj<0.001&abs(avg_log2FC)>0.5))
#stemness gatekeeper hto specific
genes_int<-c("SOCS3","HES1","JUNB","EGR1","KLF2","FOS")
VlnPlot(cbps6,genes_int,group.by="hto")

genes_exp<-merge(melt(data.table(as.matrix(cbps6@assays$SCT@data[genes_int,]),keep.rownames = "gene"),id.vars = "gene",variable.name = "bc",value.name = "expression"),
                 data.table(cbps6@meta.data,keep.rownames = "bc"))

ggplot(genes_exp)+geom_boxplot(aes(x=gene,y=expression,fill=hto))


ggplot(unique(genes_exp[,avg.expr:=mean(expression),by=.(gene,sample)],by=c("sample","gene")))+
  geom_boxplot(aes(x=hto,y=avg.expr,fill=hto))+facet_wrap("gene")

#ggplot(cbps6@meta.data)+geom_bar(aes(x=sample,fill=lineage_hmap),position = "fill")+facet_wrap("hto",scales = "free")

#altered lga stim go bp
res_go<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
res_go[p.adjust<0.1]$Description

res_kegg<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg_dn.csv")
res<-rbind(res_go[,source:="GO BP"],res_kegg[,source:="KEGG"])
res[,n.overlap:=as.numeric(Count)]
res[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]

res[,pct.overlap:=n.overlap/n.gene_set]

terms_of_interest<-c("regulation of cell growth",
                     "negative regulation of cell growth",
                     "cellular response to cytokine stimulus"   ,
                     "cellular response to chemical stimulus",
                     "cell activation" ,
                     "cellular response to external stimulus" ,
                     "response to stress"    ,
                     "regulation of cellular response to stress" ,
                     "regulation of response to stimulus" ,
                     "cellular response to stress",
                     "regulation of cellular response to stress",
                     "cellular response to growth factor stimulus",
                     "regulation of cellular response to growth factor stimulus" ,
                     "cellular response to extracellular stimulus" ,
                     "regulation of response to external stimulus" ,
                     "negative regulation of response to stimulus" ,
                     "positive regulation of response to cytokine stimulus" ,
                     "positive regulation of response to stimulus",
                     "regulation of stress-activated MAPK cascade"  ,
                     "regulation of cell activation",
                      "cellular response to starvation",
                     "negative regulation of cell activation" ,
                     "positive regulation of cell activation",
                      "regulation of cell differentiation",
                     "positive regulation of cell differentiation" ,
                     "negative regulation of cell differentiation",
                     "cell differentiation"  ,
                     "cellular response to chemical stress",
                     "response to temperature stimulus",
                     "regulation of cell differentiation" ,
                     "stem cell differentiation",
                     "hemopoiesis"  ,
                     "regulation of hemopoiesis" ,
                     "negative regulation of hemopoiesis","positive regulation of hemopoiesis",
                     "response to oxidative stress",
                     "erythrocyte differentiation","lymphocyte differentiation","myeloid cell differentiation",
                     "myeloid cell homeostasis","myeloid cell homeostasis","erythrocyte homeostasis","lymphocyte homeostasis",
                     "hematopoietic progenitor cell differentiation"   ,
                     "lymphoid progenitor cell differentiation" ,
                     "cell population proliferation",
                     "regulation of cell population proliferation",
                     "cell cycle"  ,
                     "regulation of cell cycle",
                     "MAPK signaling pathway",
                    "TNF signaling pathway",
                    "FoxO signaling pathway",
                    "NF-kappa B signaling pathway",
                    "Hippo signaling pathway",
                    "TGF−beta signaling pathway" )


ggplot(res,aes(x=pct.overlap,y=-log10(p.adjust),col=p.adjust<0.05))+geom_point()+
  geom_label_repel(aes(label=ifelse(p.adjust<0.05&Description%in%c(terms_of_interest),Description,"")),
                   max.overlaps = 3000)+
  facet_wrap("source")+
  scale_color_manual(values = c("grey",'red'))+theme_minimal()
ggsave(fp(out,"supp-points_hto_up_dn_kegg_go.pdf"))

go<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.rds")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
emapplot(pairwise_termsim(go,showCategory = 20),showCategory = 20)

#HSC, 
res_dt<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=='HSC']
#go bp up
res_go_bp_up<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1)
saveRDS(res_go_bp_up,fp(out,"res_hsc_hto_lga_go_bp_up.rds"))
res_go_bp_up_dt<-data.table(as.data.frame(res_go_bp_up))
fwrite(res_go_bp_up_dt,fp(out,"res_hsc_hto_lga_go_bp_up.csv"))
#go bp dn
res_go_bp_dn<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1)
saveRDS(res_go_bp_dn,fp(out,"res_hsc_hto_lga_go_bp_dn.rds"))
res_go_bp_dn_dt<-data.table(as.data.frame(res_go_bp_dn))
fwrite(res_go_bp_dn_dt,fp(out,"res_hsc_hto_lga_go_bp_dn.csv"))

emapplot(pairwise_termsim(res_go_bp_dn,showCategory = 20),showCategory = 20)

res_kegg_up<-enrichKEGG(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1)
res_kegg_up_dt<-data.table(as.data.frame(res_kegg_up))

fwrite(res_kegg_up_dt,fp(out,"res_hsc_hto_lga_kegg_up.csv"))


res_kegg_dn<-enrichKEGG(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1)
res_kegg_dn_dt<-data.table(as.data.frame(res_kegg_dn))

fwrite(res_kegg_dn_dt,fp(out,"res_hsc_hto_lga_kegg_dn.csv"))

res_go<-rbind(res_go_bp_up_dt[,query:="up"],res_go_bp_dn_dt[,query:="dn"])

res_kegg<-rbind(res_kegg_up_dt[,query:="up"],res_kegg_dn_dt[,query:="dn"])
res<-rbind(res_go[,source:="GO BP"],res_kegg[,source:="KEGG"])
res[,n.overlap:=as.numeric(Count)]
res[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]

res[,pct.overlap:=n.overlap/n.gene_set]


ggplot(res,aes(x=pct.overlap,y=-log10(p.adjust),col=p.adjust<0.05))+geom_point()+
  geom_label_repel(aes(label=ifelse(p.adjust<0.05&Description%in%c(terms_of_interest),Description,"")),
                   max.overlaps = 3000)+
  facet_grid(source~query)+
  scale_color_manual(values = c("grey",'red'))+theme_minimal()

#SCENIC 
library(Seurat)
cbps<-readRDS("outputs/10-SCENIC/cbps_with_regulons_activity.rds")
Idents(cbps)<-"lineage_hmap"
FeaturePlot(cbps,c("EGR1","KLF2","JUNB"),
            label=T,min.cutoff = 0.2,
            max.cutoff = "q95",ncol = 3)
