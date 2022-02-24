library(data.table)
source("scripts/utils/new_utils.R")
out<-"outputs/figures_epi_response2"
dir.create(out)

#% of the total number of analyzed CpGs fall in putative regulatory regions
cpgs_annot<-fread("outputs/02-gene_score_calculation_and_validation/cpgs_genes_annot_and_weight.csv.gz")
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
cpgs_annotf<-cpgs_annot[cpg_id%in%res_meth$cpg_id]
unique(cpgs_annotf[ensembl_reg_score>0|chromatin_feature%in%c(4,5,6)],by="cpg_id")
381075/756470#50%
cpgs_annotf[in_eQTR==T]
253927/756470

res_methg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

#2D- validation GS based on DEGs ####
res_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
resg<-unique(res_anno[order(gene,pval)],by="gene")
resg[,mlog10pval:=-log10(pval)]
meth_scores<-melt(resg[,.SD,.SDcols=c("gene","gene_score_add","min.pval","avg.pval","mlog10pval","meth.change","avg.mlog10.pval","avg.meth.change","max.dmc_score", "avg.dmc_score")],id.vars = "gene",variable.name = "meth_metric",value.name = "score")
meth_scores[score==Inf,score:=NA]
meth_scores[score==-Inf,score:=NA]
meth_scores<-meth_scores[!is.na(score)]
res_de_cl<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_all_cbps.csv.gz")

res_de_cl<-res_de_cl[!is.na(padj)]
meth_scores_de<-merge(meth_scores,res_de_cl,by=c("gene"))
unique(meth_scores_de[padj<0.05]$gene)#362 DEGS


meth_scores_de[,score_scaled:=scale(score),by="meth_metric"]
meth_scores_de$meth_metric<-factor(meth_scores_de$meth_metric,levels = c("min.pval","mlog10pval","meth.change","avg.pval","avg.mlog10.pval","avg.meth.change","max.dmc_score", "avg.dmc_score","gene_score_add"))
ggplot(meth_scores_de)+
  geom_boxplot(aes(fill=padj<0.05,y=score_scaled,x=meth_metric),outlier.shape = NA)+
  coord_cartesian(ylim = c(-3,3))+
  theme_classic()+
  scale_fill_manual(values = c("white","grey"))
ggsave(fp(out,"2D-boxplot_valid_genescore_with_degs_association.pdf"))

meth_scores_de[,pval:=wilcox.test(score[padj<0.05],score[padj>=0.05])$p.value,by="meth_metric"]
meth_scores_de[,score_diff:=mean(score_scaled[padj<0.05])-mean(score_scaled[padj>=0.05]),by="meth_metric"]
meth_scores_de[,FC:=mean(score[padj<0.05])/mean(score[padj>=0.05]),by="meth_metric"]

unique(meth_scores_de[,.(meth_metric,pval,score_diff,FC)],by="meth_metric")

#        meth_metric         pval
# 1:  gene_score_add 0.0006818586
# 2:        min.pval 0.0249488972
# 3:        avg.pval 0.7968159120
# 4:      mlog10pval 0.0757731002
# 5:     meth.change 0.0180098893
# 6: avg.mlog10.pval 0.3534897816
# 7: avg.meth.change 0.0368265523
# 8:   max.dmc_score 0.0192803029
# 9:   avg.dmc_score 0.4128018069

ggplot(unique(meth_scores_de[!meth_metric%in%c('min.pval','avg.pval')],by="meth_metric"))+
  geom_col(aes(y=-log10(pval),x=meth_metric,fill=FC))+
  scale_fill_gradient(low = "grey",high = "red")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(fp(out,"2D-barplot_valid_genescore_with_degs_association.pdf"))

#3A-hematomap ####
hmap<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/hematomap_ctrls_sans_stress/hematomap_ctrls_sans_stress.rds")
hmap<-subset(hmap,cell_type!="18")

#reorder lineage levels
Idents(hmap)<-"lineage"

levels(hmap)<-c("LT-HSC",
                "HSC",
                "MPP/LMPP",
                "Lymphoid",
                "B cell",
                "T cell",
                "Erythro-Mas",
                "Mk/Er",
                "Myeloid",
                "DC")
#3A-UMAP 
hmap[["lineage"]]<-Idents(hmap)
DimPlot(hmap,group.by = c("lineage"),label = T)
ggsave(fp(out,"3A-hematomap.pdf"))


#+% OF EACH SUBPOPULATION
round(table(Idents(hmap))/ncol(hmap)*100)

#3B-markers ####
#dotplot
markers<-fread("../singlecell/ref/hematopo_markers.csv")
markers[cell_type=="HSC"]
DefaultAssay(hmap)<-"SCT"
key_genes_lin<-c("ID1","DUSP2", #LT-HSC
           "AVP","FOS", #HSC
           "MLLT3","CDK6", #MPP
           "CD99", #LMPP
           "LTB", #CLP
           "IGHM","CD37", #B cell
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "PLEK","HBD", #Mk/Er
           "MPO","CEBPA","CTSG","AZU1", #GMP
         "CST3","CD83") #DC
DotPlot(hmap,features = key_genes_lin,
        group.by = "lineage",cols = c("white","black"))
ggsave(fp(out,"3B-dotplot_markers.pdf"))

#featureplot
ps<-FeaturePlot(hmap,
                features = c("ID1","AVP","CDK6","LTB","GATA1","MPO"),
                max.cutoff = "q99",
                combine = F)
ps[[3]]<-FeaturePlot(hmap,
                features = c("CDK6"),
                min.cutoff = "q60")
ps<-lapply(ps, function(x)x+NoAxes()+NoLegend())
wrap_plots(ps)
ggsave(fp(out,"3B-featureplot_markers.pdf"))


#3C-distr LGA / Ctrl
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_light.rds")
cbps<-subset(cbps,lineage_hmap!="18")

#dimplot
DimPlot(cbps,group.by = c("group"))
ggsave(fp(out,"3C-dimplot_group.pdf"))

#barplot
ggplot(cbps@meta.data)+geom_bar(aes(x=group,fill=lineage_hmap),position = "fill")
ggsave(fp(out,"3C-barplot_group.pdf"))

#3D- HTO effect ####
library(ComplexHeatmap)
library(Seurat)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_sct_light.rds")
cbps6<-cbps[,str_detect(cbps$batch,"[0-6]")]
#major gene expr differences
m_hto<-data.table(FindMarkers(cbps6,group.by = "hto",ident.1 = "TRUE",ident.2 ="FALSE" ),keep.rownames = "gene")
m_hto #4079 DEGs
fwrite(m_hto,fp(out,"3D-hto_markers.csv"))
avg_expr<-AverageExpression(cbps6,group.by = "sample",assays = "SCT",slot = "data")$SCT
pdf(fp(out,"3D-heatmap_HTO_effect.pdf"),width = 6,height = 6)

Heatmap(t(scale(t(avg_expr[m_hto$gene,]))),
                           top_annotation = HeatmapAnnotation(df=data.frame(unique(data.table(cbps6@meta.data),by="sample"),row.names = "sample")[colnames(avg_expr),c('hto',"batch",'group','sex')] ),
                           show_row_names = F,name = "Expression")
                   
dev.off()

#all cbps
avg_expr<-AverageExpression(cbps,group.by = "sample_hto",assays = "SCT",slot = "data")$SCT

Heatmap(t(scale(t(avg_expr[m_hto$gene,]))),
                           top_annotation = HeatmapAnnotation(df=data.frame(unique(data.table(cbps@meta.data),by="sample_hto"),row.names = "sample_hto")[colnames(avg_expr),c('hto',"batch",'group','sex')] ),
                           show_row_names = F,name = "Expression")

#cbp8 not hto response
#compared to Cbp7b?
m<-FindMarkers(cbps,group.by = "batch",ident.1 = "cbp8",ident.2 = "cbp7b")
m<-data.table(m,keep.rownames = "gene")
m[avg_log2FC>0]#IER, KLF,NFKBIA
head(m[avg_log2FC<0],100)
m[gene%in%c("SOCS3","HES1","EGR1","STAT3")]

#3D-MA plot HTO####
res_hto_dup<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")
fwrite(res_hto_dup[padj<0.05&abs(log2FoldChange)>0.5],fp(out,"3D-res_hto_signature_padj0.05_log2FC0.5.csv"))
genes_of_interest<-c("SESN2","SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

ggplot(res_hto_dup,aes(y=log2FoldChange,x=baseMean,col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +scale_x_log10()+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"3D-ma_plot_hto.pdf"))
#3E-pathway HTO ####
#gsea go bp
library(clusterProfiler)
library(enrichplot)
res_gsea_go<-readRDS("outputs/08-HTO_signature/res_gsea_go_bp.rds")

res_gsea_go_up<-res_gsea_go
res_gsea_go_up@result<-subset(res_gsea_go@result,NES>0)
pdf(fp(out,"3E-emapplot_gsea_go_bp_up_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_gsea_go_up),
         showCategory = sum(res_gsea_go_up@result$p.adjust<0.05),
         cex_label_category=0.66)
dev.off()


res_gsea_go_dn<-res_gsea_go
res_gsea_go_dn@result<-subset(res_gsea_go@result,NES<0)
pdf(fp(out,"3E-emapplot_gsea_go_bp_dn_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_gsea_go_dn),showCategory = 50,cex_label_category=0.66)
dev.off()



res_gsea_go_dt<-fread("outputs/08-HTO_signature/res_gsea_go_bp.csv.gz")
table(res_gsea_go_dt[p.adjust<0.05]$sens)
res_gsea_go_dt[p.adjust<0.05&sens=="up"]
ggplot(res_gsea_go_dt[p.adjust<0.05],aes(x=Description,y=-log10(p.adjust)))+
  geom_point(aes(size=n.enriched,col=NES))+
  facet_grid(~sens,scales = "free_x",space = "free_x")+
  scale_x_discrete(guide = guide_axis(angle=80))+
  scale_y_continuous()+
  scale_color_gradient2(low = "blue",mid="grey",high = "red")
#simplify with rrvgo
library(rrvgo)
simMatrix <- calculateSimMatrix(res_gsea_go_dt[p.adjust<0.05]$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(res_gsea_go_dt[p.adjust<0.05]$p.adjust), res_gsea_go_dt[p.adjust<0.05]$ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
redgo<-data.table(reducedTerms)

redgo[,top.term.cluster:=go==parent]
sum(redgo$top.term.cluster) #19
redgo[,ID:=go]
res_red<-merge(res_gsea_go_dt,redgo,by="ID")
res_red[,score.max:=max(score),by="cluster"]
res_red[top.term.cluster==T,top10clus:=rank(-score)<=15]

res_red[,term:=factor(term,levels=res_red[order(p.adjust)]$term)]
ggplot(res_red[top.term.cluster==T&p.adjust<0.05&sens=="up"],
       aes(x=term,y=-log10(p.adjust)))+
  geom_point(aes(size=n.enriched,col=NES))+
  scale_x_discrete(guide = guide_axis(angle=80))+
  scale_color_gradient(low = "white",high = "black")

ggsave(fp(out,"3E-dotplot_gsea_HTO_GO_BP_non_redundant_terms_padj0.05_up_only.pdf"),height=8)
fwrite(res_red[top.term.cluster==T&top10clus==T&p.adjust<0.05],fp(out,"3E-top_terms_clusters.csv"))
  

# enrichgo

res_go_up<-fread("outputs/08-HTO_signature/res_hto_signature_go_bp_up.csv")
res_go_dn<-fread("outputs/08-HTO_signature/res_hto_signature_go_bp_dn.csv")
res_go<-rbind(res_go_up[,query:="up"],res_go_dn[,query:="dn"])
res_go[,n.overlap:=as.numeric(Count)]
res_go[,n.query:=as.numeric(str_extract(GeneRatio,"[0-9]+$"))]
res_go[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]
res_go[,n.background:=as.numeric(str_extract(BgRatio,"[0-9]+$"))]

res_go[,pct.overlap.query:=n.overlap/n.gene_set]
res_go[,pct.overlap.background:=n.gene_set/n.background]
res_go[,fold.enrichment:=pct.overlap.query/pct.overlap.background]
res_go[,log2FE:=log2(fold.enrichment)]
res_go_up<-res_go[query=="up"]
fwrite(res_go_up,fp(out,"3E-res_go_bp_hto_up_padj0.05.csv"))
res_go_up[p.adjust<0.05]#243

#simplify with rrvgo
library(rrvgo)
simMatrix <- calculateSimMatrix(res_go_up[p.adjust<0.05]$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(res_go_up[p.adjust<0.05]$qvalue), res_go_up[p.adjust<0.05]$ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
head(reducedTerms)
redgo<-data.table(reducedTerms)

scatterPlot(simMatrix, reducedTerms,size = "score")

redgo[,n.by.cluster:=.N,by="cluster"]
redgo[cluster==1]
redgo[,top.term.cluster:=go==parent]
max(redgo$cluster) #31
sum(redgo$top.term.cluster) #31
redgo[,ID:=go]
res_red<-merge(res_go_up,redgo,by="ID")

res_red[,score.max:=max(score),by="cluster"]
res_red[top.term.cluster==T,top10clus:=rank(-score)<=15]
ggplot(res_red[top.term.cluster==T&top10clus==T&p.adjust<0.05],aes(x=term,y=-log10(p.adjust)))+
  geom_point(aes(size=n.overlap,col=pct.overlap.query))+
  scale_x_discrete(guide = guide_axis(angle=80),limits=res_red[query=="up"&top.term.cluster==T&top10clus==T][order(p.adjust)]$term)+
  scale_y_continuous()+
  scale_color_gradient2(low = "white",mid="grey",high = "black")
ggsave(fp(out,"3E-dotplot_pathway_HTO_GO_BP_up_top15_non_redundant_terms.pdf"),height=6)
fwrite(res_red[top.term.cluster==T&top10clus==T&p.adjust<0.05],fp(out,"3E-top_terms_clusters.csv"))
  
#or emaplot
library(clusterProfiler)
library(enrichplot)
go_up<-readRDS("outputs/08-HTO_signature/res_hto_signature_go_bp_up.rds")
pdf(fp(out,"3E-emapplot_go_bp_up_padj0.05_top50.pdf"),width = 12,height = 8)
emapplot(pairwise_termsim(go_up),showCategory = 50,cex_label_category=0.66)
dev.off()
#and save reduce terms


# terms_of_interest<-c("activation",
#                      "response",
#                      "stimulus",
#                      "stress",
#                      "differentiation",
#                      "hemopoiesis" ,
#                      "hematopoietic"   ,
#                     "proliferation",
#                      "cell cycle" )
# 
# res_go[,n.term.int:=sapply(Description,function(x)sum(str_detect(x,terms_of_interest)))]
# 
# res_go[p.adjust<0.001]
# #res_go[p.adjust<0.001][order(-fold.enrichment)]
# 
# res_go[p.adjust<0.001&n.term.int>0]
# 
# # ggplot(res_go[query=="up"],aes(x=log2FE,y=-log10(p.adjust)))+geom_point(aes(size=n.overlap,col=pct.overlap.query))+
# #   geom_label_repel(aes(label=ifelse(p.adjust<0.001&str_detect(Description,paste(terms_of_interest,collapse = "|")),Description,"")),
# #                    max.overlaps = 3000)+facet_wrap("query")
# 
# res_go[,wc:=str_count(Description,"\\S+")]
# res_go[wc>5,desc:=sapply(Description,function(x){
#   vec<-strsplit(x," ")[[1]]
#   vec_<-c(vec[1:4],"\n",vec[5:length(vec)])
#   return(paste(vec_,collapse = " "))
#   })]
# 
# res_go[wc<=5,desc:=Description]
# 
# ggplot(res_go[query=="up"][p.adjust<0.001],aes(x=desc,y=-log10(p.adjust)))+
#   geom_point(aes(size=n.overlap,col=pct.overlap.query))+
#   scale_x_discrete(guide = guide_axis(angle=80),limits=res_go[query=="up"][p.adjust<0.001][order(p.adjust)]$desc)+
#   scale_y_continuous(limits=c(0,7))+
#   scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 0.3,limits=c(0,0.7))
# 
# ggsave(fp(out,"3E-dotplot_pathway_HTO_GO_BP_up_padj0.001.pdf"),height = 6.6)
# 
# 
# ggplot(res_go[query=="up"][p.adjust<0.001&n.term.int>0],aes(x=desc,y=-log10(p.adjust)))+
#   geom_point(aes(size=n.overlap,col=pct.overlap.query))+
#   scale_x_discrete(guide = guide_axis(angle=80),limits=res_go[query=="up"][p.adjust<0.001&n.term.int>0][order(p.adjust)]$desc)+
#   scale_y_continuous(limits=c(0,5))+
#   scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 0.3,limits=c(0,0.7))
# 
# ggsave(fp(out,"3E-dotplot_pathway_HTO_GO_BP_up_padj0.001_selected.pdf"),height = 6.6)

# #kegg
# res_kegg_up<-fread("outputs/08-HTO_signature/res_hto_signature_kegg_up.csv")
# res_kegg_up[,n.overlap:=as.numeric(Count)]
# res_kegg_up[,n.query:=as.numeric(str_extract(GeneRatio,"[0-9]+$"))]
# res_kegg_up[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]
# res_kegg_up[,n.background:=as.numeric(str_extract(BgRatio,"[0-9]+$"))]
# 
# res_kegg_up[,pct.overlap.query:=n.overlap/n.gene_set]
# res_kegg_up[,pct.overlap.background:=n.gene_set/n.background]
# res_kegg_up[,fold.enrichment:=pct.overlap.query/pct.overlap.background]
# res_kegg_up[,log2FE:=log2(fold.enrichment)]
# 
# 
# ggplot(res_kegg_up[p.adjust<0.01],aes(x=Description,y=-log10(p.adjust)))+
#   geom_point(aes(size=n.overlap,col=pct.overlap.query))+
#   scale_x_discrete(guide = guide_axis(angle=80),limits=res_kegg_up[p.adjust<0.01][order(p.adjust)]$Description)+
#   scale_y_continuous(limits=c(0,8))+
#   scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 0.3,limits=c(0,0.7))
# 
# ggsave(fp(out,"3E-dotplot_pathway_HTO_KEGG_up_padj0.01.pdf"),height = 6.6)


#4-LGA vs CTRL####
#4A1-lGA_vs_Ctrl HTO / NOT
res_lc<-rbind(fread("outputs/07-LGA_vs_Ctrl_Basal/res_pseudobulkDESeq2_all_cbps.csv.gz")[,hto:=FALSE],
              fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_all_cbps.csv.gz")[,hto:=TRUE],fill=T)

ggplot(res_lc,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point(size=1)+
  facet_wrap("hto")+
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-4.5,4.5))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4A1-volcanos_all_cbps.pdf"),height = 6.6)

#4A-volcano HSC ####
res_lin<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")

res_hsc<-res_lin[lineage=="HSC"]
ggplot(res_hsc,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-4.5,4.5))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4A-volcano_HSC.pdf"),height = 6.6)

#MA plot
genes_of_interest<-c("SESN2","SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

ggplot(res_hsc,aes(y=log2FoldChange,x=baseMean,col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),

                   max.overlaps = 3000,
                   size=3)+
  scale_color_manual(values = c("grey","red")) +scale_x_log10()+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4A-ma_plot_HSC.pdf"))

#all progens
progens<-c("LT-HSC","HSC","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas")
res_lin<-res_lin[lineage%in%progens][,lineage:=factor(lineage,levels = progens)]
ggplot(res_lin,aes(y=log2FoldChange,x=baseMean,alpha=padj<0.05&abs(log2FoldChange)>0.5,
                   col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  facet_wrap("lineage")+
  scale_color_manual(values = c("grey","red")) +scale_x_log10()+
    scale_alpha_manual(values = c(0.6,1)) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(fp(out,"4A-ma_plot_all_lin.pdf"))


#4B-pathway_down ####

#gsea go bp
library(clusterProfiler)
library(enrichplot)
res_gsea_go<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.rds")
pdf(fp(out,"4B-emapplot_gsea_go_bp_dn_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_gsea_go),showCategory = 50,cex_label_category=0.66)
dev.off()


library(rrvgo)
go_dn<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res")
res_go_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
res_go_dn[p.adjust<0.05]#47
res_go_dn[,n.overlap:=as.numeric(Count)]
res_go_dn[,n.query:=as.numeric(str_extract(GeneRatio,"[0-9]+$"))]
res_go_dn[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]
res_go_dn[,n.background:=as.numeric(str_extract(BgRatio,"[0-9]+$"))]

res_go_dn[,pct.overlap.query:=n.overlap/n.gene_set]
res_go_dn[,pct.overlap.background:=n.gene_set/n.background]
res_go_dn[,fold.enrichment:=pct.overlap.query/pct.overlap.background]
res_go_dn[,log2FE:=log2(fold.enrichment)]
fwrite(res_go_dn[p.adjust<0.05],fp(out,"4B-res_go_bp_down_padj0.05.csv"))

pdf(fp(out,"4B-emapplot_go_bp_dn_padj0.05.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(go_dn),showCategory = 47,cex_label_category=0.66)
dev.off()
#and save reduce terms
simMatrix <- calculateSimMatrix(res_go_dn[p.adjust<0.05]$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(res_go_dn[p.adjust<0.05]$fold.enrichment, res_go_dn[p.adjust<0.05]$ID)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
redgo_dn<-data.table(reducedTerms)
redgo_dn[,top.term.cluster:=go==parent]
sum(redgo_dn$top.term.cluster) #15
redgo_dn[,ID:=go]

res_red_dn<-merge(res_go_dn,redgo_dn,by="ID")
res_red_dnf<-res_red_dn[p.adjust<0.05][top.term.cluster==T]
res_red_dnf#15
res_red_dnf$term
fwrite(res_red_dnf,fp(out,"4B-top_terms_clusters.csv"))

#dotplot
ggplot(res_red_dnf,aes(x=term,y=-log10(p.adjust)))+
  geom_point(aes(size=n.overlap,col=pct.overlap.query))+
  scale_x_discrete(guide = guide_axis(angle=80),limits=res_red_dnf[order(p.adjust)]$term)+
  scale_y_continuous()+
  scale_color_gradient2(low = "white",mid="grey",high = "black")
ggsave(fp(out,"4B-dotplot_pathway_HTO_GO_BP_dn_padj0.05_non_redundant_terms.pdf"),height=8)
       
# res_go_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
# res_go_dn[,n.overlap:=as.numeric(Count)]
# res_go_dn[,n.query:=as.numeric(str_extract(GeneRatio,"[0-9]+$"))]
# res_go_dn[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]
# res_go_dn[,n.background:=as.numeric(str_extract(BgRatio,"[0-9]+$"))]
# 
# res_go_dn[,pct.overlap.query:=n.overlap/n.gene_set]
# res_go_dn[,pct.overlap.background:=n.gene_set/n.background]
# res_go_dn[,fold.enrichment:=pct.overlap.query/pct.overlap.background]
# res_go_dn[,log2FE:=log2(fold.enrichment)]
# 
# terms_of_interest<-c("activation",
#                      "response",
#                      "stimulus",
#                      "stress",
#                      "differentiation",
#                      "hemopoiesis" ,
#                      "hematopoietic"   ,
#                     "proliferation",
#                      "cell cycle" )
# 
# res_go_dn[,n.term.int:=sapply(Description,function(x)sum(str_detect(x,terms_of_interest)))]
# 
# res_go_dn[p.adjust<0.05]
# #res_go_dn[p.adjust<0.001][order(-fold.enrichment)]
# 
# res_go_dn[p.adjust<0.05&n.term.int>0]
# 
# # ggplot(res_go_dn[query=="up"],aes(x=log2FE,y=-log10(p.adjust)))+geom_point(aes(size=n.overlap,col=pct.overlap.query))+
# #   geom_label_repel(aes(label=ifelse(p.adjust<0.001&str_detect(Description,paste(terms_of_interest,collapse = "|")),Description,"")),
# #                    max.overlaps = 3000)+facet_wrap("query")
# 
# res_go_dn[,wc:=str_count(Description,"\\S+")]
# res_go_dn[wc>5,desc:=sapply(Description,function(x){
#   vec<-strsplit(x," ")[[1]]
#   vec_<-c(vec[1:4],"\n",vec[5:length(vec)])
#   return(paste(vec_,collapse = " "))
#   })]
# 
# res_go_dn[wc<=5,desc:=Description]
# 
# ggplot(res_go_dn[p.adjust<0.05],aes(x=desc,y=-log10(p.adjust)))+
#   geom_point(aes(size=n.overlap,col=pct.overlap.query))+
#   scale_x_discrete(guide = guide_axis(angle=80),limits=res_go_dn[p.adjust<0.05][order(p.adjust)]$desc)+
#   scale_y_continuous(limits=c(0,3.5))+
#   scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 0.3,limits=c(0,0.7))
# 
# ggsave(fp(out,"4B-dotplot_pathway_HSC_down_GO_padj0.05.pdf"),height = 6.6,width = 12)
# 
# 
# ggplot(res_go_dn[p.adjust<0.05&n.term.int>0],aes(x=desc,y=-log10(p.adjust)))+
#   geom_point(aes(size=n.overlap,col=pct.overlap.query))+
#   scale_x_discrete(guide = guide_axis(angle=80),limits=res_go_dn[p.adjust<0.05&n.term.int>0][order(p.adjust)]$desc)+
#   scale_y_continuous(limits=c(0,3.5))+
#   scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 0.3,limits=c(0,0.7))
# 
# ggsave(fp(out,"4B-dotplot_pathway_HSC_down_GO_padj0.05_selected.pdf"),height = 6.6)

#kegg
res_kegg_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg_dn.csv")

res_kegg_dn[p.adjust<0.05]
res_kegg_dn[,n.overlap:=as.numeric(Count)]
res_kegg_dn[,n.query:=as.numeric(str_extract(GeneRatio,"[0-9]+$"))]
res_kegg_dn[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]
res_kegg_dn[,n.background:=as.numeric(str_extract(BgRatio,"[0-9]+$"))]

res_kegg_dn[,pct.overlap.query:=n.overlap/n.gene_set]
res_kegg_dn[,pct.overlap.background:=n.gene_set/n.background]
res_kegg_dn[,fold.enrichment:=pct.overlap.query/pct.overlap.background]
res_kegg_dn[,log2FE:=log2(fold.enrichment)]


# ggplot(res_kegg_dn[p.adjust<0.05],aes(x=Description,y=-log10(p.adjust)))+
#   geom_point(aes(size=n.overlap,col=pct.overlap.query))+
#   scale_x_discrete(guide = guide_axis(angle=80),limits=res_kegg_dn[p.adjust<0.05][order(p.adjust)]$Description)+
#   scale_y_continuous(limits=c(0,3.5))+
#   scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 0.3,limits=c(0,0.7))
# 
# ggsave(fp(out,"4B-dotplot_pathway_HSC_down_KEGG_padj0.05.pdf"))
# 


#4C-correl meth expr (boxplot)####
res_meth<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
res_merge<-merge(res_hsc,res_meth[,.(gene,gene_score_add)])
res_merge[,degs:=ifelse(padj<0.05&abs(log2FoldChange)>0.5,ifelse(log2FoldChange>0,"upregulated","downregulated"),"non-regulated")]
res_merge[,genes:=factor(degs,levels = c("non-regulated","upregulated","downregulated"))]
res_merge[,genes.:=paste0(genes,"\n(n=",.N,")"),by="genes"]

ggplot(res_merge)+
  geom_boxplot(aes(x=genes.,y=gene_score_add,fill=genes),outlier.shape = NA,col="bisque4")+
  scale_x_discrete(limits=unique(res_merge[order(genes)]$genes.))+
  scale_fill_manual(values = c("white","grey","black"))+
  coord_cartesian(ylim = c(0,1800))+theme_minimal()
ggsave(fp(out,"4C-boxplot_correl_meth_expr.pdf"))

res_merge[,p.up:=wilcox.test(gene_score_add[genes=='upregulated'],
                                                      gene_score_add[genes=='non-regulated'])$p.value]
res_merge[,p.down:=wilcox.test(gene_score_add[genes=='downregulated'],
                                                      gene_score_add[genes=='non-regulated'])$p.value]

unique(res_merge,by="p.down")
#  p.up       p.down
# 0.398549 2.373248e-05


#4D-MAP plot meth expr ####
quants<-quantile(res_mg$gene_score_add,1:100/100)
res_mg[,gs_percentile:=sum(gene_score_add>=quants),by="gene"]
res_meth_expr<-merge(res_mg,res_hsc_dn,by="gene")

ggplot(res_meth_expr,aes(y=log2FoldChange*-log10(padj.y),x=gs_percentile,col=padj.y<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj.y<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4D-plot_expr_and_meth_change_.pdf"))


#4E-overlap####
#pathways overlap
#[test] pathways analysis with same level
{
library(clusterProfiler)
library(org.Hs.eg.db)
# #for meth genes, take genes with genescore >= genescore where 90% of the genes have DMCs
# res_meth<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
# res_meth[,n.dmcs:=sum(pval<0.001&abs(meth.change)>25),by="gene"]
# res_mg<-unique(res_meth[order(-gene_score_add,cpg_score)],by="gene")
# res_mg[,dmcg:=n.dmcs>0]
# res_mg[,pct.dmcg_before:=sapply(1:.N,function(i)sum(res_mg[1:i]$dmcg)/i)]
# res_mg[,first.gs.thr:=max(gene_score_add[pct.dmcg_before<0.90])]
# res_mg[,last.gs.thr:=min(gene_score_add[pct.dmcg_before>0.90])]
# 
# res_mg$last.gs.thr[1] #1798.899
# nrow(res_mg[gene_score_add>1798.899])#332 genes
# 
# res_hto<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")
# 
# genes_list<-list(meth=res_mg[gene_score_add>1798.899]$gene,
#                  hsc_dn=res_hsc[padj<0.05&log2FoldChange<(-0.5)]$gene,
#                  hto_up=res_hto[padj<0.05&log2FoldChange>0.5]$gene)
# lapply(genes_list, length)
# # $meth
# # [1] 332
# # 
# # $hsc_dn
# # [1] 285
# # 
# # $hto_up
# # [1] 1075
# 
# pathways_list<-lapply(genes_list,function(genes){
#   res_enr<-enrichGO(bitr(genes,fromType = "SYMBOL",
#                              toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
#            OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
#            universe = bitr(intersect(res_mg$gene,res_hsc$gene),fromType = "SYMBOL",
#                              toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
#   res_enr<-data.table(as.data.frame(res_enr))
#   return(res_enr[p.adjust<0.05]$Description)
#   })
}
renv::install("ggvenn")
library(ggvenn)

#emaplot meth expr

meth<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")[p.adjust<0.05]
hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.csv.gz")
hsc_dn<-hsc_dn[order(pvalue)][1:50]
meth_dn<-merge(hsc_dn,meth,all.x=T,by=c("ID","Description"),suffixes = c("",".meth"))

res_hsc_dn<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.rds")

res_hsc_dnm<-res_hsc_dn
res_hsc_dnm@result<-data.frame(meth_dn,row.names = "ID")
head(res_hsc_dnm@result)
res_hsc_dnm@result$ID<-rownames(res_hsc_dnm@result)
pdf(fp(out,"4E-emapplot_overlap_methyl_on_expr_go_bp_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_hsc_dnm),
         color="p.adjust.meth",
         showCategory = 50,
         cex_label_category=0.66)
dev.off()


#venn
pathways_list<-list(
  meth=fread("outputs/03-pathway_analysis/res_gsea_go.csv")[p.adjust<0.001]$Description,
  hsc_lga=fread("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.csv.gz")[p.adjust<0.05]$Description,
  hto=fread("outputs/08-HTO_signature/res_gsea_go_bp.csv.gz")[p.adjust<0.05]$Description
  )
lapply(pathways_list, length)
# $meth
# [1] 637
# 
# $hsc_dn
# [1] 42
# 
# $hto_up
# [1] 68

ggvenn(
  pathways_list, 
  fill_color = c("grey","darkgreen","chartreuse3"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
  )
ggsave(fp(out,"4E-pathways_overlap_GO.pdf"))

#only meth and hsc dn
ggvenn(
  pathways_list[1:2], 
  fill_color = c("grey","darkgreen"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F
  )
ggsave(fp(out,"4E-pathways_overlap_without_hto_GO.pdf"))
cat(ps(Reduce(intersect,pathways_list),collapse = "\n"))


setdiff(intersect(pathways_list$meth,pathways_list$hsc_lga),
        pathways_list$hto)
# [1] "negative regulation of growth"                         
# [2] "regulation of cell growth"                             
# [3] "signal transduction by p53 class mediator"             
# [4] "regulation of protein serine/threonine kinase activity"
# [5] "regulation of body fluid levels"                       
# [6] "cellular response to chemical stress"                  
# [7] "hemostasis"


#DotPlot methylation pathway
path_comm<-Reduce(intersect,pathways_list[1:2])
path_comm
res_gsea_meth<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")
res_gsea_methf<-res_gsea_meth[Description%in%pathways_list$hsc_lga]
res_gsea_methf[,n.enriched:=length(tr(core_enrichment)),"ID"]
ggplot(res_gsea_methf,aes(x=Description,y=-log10(p.adjust)))+
  geom_point(aes(size=n.enriched,col=NES))+
  scale_x_discrete(guide = guide_axis(angle=80),limits=res_gsea_methf[order(p.adjust)]$Description)+
  scale_y_continuous()+
  scale_color_gradient2(low = "white",mid="grey",high = "black",midpoint = 1.3)
ggsave(fp(out,"4E-doplot_methylation_enrichment_go_bp_hsc_lga.pdf"),height = 8)

#dotplot heatmap DEGs DMCs, HTO
res_all<-Reduce(function(x,y)rbind(x,y,fill=T),list(
  meth=fread("outputs/03-pathway_analysis/res_gsea_go.csv")[,source:="meth"],
  hsc_lga=fread("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.csv.gz")[,source:="expr"],
  hto=fread("outputs/08-HTO_signature/res_gsea_go_bp.csv.gz")[,source:="hto"]
  ))
res_all[,n.enriched:=length(tr(core_enrichment)),"ID"]

library(rrvgo)
path_expr<-res_all[p.adjust<0.05&source=="expr"]$ID
simMatrix <- calculateSimMatrix(path_expr,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(res_all[p.adjust<0.05&source=="expr"]$p.adjust)*-log10(res_all[source=="meth"][path_expr,on="ID"]$p.adjust),
                   path_expr)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
redgo<-data.table(reducedTerms)

redgo[,top.term.cluster:=go==parent]
sum(redgo$top.term.cluster) #15
redgo[,ID:=go]

pathw_int<-redgo[(top.term.cluster)]$term
pathw_int<-intersect(pathw_int,res_all[source=="expr"&NES<0]$Description)

ggplot(res_all[p.adjust<0.05&Description%in%pathw_int&source!="hto"])+
  geom_point(aes(x=source,y=Description,size=n.enriched,col=-log10(p.adjust)))+
    scale_color_gradient2(low = "grey",mid="red",high = "darkred",midpoint = 3,limits=c(0,6))
ggsave(fp(out,"4E-doplot_compa_methylation_expression_enrichment_go_bp_hsc_lga_.pdf"),height = 8)


#SOME EXAMPLE OF KEY GENES FROM CONSERVED PATHWAYS

res_go_hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
res_go_meth<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")
path_comm
path_int<-c("negative regulation of growth",
            "negative regulation of cell growth",
            "myeloid cell differentiation",
            "regulation of protein stability")
genes_dn_bp<-res_go_hsc_dn[Description%in%path_int,tr(geneID,tradEntrezInSymbol = T),by="Description"]
colnames(genes_dn_bp)<-c("desc","gene")

genes_meth_bp<-res_go_meth[Description%in%path_int,tr(core_enrichment,tradEntrezInSymbol = T),by="Description"]
colnames(genes_meth_bp)<-c("desc","gene")
table(genes_meth_bp$desc)
genes_bp<-rbind(genes_dn_bp[,alter:="dn"],genes_meth_bp[,alter:="meth"])
genes_bp[,n_bp:=.N,by=.(gene,alter)]
genes_bp[,n_alter:=.N,by=.(gene,desc)]
genes_bp[n_alter==2]
genes_bp[,level_alter:=ifelse(.N==2,"both",alter),.(gene,desc)]
genes_bp[,level_alter:=factor(level_alter,levels = c("meth","dn","both"))]

gene_bp_mat<-dcast(genes_bp,desc~gene,value.var = "level_alter",fun.aggregate = function(x)ifelse(is.null(x),0,x))
gene_bp_mat<-as.matrix(data.frame(gene_bp_mat,row.names = "desc"))
gene_bp_mat[is.na(gene_bp_mat)]<-0
pheatmap::pheatmap(gene_bp_mat,color = c("white","purple","orange","red"),cluster_rows = T,cluster_cols = T,fontsize_col =  4)
#need filter for highly methylated genes
#at least 1 DMCs
res_meth<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
res_meth[,n.dmcs:=sum(pval<0.001&abs(meth.change)>25),by="gene"]
res_mg<-unique(res_meth[order(-gene_score_add,-cpg_score)],by="gene")
plot(density(res_mg$gene_score_add))
abline(v = 300)
res_mg[gene_score_add>300]
genes_meth_bpg<-merge(genes_meth_bp,res_mg[n.dmcs>0|gene_score_add>300])

genes_bp<-rbind(genes_dn_bp[,alter:="dn"],genes_meth_bpg[,.(desc,gene)][,alter:="meth"])
genes_bp[,n_bp:=.N,by=.(gene,alter)]
genes_bp[,n_alter:=.N,by=.(gene,desc)]
genes_bp[n_alter==2]
genes_bp[,level_alter:=ifelse(.N==2,"both",alter),.(gene,desc)]
genes_bp[,level_alter:=factor(level_alter,levels = c("meth","dn","both"))]

gene_bp_matf<-dcast(genes_bpf,desc~gene,value.var = "level_alter",fun.aggregate = function(x)ifelse(is.null(x),0,x))
gene_bp_matf<-as.matrix(data.frame(gene_bp_matf,row.names = "desc"))
gene_bp_matf[is.na(gene_bp_matf)]<-0
pdf(fp(out,"4D-heatmap_genes_common_bp_of_int_expr_meth.pdf"),width = 12,height = 2)
pheatmap::pheatmap(gene_bp_matf,color = c("white","purple","orange","red"),cluster_rows = T,cluster_cols = T,fontsize_col =  4,fontsize_row = 6)
dev.off()

#focus on negative regulation of cell growth
#plot meth / expr of this bp
quants<-quantile(res_mg$gene_score_add,1:100/100)
res_mg[,gs_percentile:=sum(gene_score_add>=quants),by="gene"]
res_meth_expr<-merge(res_mg,res_hsc_dn,by="gene")
genes_neg_growth<-FindGOGenes("negative regulation of cell growth")
res_me_ng<-res_meth_expr[gene%in%genes_neg_growth$hgnc_symbol]

ggplot(res_me_ng,aes(y=log2FoldChange*-log10(padj.y),x=gs_percentile,col=padj.y<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj.y<0.05&abs(log2FoldChange)>0.5,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4E-plot_expr_and_meth_change_neg_regul_cell_growth_genes.pdf"))
#most epigen altered genes of this bp
genes_of_int<-unique(genes_bp[desc=="negative regulation of cell growth"&level_alter=="both"]$gene)
genes_of_int#"SESN2"   "SEMA4A"  "SIRT1"   "SEMA7A"  "SERTAD3"

res_mg[gene%in%genes_of_int][order(-gene_score_add)]
res_mg[gene_score_add>=500]
res_meth[gene=="SESN2"][order(pval)][pval<0.01]
res_meth[in_eQTR==T&gene=="SESN2",tss_dist:=tss_dist-326501]
res_meth[gene=="SESN2"][order(pval)][pval<0.01]
res_meth[cpg_id==52791]
res_meth[gene=="SIRT1"][order(pval)][pval<0.01]
res_meth[gene=="SEMA7A"][order(pval)][pval<0.01]

ggplot(unique(res_meth[gene%in%c("SESN2","SIRT1","SEMA7A")&pval<0.01][order(gene_score_add,pval,-in_eQTR)],by=c("gene","cpg_id"))[order(gene_score_add,tss_dist)],
       by=c("gene","cpg_id"),)+
  geom_col(aes(x=factor(tss_dist),y=meth.change,fill=-log10(pval),col=in_eQTR))+
  facet_wrap("gene",scales = "free_x")+scale_fill_gradient(low = "white",high = "black",limits=c(2,3))+
ggsave(fp(out,"4E-barplot_meth_change_cpg_top_meth_genes_neg_regul_cell_growth.pdf"))
res_meth[gene=="HSPA1B"][order(pval)][pval<0.01]

go_ids<-lapply(path_comm,FindGO_ID) #GO:0030308
genes_bp<-lapply(go_ids,FindGOGenes)
genes_bp<-unique(genes_bp[go_id=="GO:0030308"]$hgnc_symbol)

#negative regulation of growth

res_hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]
res_meth<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
res_meth[,n.dmcs:=sum(pval<0.001&abs(meth.change)>25),by="gene"]
res_mg<-unique(res_meth[order(-gene_score_add,cpg_score)],by="gene")
res_mg[,gs_scaled:=scale(gene_score_add)]
res_merge<-merge(res_hsc_dn,res_mg,by="gene")
res_merge[,bp_of_int:=gene%in%genes_bp]

res_merge[,avg.dmcs:=mean(n.dmcs),by="bp_of_int"]

ggplot(res_merge)+geom_boxplot(aes(fill=avg.dmcs,y=gs_scaled,x=bp_of_int),outlier.shape = NA)+coord_cartesian(ylim = c(-1,3))

ggplot(res_merge)+geom_boxplot(aes(y=-log10(padj.x)*abs(log2FoldChange),x=bp_of_int),outlier.shape = NA)+coord_cartesian(ylim = c(0,0.5))



res_merge_bp<-res_merge[gene%in%genes_bp]
ggplot(res_merge_bp)+geom_point(aes(x=log2FoldChange,y=-log10(padj.x),size=n.dmcs,col=gs_scaled))

ggplot(res_merge_bp)+geom_point(aes(x=log2FoldChange*-log10(padj.x),y=gs_scaled,size=n.dmcs,col=padj.x<0.05&abs(log2FoldChange)>0.5))

ggplot(res_merge_bp[padj.x<0.05|n.dmcs>0])+geom_point(aes(x=gene,y=gs_scaled,size=n.dmcs,col=log2FoldChange))

#4E-genes negative regul cell growth
neg_cell_growth_genes<-c("SIRT1", "SESN2", "CDKN1A") 
res_meth_int<-res_meth[gene%in%neg_cell_growth_genes]
res_meth_int[pval<0.01]
unique(res_meth_int[order(-gene_score_add,pval)],by="gene")

# [test]-DEGs/DMGs overlap+ pathway enrichment
# res_meth<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
# res_meth[,n.dmcs:=sum(pval<0.001&abs(meth.change)>25),by="gene"]
# res_mg<-unique(res_meth[order(-gene_score_add,cpg_score)],by="gene")
# res_mg[,dmcg:=n.dmcs>0]
# res_mg[,pct.dmcg_before:=sapply(1:.N,function(i)sum(res_mg[1:i]$dmcg)/i)]
# res_mg[,first.gs.thr:=max(gene_score_add[pct.dmcg_before<0.75])]
# res_mg[,last.gs.thr:=min(gene_score_add[pct.dmcg_before>0.75])]
# 
# res_mg$first.gs.thr[1] #1084.112
# nrow(res_mg[gene_score_add>1084.112])#1408 genes
# 
# res_hto<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")
# 
# genes_list<-list(meth=res_mg[gene_score_add>1084.112]$gene,
#                  hsc_dn=res_hsc[padj<0.05&log2FoldChange<(-0.5)]$gene,
#                  hto_up=res_hto[padj<0.05&log2FoldChange>0.5]$gene)
# #
# ggvenn(
#   genes_list, 
#   fill_color = c("grey","darkgreen","chartreuse3"),
#   stroke_size = 0.5, set_name_size = 4,show_percentage = F
#   )
# ggsave(fp(out,"4E-genes_overlap.pdf"))
# 
# #pathway DMGs DEGs
# inter<-intersect(res_mg[gene_score_add>1084.112]$gene,res_hsc[padj<0.05&log2FoldChange<(-0.5)]$gene)
# inter#34
# res_enr<-enrichGO(bitr(inter,fromType = "SYMBOL",
#                            toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
#          OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
#          universe = bitr(intersect(res_mg$gene,res_hsc$gene),fromType = "SYMBOL",
#                            toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
# 
# res_enr<-data.table(as.data.frame(res_enr))
# res_enr[p.adjust<0.1]

#5 : Regulons####
#n_regulons
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")
length(regulons_list)#250
sum(!str_detect(names(regulons_list),"e$")) #106
#5B : TF activity alteration ####
cbps<-readRDS("outputs/10-SCENIC/cbps_with_regulons_activity.rds")
res_tf_diff<-fread("outputs/10-SCENIC/regulon_activity_lga_vs_ctrl_HTO_by_lineage.csv.gz")
res_tf_sig<-res_tf_diff[p_val_adj<0.001&lineage=="HSC"&hto==T&!str_detect(regulon,"e$")]
ggplot(res_tf_sig)+geom_density(aes(x=abs(avg_diff)))


res_tf_sig<-res_tf_diff[p_val_adj<0.001&abs(avg_diff)>0.03&lineage=="HSC"&hto==T&!str_detect(regulon,"e$")]

fwrite(res_tf_sig,fp(out,"5B-res_tf_activity_alteration_lga_hto_hsc_padj0.001_avg_diff0.3.csv"))

ggplot(res_hsc_htof)+
  geom_col(aes(x=regulon,y=avg_diff,fill=-log10(p_val_adj)))+
 scale_x_discrete(limits = res_hsc_htof[order(p_val_adj)]$regulon)

ggsave(fp(out,"5B-barplot_tf_change_hsc_lga_vs_ctrl_hto.pdf"))

mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tfs_alt<-res_hsc_htof$regulon
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[tfs_alt,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap=="HSC"&hto==T&regulon!="KLF2e"])+
  geom_boxplot(aes(x=regulon,y=activity,fill=group))+theme_bw()
ggsave(fp(out,"5B-boxplot_tf_activity_change_padj0.001_avg_diff0.3_lga_vs_ctrl_hsc_hto.pdf"))

#regulons activation in response to hto stress
res_hto_sig<-fread(fp(out,"3D-res_hto_signature_padj0.05_log2FC0.5.csv"))
res_hto_sig[gene%in%res_tf_sig$regulon] #all are upreg

res_tf_hto<-fread("outputs/10-SCENIC/regulon_activity_HTO_vs_not_by_lineage.csv.gz")
fwrite(res_tf_hto[avg_diff>0.03&p_val_adj<0.001],fp(out,"5B-res_hto_regulons_activation_by_lineage_padj0.001_avg_diff0.03.csv"))

res_tf_hto_sig<-res_tf_hto[avg_diff>0.03&p_val_adj<0.001]
res_tf_hto_sig[regulon%in%res_tf_sig$regulon&lineage=="HSC"] #all are upreg

mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[res_tf_sig$regulon,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap=="HSC"])+
  geom_boxplot(aes(x=regulon,y=activity,fill=hto))+theme_bw()
ggsave(fp(out,"5Bsupp-boxplot_hto_regulons_of_int_activation_in_hsc.pdf"))

tf_alt_act<-data.table(t(as.matrix(cbps@assays$SCT@data[res_tf_sig$regulon,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "expression" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap=="HSC"])+
  geom_boxplot(aes(x=regulon,y=expression,fill=hto))+theme_bw()
ggsave(fp(out,"5Bsupp-boxplot_hto_tf_of_int_activation_in_hsc.pdf"))

#all lin
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[res_tf_sig$regulon,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap%in%progens])+
  geom_boxplot(aes(x=regulon,y=activity,fill=hto))+
  facet_wrap("lineage_hmap")+
  theme_bw()
ggsave(fp(out,"5Bsupp-boxplot_hto_regulons_of_int_activation_by_lin.pdf"))

#5Bsupp2_regulons_annotation
library(ComplexHeatmap)
inp<-"outputs/08-HTO_signature/regulon_level/"
reg_enrich<-readRDS(fp(inp,"reg_enrich_gprofil.rds"))

res_bp<-fread(fp(inp,"res_regulons_enrich_for_BP_of_interest.csv"),select = c(1,11,3),col.names = c("regulon","term","padj"))
res_hsc<-fread(fp(inp,"res_regulons_enrich_hsc_state_signatures_geneset.csv"),select = c(9,1,8),col.names = c("regulon","term","padj"))
res_hscf<-res_hsc[term!='5FU_treatment']
terms_of_interest<-c("cellular response to stress",
                    "regulation of cellular response to stress" ,
                     "cellular response to growth factor stimulus",
                     "regulation of cellular response to growth factor stimulus" ,
                      "regulation of cell differentiation",
                     "regulation of hemopoiesis" ,
                     "hematopoietic progenitor cell differentiation"   ,
                    "cell population proliferation",
                     "regulation of cell population proliferation",
                     "regulation of cell cycle" )

enrich_mat<-as.matrix(data.frame(dcast(Reduce(rbind,list(res_bp[term%in%terms_of_interest],res_hscf)),formula = term~regulon),row.names = "term"))
enrich_mat<-(-log10(enrich_mat))
dim(enrich_mat)


#regulon order ~ gene overlap
# overlap_reg<-sapply(names(regulons_list),function(reg){
#   reg1<-regulons_list[[reg]]
#   return(sapply(regulons_list,function(reg2)length(intersect(reg1,reg2))/length(union(reg1,reg2))))
#   })
# clus<-hclust(as.dist(1-overlap_reg))
# reg_order<-rownames(overlap_reg)[clus$order]
# library(circlize)

col_fun4 = colorRamp2(c(0,3, 8), c("white", "grey","black"))
pdf(fp(out,"5Bsupp2-regulon_of_int_annotation.pdf"),height = 12)
Heatmap(t(enrich_mat[,res_tf_sig$regulon]), 
                col=col_fun4,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                cluster_rows = T,
                cluster_columns=T,
                name="-log10(djusted p-value)",
                row_title = "Regulon",
                column_title = "Gene sets")
dev.off()


#XX-degs dn pathways ####
degs<-fread("outputs/08-HTO_signature/by_lineage/res_pseudobulk_DESeq2_3replicates.csv.gz")
degs<-degs[padj<0.05&abs(log2FoldChange)>0.5]
table(degs$lineage)
degs[,ndegs:=.N,by="lineage"]
degs_pro<-degs[ndegs>60& lineage!="DC"]
degs_pro[,common_degs:=.N==5,by="gene"]
unique(degs_pro[common_degs==T]$gene)

res_kegg_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg_dn.csv")
res_kegg_dn<-res_kegg_dn[p.adjust<0.05]

res_go_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
res_go_dn[order(p.adjust)][p.adjust<0.05]$Description
res_go_dn<-res_go_dn[p.adjust<0.05]
terms_of_int<-rbind(res_go_dn[,source:="GO:BP"],res_kegg_dn[,source:="KEGG"])

terms_of_int[,intersection_size:=Count]
terms_of_int[,GeneRatio:=as.numeric(str_extract(GeneRatio,"[0-9]+"))/as.numeric(str_extract(GeneRatio,"[0-9]+$"))]

terms_of_int[,top20:=rank(p.adjust)<=20,by="source"]
terms_of_int$Description <- reorder(terms_of_int$Description, terms_of_int$p.adjust)
ggplot(terms_of_int[top20==T])+
  geom_point(aes(x=Description,y=-log10(p.adjust),size=Count,col=GeneRatio))+
  scale_x_discrete(guide = guide_axis(angle = 80))+facet_grid(~source,scales = "free",space = "free")

ggsave(fp(out,"dotplot_go_kegg_enrichment_HSC_LGA_down_top20_padj0.05.pdf"),width = 12,height = 8)



#figure ATAC
#
#a) figure UMAP predicted id####
library(Seurat)
library(Signac)
inp<-"outputs/14-DMCs_atac_integr/"
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")

DimPlot(atacs, group.by = "predicted.id", label = TRUE,reduction = "humap")
ggsave(fp(out,"Aa-cbps_atac_umap_predicted_lineage.pdf"))


#b) Lineage spe TF####
res_motif_lin<-fread(fp(inp,"enriched_motifs_by_lineage_specific_peaks.csv"))
#rm redundant motif
res_motif_lin[,top15:=rank(pvalue)<=15,by="lineage"]
res_motif_lin[,pvalue:=pvalue+10^-316]

res_motif_lin[,]
ggplot(res_motif_lin[pvalue<10^-50&lineage!="MPP/LMPP"])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  facet_grid(~lineage,scales = "free_x",space="free_x")+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))
ggsave(fp(out,"Ab-lineage_spe_tf_mlog10pval50.pdf"))

ggplot(res_motif_lin[pvalue<10^-50&lineage!="MPP/LMPP"])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  facet_wrap("lineage",scales = "free_x",ncol = 2)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))
ggsave(fp(out,"Ab-lineage_spe_tf_mlog10pval50_v2.pdf"))
fwrite(res_motif_lin,fp(out,"enriched_motifs_by_lineage_specific_peaks.csv"))

#c) DMCs enrichment in HSC peaks####
res_dmcs_peaks_enr<-fread(fp(inp,"res_dmcs_peaks_enrichment_for_lineage_specific_peaks.csv"))
res_dmcs_peaks_enr[,percent.observed:=precision*100]
ggplot(res_dmcs_peaks_enr)+
  geom_point(aes(x=query,size=n.overlap,col=-log10(padj),y=percent.observed))+
  scale_color_gradient(low = "grey",high = "red")+
    scale_y_continuous(expand = c(0,1))+scale_x_discrete(guide = guide_axis(angle = 0))+
    theme(axis.title.x =element_blank())

ggsave(fp(out,"Ac-dmcs_enrichment_in_lineage_spe_peaks.pdf"))


#d) TF DMCs HSC peaks vs CpG peaks####
res_motif_hsc_dmcs<-fread(fp(inp,"res_tf_motif_enrichment_in_DMCs_containing_vs_CpGs_containing_hsc_peaks.csv"))
res_motif_hsc_dmcs[,pvalue:=pvalue+10^-316]
ggplot(res_motif_hsc_dmcs[1:15])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60),
                                                                       limits=res_motif_hsc_dmcs[1:15][order(pvalue,-fold.enrichment)]$motif.name)+
  scale_size(limits = c(1,2.7))+
  theme(axis.text=element_text(size=9))

ggsave(fp(out,"Ad-dotplot_tf_motif_dmcs_vs_cpg_hsc_peaks_top15.pdf"))


#e)valid trs alteration méca####
res_or_degs_dmcs_tf<-fread(fp(inp,"res_degs_enrichment_in_dmcs_tf_hsc_peaks.csv"))

res_or_degs_dmcs_tf[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_degs_dmcs_tf[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_degs_dmcs_tf[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_degs_dmcs_tf[,degs_peaks:=paste0(query,"\n(n=",n.query,")")]

ggplot(res_or_degs_dmcs_tf)+geom_col(aes(y=log2(fold.enrichment),x=candidat_peaks,fill=padj_b),position = "dodge")+
  facet_wrap("degs_peaks")+
  scale_fill_manual(values=c("grey","orange","red","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))
ggsave(fp(out,"Ae-barplot_dmcs_tf_peaks_enrichment_in_degs_peaks.pdf"))

#f) diff accessibility lga vs ctrl HSC####
inp<-"outputs/15-chromatin_change_LGA_vs_Ctrl/"
#volcano

#tf
hsc_lga_tf_all<-fread(fp(inp,"motif_enrichment_in_peaks_up_and_down_lga_vs_ctrl_hsc.csv.gz"))

hsc_lga_tf_all[,top20:=rank(pvalue)<=20,by="accessibility"]
hsc_lga_tf_all_20<-hsc_lga_tf_all[top20==T]
hsc_lga_tf_all_20[,motif.name:=factor(motif.name,levels =hsc_lga_tf_all_20[order(pvalue)]$motif.name,ordered = T )]
ggplot(hsc_lga_tf_all_20)+
  geom_point(aes(x=motif.name,size=fold.enrichment,y=-log10(pvalue),col=percent.observed))+
  facet_wrap("accessibility",scales = "free_x",ncol = 2)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(0,100),midpoint = 50)+
  scale_y_continuous(limits=c(0,55),expand = c(0,0))+
  scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))
ggsave(fp(out,"Ag-dotplot_TF_motif_enrichment_in_up_or_down_peaks_LGA_vs_Ctrl_HSC.pdf"))

#DMCs / DEGs enrichment
res_or_da_peaks_dmcs_degs<-fread(fp(inp,"res_da_peaks_enrichment_in_dmcs_degs_hsc_peaks.csv"))

res_or_da_peaks_dmcs_degs[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_da_peaks_dmcs_degs[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_da_peaks_dmcs_degs[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_da_peaks_dmcs_degs[,da_peaks:=paste0(query,"\n(n=",n.query,")")]
res_or_da_peaks_dmcs_degs[,candidat_peaks:=factor(candidat_peaks,levels =res_or_da_peaks_dmcs_degs[query=="down_peaks"][order(fold.enrichment)]$candidat_peaks )]

ggplot(res_or_da_peaks_dmcs_degs)+geom_col(aes(y=log2(fold.enrichment),x=candidat_peaks,fill=padj_b),position = "dodge")+
  facet_wrap("da_peaks")+
  scale_fill_manual(values=c("grey","orange","red","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))
ggsave(fp(out,"Ah-barplot_dmcs_degs_peaks_enrichment_in_hsc_da_peaks.pdf"))

#SUPP####
#distrib predic lineage
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

ggplot(atacs@meta.data)+geom_bar(aes(x=predicted.id))

#can remove badly predict lineage ?
ggplot(atacs@meta.data)+geom_bar(aes(x=predicted.id))



