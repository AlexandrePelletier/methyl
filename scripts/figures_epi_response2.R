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
fwrite(markers[order(cell_type)],fp(out,"markers_hematopo_lineage.csv"))

markers[gene%in%c("CTSG","AZU1")]
markers[gene%in%c("CTSG","AZU1")]


markers[cell_type=="ErP"]

markers<-fread("outputs/05-make_hematomap/hematomap_ctrls_sans_stressSCT_Leiden_res0.6_markers.csv.gz")
markers[gene%in%c("CDK6","MLLT3")]
m_of_int<-markers[score>=3][MarqueurPops!=""][order(lineage)][cell_type!="18"]


m_anno<-fread("outputs/05-make_hematomap/markers_subpop")
setdiff(m_anno$marker,m_of_int$gene) #"ID2"   "DUSP2"  "FOS" "CDK6"
markers[gene%in%setdiff(m_anno$marker,m_of_int$gene)]
m_of_int<-rbind(m_of_int,markers)

m_of_int

fwrite(,fp(out,"markers_clusters_subpop.csv"))

m_lin<-fread("outputs/05-make_hematomap/markers_lineage_annotated.csv.gz")
m_lin


m_of_int<-rbind(m_lin[cluster!="18"][order(-score,p_val_adj)][,.SD[1:5],by="cluster"][,topspe:=T],m_lin[cluster!="18"][order(-score,p_val_adj)][MarqueurPops!=""][,.SD[1:5],by="cluster"][,topanno:=T],fill=T)[order(cluster)]
m_of_int<-unique(m_of_int[,topspe:=gene%in%gene[topspe==T],by="cluster"][,topanno:=gene%in%gene[topanno==T],by="cluster"])
m_of_int<-m_of_int[!is.na(gene)]
m_of_int[,MarqueurPops:=str_replace_all(MarqueurPops,";","|")]
fwrite(m_of_int,fp(out,"3B-markers_lineage.csv"))


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
DimPlot(hmap,group.by = "cell_type",label=T)

FeaturePlot(hmap,c("IGLL1"))


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

#pca
avg_expr<-AverageExpression(cbps,group.by = "sample_hto",assays = "SCT",slot = "data")$SCT
pca_expr<-prcomp(t(avg_expr))

mtd<-data.table(cbps@meta.data)

mtd[,seq_depth:=sum(nCount_RNA),by="sample_hto"]
vars_of_int<-c("group","hto","group_hto","sex","batch","seq_depth")
res_pc_mtd<-merge(data.table(pca_expr$x,keep.rownames = "sample_hto"),
                  unique(mtd,by="sample_hto"))

p1<-ggplot(res_pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=group))+
  theme_classic()
p2<-ggplot(res_pc_mtd)+geom_point(aes(x=PC1,y=PC2,col=hto))+
  scale_color_manual(values = c("grey","black"))+
  theme_classic()

p1+p2+plot_layout(guides = "collect")
ggsave(fp(out,"3D-pcplot_hto_effect.pdf"))


res_pc_mtd2<-merge(melt(data.table(pca_expr$x,keep.rownames = "sample_hto"),id.vars = "sample_hto",variable.name = "PC",value.name = "coordinate"),
                  unique(mtd,by="sample_hto"))

res_ps_pcs<-Reduce(rbind,lapply(paste0("PC",1:10),function(pc){
  pc_vals<-res_pc_mtd2[PC==pc]$coordinate
  vars_list<-res_pc_mtd2[PC==pc][,.SD,.SDcols=vars_of_int]
  pvals<-sapply(vars_list,function(var)anova(lm(pc_vals~var))$`Pr(>F)`[1])
  r2s<-sapply(vars_list,function(var)summary(lm(pc_vals~var))$r.squared)
  res<-data.table(PC=pc,var=vars_of_int,pval=pvals,r2=r2s)
  return(res)
  }))

res_hto_group_pcs<-Reduce(rbind,lapply(paste0("PC",1:10),function(pc){
  pc_vals<-res_pc_mtd2[PC==pc]$coordinate
  pval<-anova(lm(pc_vals~res_pc_mtd2[PC==pc]$group:res_pc_mtd2[PC==pc]$hto))$`Pr(>F)`[1]
  r2<-summary(lm(pc_vals~res_pc_mtd2[PC==pc]$group:res_pc_mtd2[PC==pc]$hto))$r.squared
  res<-data.table(PC=pc,var="group:hto",pval=pval,r2=r2)
  return(res)
  }))
res_ps_pcs<-rbind(res_ps_pcs,res_hto_group_pcs)

res_ps_pcs[,pct.pc:=pctPC(pca_expr)[as.numeric(str_extract(PC,"[0-9]+"))]]
res_ps_pcs[,PCpct:=paste0(PC," (",round(pct.pc*100,1),"%)")]

ggplot(res_ps_pcs[var!="group_hto"])+
  geom_point(aes(x=PCpct,y=var,col=-log10(pval),size=r2))+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",midpoint = 5)+
  scale_x_discrete(limits=unique(res_ps_pcs[order(-pct.pc)]$PCpct))+
  theme_minimal()
ggsave(fp(out,"3D-dotplot_hto_effect.pdf"))


res_ps_pcs[var%in%c("group_hto","group:hto")&PC=="PC1"]

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
fwrite(res_gsea_go_dt[NES>0][p.adjust<0.05],fp(out,'3E-res_hto_up_gsea_go_padj0.05.csv'))

#kegg
library(clusterProfiler)
library(enrichplot)
res_gsea_kegg<-readRDS("outputs/08-HTO_signature/res_gsea_kegg.rds")

res_gsea_kegg_up<-res_gsea_kegg
res_gsea_kegg_up@result<-subset(res_gsea_kegg@result,NES>0)
pdf(fp(out,"3E-emapplot_gsea_kegg_up_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_gsea_kegg_up),
         showCategory = 50,
         cex_label_category=0.66)
dev.off()

res_gsea_kegg_dt<-fread("outputs/08-HTO_signature/res_gsea_kegg.csv.gz")

fwrite(res_gsea_kegg_dt[NES>0][order(pvalue)][1:50],fp(out,'3E-res_hto_up_gsea_kegg_top50.csv'))


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
res_hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]
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

#emaplot meth expr

meth<-fread("outputs/03-pathway_analysis/res_gsea_go_bp_all.csv")[p.adjust<0.05]

meth[,mlog10padj.meth:=-log10(p.adjust)]
meth[mlog10padj.meth>10,mlog10padj.meth:=10]

hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.csv.gz")
hsc_dn<-hsc_dn[order(pvalue)][1:50]
meth_dn<-merge(hsc_dn,meth,all.x=T,by=c("ID","Description"),suffixes = c("",".meth"))

res_hsc_dn<-readRDS("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.rds")

res_hsc_dnm<-res_hsc_dn
res_hsc_dnm@result<-data.frame(meth_dn,row.names = "ID")
res_hsc_dnm@result$ID<-rownames(res_hsc_dnm@result)

pdf(fp(out,"4E-emapplot_overlap_methyl_on_expr_go_bp_top50.pdf"),width = 10,height = 6)
emapplot(pairwise_termsim(res_hsc_dnm),
         color="mlog10padj.meth",
         showCategory = 50,
         cex_label_category=0.66)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",midpoint = 5,limits=c(-log10(0.05),10))
dev.off()


#venn
#renv::install("ggvenn")
library(ggvenn)

pathways_list<-list(
  meth=fread("outputs/03-pathway_analysis/res_gsea_go_bp_all.csv")[p.adjust<0.001]$Description,
  hsc_lga=fread("outputs/09-LGA_vs_Ctrl_Activated/res_gsea_go_bp.csv.gz")[p.adjust<0.05]$Description,
  hto=fread("outputs/08-HTO_signature/res_gsea_go_bp.csv.gz")[p.adjust<0.05]$Description
  )
lapply(pathways_list, length)
# $meth
# [1] 715/ 6059
# 
# $hsc_lga
# [1] 46 / 4999
# 
# $hto
# [1] 90 / 5173

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
10/46#22%
phyper(q = 10-1,m = 715,n = 6059-715,k = 46,lower.tail = F)
ggsave(fp(out,"4E-pathways_overlap_meth_expr.pdf"))
cat(paste0(Reduce(intersect,pathways_list[1:2]),collapse = "\n"))
# regulation of growth
# regulation of cell growth
# negative regulation of growth
# regulation of mitotic cell cycle
# regulation of protein serine/threonine kinase activity
# signal transduction by p53 class mediator
# regulation of body fluid levels
# cellular response to chemical stress
# regulation of DNA-binding transcription factor activity
# negative regulation of apoptotic signaling pathway


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
res_go_meth<-fread("outputs/03-pathway_analysis/res_gsea_go_bp_all.csv")
path_comm
path_int<-c("negative regulation of growth",
            "negative regulation of cell growth",
            "regulation of growth",
            "regulation of mitotic cell cycle",
            "cellular response to chemical stress")
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
genes_meth_bpg<-merge(genes_meth_bp,res_mg[n.dmcs>0|gene_score_add>500])

genes_bp<-rbind(genes_dn_bp[,alter:="dn"],genes_meth_bpg[,.(desc,gene)][,alter:="meth"])
genes_bp[,n_bp:=.N,by=.(gene,alter)]
genes_bp[,n_alter:=.N,by=.(gene,desc)]
genes_bp[n_alter==2]
genes_bp[,level_alter:=ifelse(.N==2,"both",alter),.(gene,desc)]
genes_bp[,level_alter:=factor(level_alter,levels = c("meth","dn","both"))]

gene_bp_matf<-dcast(genes_bp,desc~gene,value.var = "level_alter",fun.aggregate = function(x)ifelse(is.null(x),0,x))
gene_bp_matf<-as.matrix(data.frame(gene_bp_matf,row.names = "desc"))
gene_bp_matf[is.na(gene_bp_matf)]<-0
pdf(fp(out,"4D-heatmap_genes_common_bp_of_int_expr_meth.pdf"),width = 12,height = 2)
pheatmap::pheatmap(gene_bp_matf,color = c("white","purple","orange","red"),cluster_rows = T,cluster_cols = T,fontsize_col =  4,fontsize_row = 6)
dev.off()

#focus of regulation of growth
#plot meth / expr of this bp
quants<-quantile(res_mg$gene_score_add,1:100/100)
res_mg[,gs_percentile:=sum(gene_score_add>=quants),by="gene"]
res_meth_expr<-merge(res_mg,res_hsc_dn,by="gene")
reg_growth_modul<-c("regulation of growth","regulation of cell growth","negative regulation of growth","negative regulation of cell growth")
genes_reg_growth<-FindGOGenes(reg_growth_modul)
genes_reg_growth#270
res_me_ng<-res_meth_expr[gene%in%genes_reg_growth$hgnc_symbol]
res_me_ng#153
ggplot(res_me_ng,aes(y=log2FoldChange*-log10(padj.y),x=gs_percentile,col=padj.y<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj.y<0.05&abs(log2FoldChange)>0.5,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4E-plot_expr_and_meth_change_genes_of_regul_growth_modul.pdf"))


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
#5A-TF motif enricment in DMCS (homer)
res<-fread("outputs/03B-motif_analysis/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[padj<0.05]
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]

res_perm<-fread("outputs/03B-motif_analysis/res_known_motif_all_perm_bgrandom.csv")
res_perm[,motif:=str_remove(motif,"/Homer")]

res[,permut:=0]
res_perm_merge<-merge(res,res_perm,all=T)
res_perm_merge[,p.perm:=sum(pct_dmc_with_motif[permut==0]<=pct_dmc_with_motif[permut!=0])/sum(permut!=0,na.rm = T),by=.(motif)]

res_perm_merge[padj<=0.05&n_dmc_with_motif>30&p.perm<0.05&permut==0]
res_perm_<-res_perm_merge[permut==0]

res_perm_[padj<=0.05] #26
res_perm_[padj<=0.05&n_dmc_with_motif>30] #26
res_perm_[padj<=0.05&n_dmc_with_motif>30&p.perm<0.01] #26

ggplot(res_perm_[padj<=0.05&pct_dmc_with_motif>1&p.perm<0.01])+
  geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_color_gradient(high ="brown1" ,low = scales::muted("red"))+
  scale_x_discrete(limits=res_perm_[padj<=0.05&pct_dmc_with_motif>1&p.perm<0.01][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust=1))


res_perm_[,motif.name:=sapply(motif,function(x)strsplit(x,"/")[[1]][1])]
ggplot(res_perm_[padj<=0.05&p.perm<0.01&fold.enrichment>1.30][order(pval)])+
  geom_point(aes(x=motif.name,y=-log10(pval),col=fold.enrichment,size=n_dmc_with_motif))+
  scale_x_discrete(limits=res_perm_[padj<=0.05&p.perm<0.01&fold.enrichment>1.30][order(pval)]$motif.name)+
  scale_color_gradient2(low = "white",mid = "grey",high = "black",midpoint = 1.25,limits=c(1,2.1))+
  scale_y_continuous(limits = c(0,25))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust=1))

ggsave(fp(out,"5A-motif_enrichment_homer_own_background_no_norm_padj0.05_p.perm0.01_fold.enrichment1.3.pdf"))

setdiff(res_perm_[padj<=0.05&p.perm<0.01]$motif,res_perm_[padj<=0.05&pct_dmc_with_motif>1&p.perm<0.01]$motif)
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

ggplot(res_tf_sig)+
  geom_col(aes(x=regulon,y=avg_diff,fill=-log10(p_val_adj)))+
 scale_x_discrete(limits = res_tf_sig[order(p_val_adj)]$regulon)+
  scale_fill_gradient(low="white",high="black",limits=c(10,60))

ggsave(fp(out,"5B-barplot_tf_change_hsc_lga_vs_ctrl_hto.pdf"))

mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tfs_alt<-res_hsc_htof$regulon
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[tfs_alt,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap=="HSC"&hto==T&regulon!="KLF2e"])+
  geom_boxplot(aes(x=regulon,y=activity,fill=group))+theme_bw()
ggsave(fp(out,"5B-boxplot_tf_activity_change_padj0.001_avg_diff0.3_lga_vs_ctrl_hsc_hto.pdf"))

#calc FC diff threshold
mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tfs<-unique(res_tf_diff$regulon)
tf_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[tfs,])),keep.rownames = "cell")
tf_act<-melt(tf_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_act<-merge(tf_act,mtd)
tf_act_hschto<-tf_act[lineage_hmap=="HSC"&hto==T]
tf_act_hschto[,FoldChange:=mean(activity[group=="lga"])/mean(activity[group=="ctrl"]),by="regulon"]

max(tf_act_hschto[regulon%in%unique(res_tf_sig$regulon)]$FoldChange) #10%
max(tf_act$activity)


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
                name="-log10(adjusted p-value)",
                row_title = "Regulon",
                column_title = "Gene sets")
dev.off()



#figure 6####
#
#6A : figure UMAP lineage####
library(Seurat)
library(Signac)
inp<-"outputs/14-DMCs_atac_integr/"
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")

DimPlot(atacs, group.by = "predicted.id", label = TRUE,reduction = "humap")
DimPlot(subset(atacs,predicted.id=="MPP/LMPP"&prediction.score.max<0.55,invert=T), group.by = "predicted.id", label = FALSE,reduction = "humap")

DimPlot(atacs, group.by = "seurat_clusters", label = TRUE,reduction = "humap")
mtd<-data.table(atacs@meta.data,keep.rownames = "bc")
mtd_pr<-melt(mtd,id.vars = c("bc","seurat_clusters"),value.name = "prediction.score",variable.name = "lineage",
             measure.vars = colnames(mtd)[str_detect(colnames(mtd),"prediction.score")])

mtd_pr[,lineage:=str_remove(lineage,"prediction.score.")]
mtd_pr<-mtd_pr[lineage!="max"]
ggplot(mtd_pr)+geom_boxplot(aes(y=prediction.score,x=seurat_clusters))+
facet_wrap("lineage")
head(atacs[[]])

VlnPlot(atacs,features = "prediction.score.max",group.by = "predicted.id",pt.size = 0)

DefaultAssay(atacs)<-"ATAC"
atacs <- FindClusters(object = atacs, verbose = FALSE, algorithm = 3,resolution = 0.8 )
DimPlot(object = atacs, reduction = "humap",label = TRUE) + NoLegend()
ggplot(mtd)+geom_bar(aes(x=seurat_clusters,fill=predicted.id),position = 'fill')

Idents(atacs)<-"seurat_clusters"

new.id<-c("MPP/LMPP",
  "HSC",
  "HSC",
  "MPP/LMPP",
  "MPP/LMPP",
  "Myeloid",
  "Erythro-Mas",
  "Erythro-Mas",
  "MPP/LMPP",
  "Lymphoid",
  "MPP/LMPP",
  "18",
  "Lymphoid",
  "B cell",
  "DC",
  "T cell")
names(new.id)<-levels(atacs)
atacs<-RenameIdents(atacs,new.id)
atacs[["clusters_anno"]]<-Idents(atacs)
atacs<-subset(atacs,seurat_clusters!="11")
DimPlot(atacs, group.by =c("seurat_clusters","predicted.id") ,label = TRUE,reduction = "humap",pt.size = 0.2,label.size = 3)

ggsave(fp(out,"6A-cbps_atac_umaps.pdf"))

DimPlot(atacs, group.by = "dataset", label = TRUE,reduction = "humap")

#n peaks spe by lin

res_peaks<-fread("outputs/14-DMCs_atac_integr/peaks_markers_lineage.csv.gz")
res_peaks[p_val_adj<0.001]

progens<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")
res_peaksf<-res_peaks[cluster%in%progens]
res_peaksf[,lineage:=factor(cluster,levels =progens )]
table(res_peaksf[p_val_adj<0.001]$lineage)
    # LT-HSC         HSC    MPP/LMPP     Myeloid    Lymphoid Erythro-Mas 
    #       2         465           1         359         622        1932 
ggplot(res_peaksf[p_val_adj<0.001])+geom_bar(aes(x=lineage,fill=lineage))
ggsave(fp(out,"6Asupp_barplot_n_peak_spe_by_lineage_padj0.001.pdf"))
fwrite(res_peaksf[p_val_adj<0.001],"outputs/6Asupp-res_peak_spe_by_lineage_padj0.001.csv.gz")


#6B-Lineage spe TF####
res_motif_lin<-fread(fp(inp,"enriched_motifs_by_lineage_specific_peaks.csv"))
resf<-res_motif_lin[pvalue<10^-50&lineage!="MPP/LMPP"]

resf[,motif.name:=factor(motif.name,levels=unique(motif.name[order(pvalue)]))]

res_motif_lin[,]
ggplot(resf)+
  geom_point(aes(x=motif.name,col=fold.enrichment,y=-log10(pvalue),size=observed))+
  facet_grid(~lineage,scales = "free",space="free_x")+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred")+
  scale_y_continuous()+scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))

ps<-lapply(unique(resf$lineage),function(lin)ggplot(resf[lineage==lin])+
  geom_point(aes(x=motif.name,col=fold.enrichment,y=-log10(pvalue),size=observed))+
  scale_color_gradient2(low = "white",mid = "darkgrey",high = "black",midpoint = 2,limits=c(1,8))+
  scale_y_continuous(expand = c(0.2,0.2))+
    scale_x_discrete(guide = guide_axis(angle = 60))+
    scale_size_continuous(limits = c(100,1200))+
  theme(axis.text=element_text(size=8))+ggtitle(lin))

wrap_plots(ps)+plot_layout(guides = 'collect')
ggsave(fp(out,"6B-lineage_spe_tf_mlog10pval50.pdf"),height = 8)


# DMCs enrichment in HSC peaks####
res_dmcs_peaks_enr<-fread(fp(inp,"res_dmcs_peaks_enrichment_for_lineage_specific_peaks.csv"))
res_dmcs_peaks_enr[,percent.observed:=precision*100]
ggplot(res_dmcs_peaks_enr)+
  geom_point(aes(x=query,size=n.overlap,col=-log10(padj),y=percent.observed))+
  scale_color_gradient(low = "grey",high = "red")+
    scale_y_continuous(expand = c(0,1))+scale_x_discrete(guide = guide_axis(angle = 0))+
    theme(axis.title.x =element_blank())

ggsave(fp(out,"6C-dmcs_enrichment_in_lineage_spe_peaks.pdf"))



# TF DMCs HSC peaks vs CpG peaks####
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


#6C: valid trs alteration méca####
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

#6D diff accessibility lga vs ctrl HSC####
inp<-"outputs/15-chromatin_change_LGA_vs_Ctrl/"
#volcano
res_da<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_hsc_logFC0.csv.gz")
res_da_xy<-res_da[!str_detect(peak,"X|Y")]
peak_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peak_anno[,peak:=query_region]
res_da_xy_anno<-merge(res_da_xy,peak_anno)
genes_of_interest<-c("LMNA","FOS","KLF13","NFKBIA","CDCA4","JUNB",'SAMD1',"GADD45B","SP3","GSK3A","WNT10A","FOXO3","E2F1",
                     "SGK1","SENP2","KIF2","GFPT2")

res_da_xy_anno[p_val_adj<0.001&avg_log2FC<(-0.25)]$gene_name
res_da_xy_anno[p_val_adj<0.001&avg_log2FC>0.25][order(p_val_adj)]$gene_name

res_da_xy_anno[,peak_gene:=paste0(peak,"(",gene_name,")")]
ggplot(res_da_xy_anno,aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.001&abs(avg_log2FC)>0.25))+
  geom_point(size=1) +
  geom_label_repel(aes(label=ifelse(p_val_adj<0.001&abs(avg_log2FC)>0.25&gene_name%in%genes_of_interest,peak_gene,"")),
                    label.size = NA,
                   max.overlaps = 3000,
                   size=2)+
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-1,1))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"6D-volcano_da_peaks_LGA_vs_Ctrl_HSC.pdf"),height = 6.6)
fwrite(res_da_xy_anno[p_val_adj<0.001&abs(avg_log2FC)>0.25],fp(out,"6D-res_da_peaks_LGA_vs_Ctrl_padj0.001_log2FC0.25.csv"))


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
ggsave(fp(out,"6D-dotplot_TF_motif_enrichment_in_up_or_down_peaks_LGA_vs_Ctrl_HSC_top20.pdf"))
fwrite(hsc_lga_tf_all[pvalue<10^-6],fp(out,"6D-res_TF_motif_enrichment_in_up_or_down_peaks_LGA_vs_Ctrl_HSC_pval10m6.csv"))

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


#6E : network final
#genescore 500 = top cb ?
res_m<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
quantile(res_m$gene_score_add,0.8+(1:10/100))
res_m[gene_score_add>500]
4289/24857 #0.172
1-0.172 #0.828
tfstarget<-fread("outputs/16-GRN_final/egr1_KLF2_KLF4_network_tf_target_interactions.csv")
res_mf<-res_m[gene%in%unique(union(tfstarget$target,tfstarget$tf))]

setdiff(res_mf[gene_score_add>500]$gene,res_mf[gene_score_add>572]$gene)
setdiff(res_mf[gene_score_add>quantile(gene_score_add,0.8)]$gene,res_mf[gene_score_add>500]$gene)

#82% top18% 

#6F : genome track JUNB, SOCS3, GADD45B, EGR1, LMNA####
#try with JUNB 
library(Seurat)
library(Signac)
inp<-"outputs/14-DMCs_atac_integr/"
inp2<-"outputs/15-chromatin_change_LGA_vs_Ctrl/"

atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
progens<-c("LT-HSC","HSC","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas")

peaks_da<-fread(fp(inp2,"differential_peaks_accessibility_lga_vs_ctrl_hsc_without_xy.csv.gz"))
peaks_da_anno<-ClosestFeature(atacs,peaks_da[p_val_adj<0.001&abs(avg_log2FC)>0.25]$peak)
peaks_da_anno<-data.table(peaks_da_anno)
peaks_da_anno[,peak:=query_region]
peaks_da_anno<-merge(peaks_da,peaks_da_anno,by="peak")

peaks_da_anno[gene_name=="JUNB"]

#showing peaks accessibility in CTRL LGA, gene anno
atacs[["group"]]<-ifelse(atacs$dataset%in%c("cbp.atac1","cbp.atac3"),"ctrl","lga")
atacs[["lineage.group"]]<-paste(atacs$predicted.id,atacs$group,sep=".")

annotations<-readRDS("../atac/ref/gene_annotations_hg38_GRanges.rds")
Annotation(atacs) <- annotations

atacs@assays$lin_peaks@motifs<-readRDS("outputs/14-DMCs_atac_integr/atacs_cbps_lin_peaks_motif_object.rds")
atacs_hsc<-subset(atacs,predicted.id=="HSC")
Idents(atacs_hsc)<-"group"
DefaultAssay(atacs_hsc)<-"ATAC"

cov_junb<-CoveragePlot(
  object = atacs_hsc,
  region = "chr19-12791128-12792822",
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_junb <- AnnotationPlot(
  object = atacs_hsc,
  region = "chr19-12791128-12792822"
)

peak_junb <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = "chr19-12791128-12792822"
)

tile_junb <- TilePlot(
  object = atacs_hsc,
  order.by = "random",
  region = "chr19-12791128-12792822",
  tile.cells = 100
)
tile_junb<-tile_junb+scale_fill_gradient(low = "white",high = "darkred",breaks=c(0,2,4,6))
cbps_hsc<-subset(readRDS("outputs/06-integr_singlecell_cbps/cbps_sct_light.rds"),lineage_hmap=="HSC")
Idents(cbps_hsc)<-"group"
expr_junb <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = "JUNB",
  assay = "SCT"
)

CombineTracks(
  plotlist = list(cov_junb, tile_junb, peak_junb, gene_junb),
  expression.plot = expr_junb,
  heights = c(10, 6, 1, 3),
  widths = c(9, 2)
)

#add MethInfo
res_meth<-fread("outputs/14-DMCs_atac_integr/res_cpgs_hg38.cs.gz")
res_meth_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
#chr19-12791128-12792822
start.pos <- 12791128
  end.pos <- 12792822
  chromosome <- "chr19"
  
res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
res_meth_reg[,start:=pos][,end:=pos+1]
p<-ggplot(data = res_meth_reg) + geom_segment(aes(x = start, y = 0, 
      xend = end, yend = logFC,col=-log10(P.Value)), size = 2, data = res_meth_reg)+
  scale_color_gradient(low = "white",high = "black")

meth_junb<-p+ theme_classic() + ylab(label = "Methylation change") + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))

CombineTracks(
  plotlist = list(cov_junb, meth_junb, peak_junb, gene_junb),
  expression.plot = expr_junb,
  heights = c(10, 6, 1, 3),
  widths = c(9, 2)
)

#addTFMotifInfo
start<-function(x)as.numeric(strsplit(x,"-")[[1]][2])
end<-function(x)as.numeric(strsplit(x,"-")[[1]][3])
seqi<-function(x)strsplit(x,"-")[[1]][1]
peak<-"chr19-12791128-12792822"
start.pos <- start(peak)
end.pos <- end(peak)
chromosome <- seqi(peak)

egr1_ranges<-atacs_hsc@assays$lin_peaks@motifs@positions[[GetMotifIDs(atacs_hsc,motif.names = "EGR1")]]


egr1_dt <- data.table(as.data.frame(egr1_ranges))
egr1_dt_reg<-egr1_dt[seqnames==chromosome&start>start.pos&end<end.pos]

p<-ggplot(data = egr1_dt_reg) + geom_segment(aes(x = start, y = 0, 
      xend = end, yend = 0),col="black", size = 2, data = egr1_dt_reg)

egr1_junb<-p+ theme_classic() + ylab(label = "EGR1") + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))


tf<-"KLF4"
tf_ranges<-atacs_hsc@assays$lin_peaks@motifs@positions[[GetMotifIDs(atacs_hsc,motif.names = tf)]]


tf_dt <- data.table(as.data.frame(tf_ranges))
tf_dt_reg<-tf_dt[seqnames==chromosome&start>start.pos&end<end.pos]

p<-ggplot(data = tf_dt_reg) + geom_segment(aes(x = start, y = 0, 
      xend = end, yend = 0),col="black", size = 2, data = tf_dt_reg)

tf_junb<-p+ theme_classic() + ylab(label = tf) + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
    xlab(label = paste0(chromosome, " position (bp)")) + 
    xlim(c(start.pos, end.pos))

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = peak,
                      motif.names= c("EGR1","KLF4"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_junb, meth_junb, peak_junb,tfs_plot, gene_junb),
  expression.plot = expr_junb,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,"6F-genome_track_JUNB.pdf"),plot = plots,height = 7)



#with all others
genes_of_int<-c("SOCS3","GADD45B", "LMNA","KLF2")

#SOCS3
gene<-"SOCS3"
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peaks_hsc_anno[gene_name==gene]
region<-"chr17-78352543-78361026"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

meth_plot<-MethChangePlot(res_meth,region = region)

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("EGR1","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,".pdf")),plot = plots,height = 7)

#GADD45B
gene<-"GADD45B"
peaks_da_anno[gene_name==gene]
peaks_hsc_anno[gene_name==gene]
region<-"chr19-2487779-2489556"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

meth_plot<-MethChangePlot(res_meth,region = region)

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("KLF4","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)
plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,".pdf")),plot = plots,height = 7)

#LMNA
gene<-"LMNA"
peaks_da_anno[gene_name==gene]
peaks_hsc_anno[gene_name==gene]
region<-"chr1-156113889-156124032"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

meth_plot<-MethChangePlot(res_meth,region = region)

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("KLF4","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)
plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,".pdf")),plot = plots,height = 7)


#KLF2
gene<-"KLF2"
peaks_da_anno[gene_name==gene]
peaks_hsc_anno[gene_name==gene]
region<-"chr19-16319352-16330749"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

meth_plot<-MethChangePlot(res_meth,region = region)

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("KLF4","KLF2","JUNB"),
                      assay = "lin_peaks",pad=20,
                      size=4,alpha = 0.6)
plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,".pdf")),plot = plots,height = 7)


#KLF13
gene<-"KLF13"
peaks_da_anno[gene_name==gene]
peaks_hsc_anno[gene_name==gene]
region<-"chr15-31325484-31332622"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

meth_plot<-MethChangePlot(res_meth,region = region)

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)
plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,".pdf")),plot = plots,height = 7)


#JUNB, KLF2, SOCS3 with same mehthchange pval col
#JUNB
res_meth_anno[gene=="JUNB"][order(pval)]
gene<-"JUNB"
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peaks_hsc_anno[gene_name==gene]
region<-"chr19-12790000282794822"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[-log10(P.Value)>5]
res_meth_reg[order(P.Value)]
meth_plot<-MethChangePlot(res_meth,region = region,limits = c(0,4))

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("EGR1","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)


#SOCS3
gene<-"SOCS3"
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peaks_hsc_anno[gene_name==gene]
region<-"chr17-78352543-78361026"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)
res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[-log10(P.Value)>5]
res_meth_reg[order(P.Value)]
meth_plot<-MethChangePlot(res_meth,region = region,limits = c(0,4))

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("EGR1","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)

#KLF2
gene<-"KLF2"
peaks_da_anno[gene_name==gene]
peaks_hsc_anno[gene_name==gene]
region<-"chr19-16319352-16330749"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)
res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[-log10(P.Value)>5]
res_meth_reg[order(P.Value)]
meth_plot<-MethChangePlot(res_meth,region = region,limits = c(0,4))

tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("KLF4","KLF2","JUNB"),
                      assay = "lin_peaks",pad=20,
                      size=4,alpha = 0.6)
plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)

#pval<0.05
#JUNB
res_meth_anno[gene=="JUNB"][order(pval)]
gene<-"JUNB"
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
peaks_hsc_anno[gene_name==gene]
region<-"chr19-12790128-12793822"

cov_plot<-CoveragePlot(
  object = atacs_hsc,
  region = region,
  group.by = "group",
  annotation = F,
  peaks = F
)

gene_plot <- AnnotationPlot(
  object = atacs_hsc,
  region = region
)

peak_plot <- PeakPlot(
  object = atacs_hsc,assay = "lin_peaks",
  region = region
)

expr_plot <- ExpressionPlot(
  object = subset(cbps_hsc,hto==T),
  features = gene,
  assay = "SCT"
)

res_meth_reg<-MethChangeReg(res_meth,region)
res_meth_reg[,pval:=sapply(P.Value,function(x)ifelse(x<0.001,"***",ifelse(x<0.01,"**",ifelse(x<0.05,"*","ns"))))]
res_meth_reg[,pval:=factor(pval,levels = c("ns","*","**","***"))]
meth_plot<-ggplot(data = res_meth_reg) + geom_segment(aes(x = start, y = 0, 
        xend = end, yend = logFC,col=pval), size = 2, data = res_meth_reg)+
    scale_color_manual(values = c("grey","yellow","orange","red"))
  
meth_plot<-meth_plot+ theme_classic() + ylab(label = "Methylation change") + 
    xlab(label = paste0(seqid(region), " position (bp)")) + 
    xlim(c(start(region), end(region)))


tfs_plot<-TFsMotifPlot(atacs_hsc,
                      region = region,
                      motif.names= c("EGR1","KLF2"),
                      assay = "lin_peaks",
                      size=4,alpha = 0.6)

plots<-CombineTracks(
  plotlist = list(cov_plot, meth_plot, peak_plot,tfs_plot, gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1,2,3),
  widths = c(9, 2)
)
plots
ggsave(fp(out,ps("6F-genome_track_",gene,"_v2.pdf")),plot = plots,height = 7)




#7G: pseudotime bias rna velo based
mtdff<-fread("outputs/20-RNA_velocity/pseudo_bias_rna_velo_based_lineages.csv.gz")
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")


mtdff[,lineage_hmap:=factor(lineage_hmap,levels=c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas"))]
ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group,group=sample),outlier.shape = NA)+
  facet_wrap("lineage_hmap")

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))+
  coord_cartesian(ylim = c(-25,100))
ggsave(fp(out,"7G-boxplot_differentiation_bias_cell_level.pdf"))

mtdff[,pval:=t.test(pseudo_bias[group=="lga"],
                           pseudo_bias[group=="ctrl"])$p.value,
        by="lineage_hmap"]
unique(mtdff[,.(lineage_hmap,pval)])

#    lineage_hmap         pval
# 1:     MPP/LMPP 4.843362e-01
# 2:          HSC 9.008333e-11
# 3:  Erythro-Mas 5.701923e-04
# 4:     Lymphoid 4.429514e-04
# 5:      Myeloid 3.395650e-07
# 6:       LT-HSC 4.602892e-01

# ggplot(mtdff[lineage_hmap%in%lins])+
#   geom_boxplot(aes(x=group,y=pseudo_bias,fill=group),outlier.shape = NA)+
#   facet_wrap("lineage_hmap")+
#   coord_cartesian(ylim = c(0,50))

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

mtdffl<-mtdff[lineage_hmap%in%lins]
mtdffl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdffl[,avg_pseudobias:=mean(pseudo_bias),by=.(lineage_hmap,sample)]
mtdffls<-unique(mtdffl,by=c("sample","lineage_hmap"))
ggplot(mtdffls)+geom_boxplot(aes(x=lineage_hmap,y=avg_pseudobias,fill=group))
ggsave(fp(out,"7G-boxplot_differentiation_bias_sample_level.pdf"))
mtdffls[,pval:=t.test(avg_pseudobias[group=="lga"],
                           avg_pseudobias[group=="ctrl"])$p.value,
        by="lineage_hmap"] 
unique(mtdffls[,.(lineage_hmap,pval)])

#    lineage_hmap       pval
# 1:     MPP/LMPP 0.41856808
# 2:          HSC 0.03910725
# 3:  Erythro-Mas 0.91731411
# 4:     Lymphoid 0.75722049
# 5:      Myeloid 0.32957624
# 6:       LT-HSC 0.27182012


#unique(mtdffl[lineage_hmap=="MPP/LMPP"],by=c("sample","lineage_hmap"))

#SUPP####
#distrib predic lineage
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

ggplot(atacs@meta.data)+geom_bar(aes(x=predicted.id))

#can remove badly predict lineage ?
ggplot(atacs@meta.data)+geom_bar(aes(x=predicted.id))



