
source("scripts/utils/new_utils.R")
out0<-"outputs/figures_epi_response"
null<-0
#figure 1 : Epigenetic Memory ####
out<-fp(out0,"figure1")
dir.create(out,recursive = T)

#supp : PC plot et cohort
#PC plot
library(pheatmap)
source("scripts/utils/methyl_utils.R")
pval_cov<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/pval_correl_covar_pcs.csv")

pval_covf<-pval_cov[!covar%in%c("groupbatch_complexity_fac","seq.depth","library_complexity","group_complexity","groupbatch_complexity")]
pval_covf[covar=="group_complexity_fac",covar:="library_complexity"]
pca<-readRDS("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/pca_lgactrl.rds")

pval_mat<-as.matrix(data.frame(pval_covf,row.names = 1))

pct.varPCs<-pctPC(pca,1:10)*100
break_cols<-c(1,0.1,0.05,0.01,0.001,0.0001,0.00001)

pdf(fp(out,"supp10-pcs_covars_plot.pdf"))
pheatmap(-log10(pval_mat),cluster_cols = F,
           labels_col= paste0("PC",1:10,"(",round(pct.varPCs[as.character(1:10)]),"%)"),
           display_numbers = T,
           color = colorRampPalette(c("white", "red"))(length(break_cols)+1),
         breaks =sort(-log10(break_cols)))
dev.off()

#cohort
res_coh<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma_cohorts.tsv.gz")
res_coh<-res_coh[compa=="C.L"]
table(res_coh[adj.P.Val<=0.1,.(batch)])
ggplot(res_coh)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>25),size=0.1)+
  facet_wrap("batch")+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()
ggsave(fp(out,"supp10-volcanos_limma_cohorts_C.L.pdf"),width = 12,height = 5)



#1A : volcanoplot 
res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
res[,meth.change:=logFC][,p_val_adj:=adj.P.Val][,p_val:=P.Value]

res[p_val_adj<0.05&abs(meth.change)>25] #43 DMCs
res[p_val_adj<0.1&abs(meth.change)>25] #1255 DMCs
res[p_val_adj<0.1&meth.change>25] #1250 DMCs

res[p_val<0.001&abs(meth.change)>25] #4815 DMCs
res[p_val<0.001&meth.change>25] #4787 
4815-4787 
#28 hypo
cpgs<-fread("ref/all_HELP_tagging_cpgs_hg19_pos.csv")
res<-merge(res,cpgs)
fwrite(res[p_val<0.001&abs(meth.change)>25],fp(out,"supptabl2-res_limma_pval0.001_meth.change25.csv"))

ggplot(res)+
  geom_point(aes(x=meth.change,y=-log10(p_val),col=p_val<0.001&abs(meth.change)>25),size=0.1)+
  scale_color_manual(values = c("grey","red"))+scale_x_continuous(limits = c(-100,100))+theme_minimal()
ggsave(fp(out,"1A-volcano_plot_hypermet_LGA.pdf"))

#table : genescore, gene
resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
fwrite(resg[,.(gene,gene_score_add,n.cpg.gene,n.cpg.sig.gene,chr,pos,pval,padj)],fp(out,"table-genes_methylation_genescore.csv"))
summary(resg[n.cpg.sig.gene>3]$gene_score_add)

#suppfig: valid genescore
res_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")
resg<-merge(resg,unique(res_anno[order(gene,pval)][,.(gene,chromatin_feature,ensembl_reg_score,in_eQTR)],by="gene"))
res_mod<-summary(lm(gene_score_add~n.cpg.gene+n.cpg.sig.gene+pval+meth.change+chromatin_feature+ensembl_reg_score+in_eQTR+abs(tss_dist),data = resg)) #best gene_score_add than gene_score
res_mod_dt<-data.table(covariate=rownames(res_mod$coefficient),pval=res_mod$coefficients[,4])
ggplot(res_mod_dt[covariate!='(Intercept)'])+geom_col(aes(x=covariate,y=-log10(pval)))+scale_y_log10()
ggsave(fp(out,"suppfig1-covariates_influence_on_genescore.pdf"))

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
ggsave(fp(out,"suppfig1-gene_score_pred_expression_change_compared_9classical_metrics.pdf"))

meth_scores_de[,pval:=wilcox.test(score[padj<0.05],score[padj>=0.05])$p.value,by="meth_metric"]
unique(meth_scores_de[,.(meth_metric,pval)],by="meth_metric")

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

ggplot(unique(meth_scores_de,by="meth_metric"))+geom_col(aes(y=-log10(pval),x=meth_metric))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


#1B : pathway GSEA 
library(enrichplot)
res_go<-readRDS("outputs/03-pathway_analysis/res_gsea_go.rds")
max(data.table(as.data.frame(res_go))$p.adjust) #699 go term with padj < 0,001

pdf(fp(out,"1B-emapplot_gsea_go.pdf"),width = 14,height = 8)
emapplot(pairwise_termsim(res_go,showCategory = 40),showCategory = 40)
dev.off()

gsea_go<-fread("outputs/03-pathway_analysis/res_gsea_go.csv")

pdf(fp(out,"1B-dotplot_gsea_go_x_avg_genescore.pdf"),width = 8,height = 10)
dotplot(res_go,x=gsea_go[order(p.adjust)]$gene_score.avg[1:40],showCategory=40)
dev.off()


#1C : TF motif enrichment
#with own background no norm, 1000 perm
res<-fread("outputs/03B-motif_analysis/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[padj<0.05]
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]
ggplot(res[padj<=0.05&  n_dmc_with_motif>30][order(pval)])+geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_x_discrete(limits=res[padj<=0.05&  n_dmc_with_motif>30][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust=1))

ggsave(fp(out,"1C-motif_enrichment_homer_own_background.pdf"))

res_perm<-fread("outputs/03B-motif_analysis/res_known_motif_all_perm_bgrandom.csv")
res_perm[,motif:=str_remove(motif,"/Homer")]

res[,permut:=0]
res_perm_merge<-merge(res,res_perm,all=T)
res_perm_merge[,p.perm:=sum(pct_dmc_with_motif[permut==0]<=pct_dmc_with_motif[permut!=0])/sum(permut!=0,na.rm = T),by=.(motif)]

res_perm_merge[padj<=0.05&n_dmc_with_motif>30&p.perm<0.05&permut==0]
res_perm_<-res_perm_merge[permut==0]

res_perm_[padj<=0.05] #26
res_perm_[padj<=0.05&n_dmc_with_motif>30] #26
 #26
res_perm_[padj<=0.05&n_dmc_with_motif>30&p.perm<0.01] #26

ggplot(res_perm_[padj<=0.05&pct_dmc_with_motif>1&p.perm<0.01])+
  geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_color_gradient(high ="brown1" ,low = scales::muted("red"))+
  scale_x_discrete(limits=res_perm_[padj<=0.05&pct_dmc_with_motif>1&p.perm<0.01][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 80, hjust=1))

ggsave(fp(out,"1C-motif_enrichment_homer_own_background_no_norm_padj0.05_pctDMC1_p.perm0.01.pdf"),height = 6.5)

fwrite(res_perm_,fp(out,"1C-res_homer_own_background_no_norm_padj0.05_pctDMC1_p.perm0.01.csv"))

res_perm_[p.perm<0.01] #32

#supp : cohort validation 

res_coh<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma_cohorts.tsv.gz")

res_coh[P.Value<0.001&abs(logFC)>25&batch==1&compa=="C.L"] #2220

nrow(res_coh[P.Value<0.001&logFC>25&batch==1&compa=="C.L"])/
  nrow(res_coh[P.Value<0.001&abs(logFC)>25&batch==1&compa=="C.L"] )#94%

res_coh[P.Value<0.001&logFC>25&batch==2&compa=="C.L"] #2250
nrow(res_coh[P.Value<0.001&logFC>25&batch==2&compa=="C.L"])/
  nrow(res_coh[P.Value<0.001&abs(logFC)>25&batch==2&compa=="C.L"] )#96%
length(intersect(res_coh[P.Value<0.001&abs(logFC)>25&batch==2&compa=="C.L"]$cpg_id,
                 res_coh[P.Value<0.001&abs(logFC)>25&batch==1&compa=="C.L"]$cpg_id))


#supp : gwas 
res_gwas<-readRDS("outputs/03-pathway_analysis/res_gsea_gwas.rds")
max(data.table(as.data.frame(res_gwas))$p.adjust) #86 go term with padj < 0,01

pdf(fp(out,"1B-emapplot_gsea_gwas.pdf"),width = 14,height = 8)
emapplot(pairwise_termsim(res_gwas,showCategory = 86),showCategory = 86)
dev.off()

pdf(fp(out,"1B-dotplot_gsea_gwas.pdf"),width = 8,height = 12)
dotplot(res_gwas,showCategory = 86)
dev.off()


# supp : kegg and go together [to update]
df_kegg_go<-merge(res_kegg,res_go,all=T) #merge in a dataframe
head(df_kegg_go)
#need transfoom in enrichment results :
# renv::install("bioc::ComplexHeatmap") #need this dependensies 
# renv::install("jmw86069/multienrichjam")
library(multienrichjam)
res_kegg_go<-enrichDF2enrichResult(df_kegg_go,keyColname = "ID",pvalueColname = "p.adjust",geneColname = "leading_edge")
data.table(as.data.frame(res_kegg_go))[!str_detect(ID,"GO")]

emapplot(pairwise_termsim(res_kegg_go,showCategory = 80),showCategory = 80)#not well informative


#supp : tf motif with auto background [to update]
res<-fread("outputs/03B-motif_analysis/with_auto_background/knownResults.txt",
           select = c(1,2,3,5,6,7,8,9),
           col.names = c("motif","consensus","pval","padj","n_dmc_with_motif","pct_dmc_with_motif","n_background_with_motif","pct_background_with_motif"))
res[,pct_dmc_with_motif:=as.numeric(str_remove(pct_dmc_with_motif,"%"))]
res[,pct_background_with_motif:=as.numeric(str_remove(pct_background_with_motif,"%"))]
res[,motif:=str_remove(motif,"/Homer")]
res[,fold.enrichment:=pct_dmc_with_motif/pct_background_with_motif]
ggplot(res[padj<=0.2&  n_dmc_with_motif>30][order(pval)])+geom_point(aes(x=motif,col=-log10(pval),size=fold.enrichment,y=n_dmc_with_motif))+
  scale_x_discrete(limits=res[padj<=0.2&  n_dmc_with_motif>30][order(pval)]$motif)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 70, hjust=1))

ggsave("outputs/figures_epi_response/figure1/1C-motif_enrichment_homer_auto_background.pdf",width = 16)


#supp : Network [to optimize]
networks<-readRDS("outputs/04-comethylation_analysis/megena_spearman_pfn_mca_outputs.rds")
networks

networks.sum<-readRDS("outputs/04-comethylation_analysis/summary_megena_outputs.rds")
head(networks.sum$module.table,100)

library(ggplot2)
library(ggraph)

g <- readRDS("outputs/04-comethylation_analysis/graph.rds")
pnet.obj <- plot_module(output.summary = networks.sum,PFN = g,subset.module = "c1_29",
	layout = "kamada.kawai",label.hubs.only = FALSE,
	gene.set = NULL,color.code =  "grey",
	output.plot = FALSE,out.dir = "modulePlot",col.names = c("magenta","green","cyan"),label.scaleFactor = 20,
	hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)
print(pnet.obj[[1]]) #not informative so next



#figure 2 : Hematomap ####
source("scripts/utils/new_utils.R")
library(Seurat)
hmap<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/hematomap_ctrls_sans_stress/hematomap_ctrls_sans_stress.rds")
out<-fp(out0,"figure2")
dir.create(out)

#table : number of cells / ech
#basal
#duplicates
#hto

#2A : UMAP 
#rm cluster 18
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
hmap[["lineage"]]<-Idents(hmap)
DimPlot(hmap,group.by = c("lineage"),label = T)
ggsave("outputs/figures_epi_response/figure2/2A-hematomap.pdf")


#2B: key genes/lin :
DefaultAssay(hmap)<-"SCT"
key_genes_lin<-c("ID1","ID2","DUSP2", #LT-HSC
           "EGR1","AVP", #HSC
           "MLLT3","CDK6", #MPP
           "SELL","CD99", #LMPP
           "LTB", #CLP
         "VPREB1","IGLL1", # proB
           "IGHM","CD37", #B cell
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "MKI67","TOP2A", #ErP-cycle
           "PLEK","HBD", #Mk/Er
           "MPO","CEBPA","CTSG","AZU1", #GMP
         "CST3","CD83") #DC
DotPlot(hmap,features = key_genes_lin,group.by = "lineage")
ggsave("outputs/figures_epi_response/figure2/2B-key_genes_by_lin.pdf",width = 19)

#supp : umap key genes
ps<-FeaturePlot(hmap,features = c("ID1","EGR1","AVP","GATA1","HDC","MPO","CST3","LTB","VPREB1"),combine = F)
ps<-lapply(ps, function(x)x+NoAxes()+NoLegend())
wrap_plots(ps)
ggsave(fp(out,"2supp-umap_key_genes_by_lin.pdf"),width = 12,height = 12)

#supp : key genes /ct
#reorder 
Idents(hmap)<-"cell_type"
levels(hmap)<-c("LT-HSC",
                "HSC-1",
                "HSC-2",
                "HSC-3",
                "HSC-4",
                "MPP",
                "LMPP",
                "CLP",
                "proB",
                "B cell",
                "T cell",
               "MPP-Ery",
               "EMP",
               "EMP-cycle",
               "ErP-cycle",
               "Mk/Er",
               "GMP",
               "GMP-cycle",
               "DC")
hmap[["cell_type"]]<-Idents(hmap)

#umap
DimPlot(hmap,group.by = "cell_type",label = T)
ggsave("outputs/figures_epi_response/figure2/suppfig2-umap_celltype.pdf")

#key genes
key_genes_ct<-c("ID1","ID2","DUSP2", #LT-HSC
           "EGR1","AVP", #HSC-1
         "CXCL8","ZFP36", #HSC-2
         "IRF1","STAT1", #HSC-3/MkP
          "KLF2","TSC22D3", #HSC-4/proT
           "MLLT3","CDK6", #MPP
           "SELL","CD99", #LMPP
           "LTB", #CLP
         "VPREB1","IGLL1", # proB
           "IGHM","CD37", #B cell
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "BIRC5","MKI67","TOP2A", #ErP-cycle
           "HDC", #Mast
           "PLEK","HBD", #Mk/Er
           "MPO","AZU1", #GMP
         "CST3","CD83") #DC

DotPlot(hmap,features = key_genes_ct,group.by = "cell_type")
ggsave("outputs/figures_epi_response/figure2/supp-key_genes_by_celltype.pdf",width = 23)

#supp : LGA vs CTRL - Basal condition
#same celltype distrib

mtd<-fread("outputs/07-LGA_vs_Ctrl_Basal/metadata_basal.csv")

ggplot(unique(mtd,by=c("sample","lineage_hmap")))+
  geom_boxplot(aes(x=lineage_hmap,y=pct.lin,fill=group))
ggsave(fp(out,"supp-distribution_lineage_control_lga_basal.pdf"))

#no DEGs
res_lin<-fread("outputs/07-LGA_vs_Ctrl_Basal/res_pseudobulkDESeq2_by_lineage.csv.gz")
res_lin$lineage<-factor(res_lin$lineage,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))

#volcano by lin
genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")


ggplot(res_lin[lineage%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],
       aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.05&
                                        abs(log2FoldChange)>0.5,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("lineage")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"supp-pseudo_bulk_deseq2_by_lineage_lga_vs_ctrl_basal.pdf"))


#figure 3 : LGA stimulation response   ####
source("scripts/utils/new_utils.R")
out<-fp(out0,"figure3")
dir.create(out)

#supp : activation signature 
res_hto_dup<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")

hto_signature<-res_hto_dup[padj<0.05&abs(log2FoldChange)>0.5]$gene
length(hto_signature) #1518

fwrite(res_hto_dup[padj<0.05&abs(log2FoldChange)>0.5],fp(out,"table-stimulation_signature_genes.csv"))

res_hto_dup[padj<0.05&log2FoldChange>0.5] 
res_hto_dup[padj<0.05&abs(log2FoldChange)>0.5&gene%in%c("SOCS3","HES1","JUN","EGR1")]
genes_of_interest<-c("SESN2","SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

ggplot(res_hto_dup,aes(x=log2FoldChange,y=-log10(padj),col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),

                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(fp(out,"supp-volcano_antibody_activation_signature_3_ctrl_same_samples.pdf"))


#supp : pathway activation 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#go term
res_go_up<-fread("outputs/08-HTO_signature/res_hto_signature_go_bp_up.csv")
res_go_dn<-fread("outputs/08-HTO_signature/res_hto_signature_go_bp_dn.csv")
res_go<-rbind(res_go_up[,query:="up"],res_go_dn[,query:="dn"])
res_kegg_up<-fread("outputs/08-HTO_signature/res_hto_signature_kegg_up.csv")
res_kegg_dn<-fread("outputs/08-HTO_signature/res_hto_signature_kegg_dn.csv")
res_kegg<-rbind(res_kegg_up[,query:="up"],res_kegg_dn[,query:="dn"])
res<-rbind(res_go[,source:="GO BP"],res_kegg[,source:="KEGG"])
res[,n.overlap:=as.numeric(Count)]
res[,n.gene_set:=as.numeric(str_extract(BgRatio,"[0-9]+"))]

res[,pct.overlap:=n.overlap/n.gene_set]

terms_of_interest<-c("cell activation" ,
                    "regulation of cell activation",
                     "cellular response to external stimulus" ,
                      "regulation of response to external stimulus" ,
                     "cellular response to stress",
                    "regulation of cellular response to stress" ,
                     "cellular response to growth factor stimulus",
                     "regulation of cellular response to growth factor stimulus" ,
                      "regulation of cell differentiation",
                     "stem cell differentiation",
                     "regulation of hemopoiesis" ,
                     "hematopoietic progenitor cell differentiation"   ,
                    "cell population proliferation",
                     "regulation of cell population proliferation",
                     "regulation of cell cycle" ,
                    "MAPK signaling pathway",
                    "TNF signaling pathway",
                    "FoxO signaling pathway",
                    "NF-kappa B signaling pathway",
                    "Hippo signaling pathway",
                    "TGF−beta signaling pathway")
terms_of_interest<-c("cellular response to cytokine stimulus"   ,
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
  facet_grid(source~query)+
  scale_color_manual(values = c("grey",'red'))+theme_minimal()
ggsave(fp(out,"supp-points_hto_up_dn_kegg_go.pdf"))

fwrite(res,fp(out,"res_enrich_hto_up_dn__kegg_gobp.csv"))

res_sig_up_dt<-data.table(as.data.frame(res_sig_up))
res_sig_up_dt[p.adjust<0.05]

pdf(fp(out,"supp4C-emapplot_pathways_activation_signature_up_padj0.05.pdf"),width = 12)
emapplot(pairwise_termsim(res_sig_up),showCategory = 49)
dev.off()

pdf(fp(out,"supp4C-emapplot_pathways_activation_signature_up_padj0.05_top25.pdf"),width = 12)
emapplot(pairwise_termsim(res_sig_up),showCategory = 25)
dev.off()

pdf(fp(out,"supp4C-dotplot_pathways_activation_signature_up_padj0.05.pdf"),width = 8,height =12 )
dotplot(res_sig_up,showCategory = 39)
dev.off()

res_or_sig_up<-fread("outputs/08-HTO_signature/res_hto_signature_kegg_up.csv")

fwrite(res_or_sig_up[p.adjust<0.1],fp(out,"supp-table_pathway_enrichment_hto_signature_up_padj0.1.csv"))


#supp : signature by lineage
#nDEGs
res_hto_lin<-fread("outputs/08-HTO_signature/by_lineage/res_pseudobulk_DESeq2_3replicates.csv.gz")
table(res_hto_lin[padj<0.05&abs(log2FoldChange)>0.5]$lineage)#1184
    # B cell          DC Erythro-Mas         HSC 
    #      51         109         161        1264 
    #  LT-HSC    Lymphoid    MPP/LMPP     Myeloid 
    #      14          61         341         188 

res_hto_lin[,pct.degs:=sum(padj<0.05&abs(log2FoldChange)>0.5&gene%in%hto_signature,na.rm=T)/length(intersect(hto_signature,gene)),by="lineage"]

res_hto_lin[,n.genes:=length(intersect(hto_signature,gene)),by="lineage"]
res_hto_lin[,n.degs:=sum(padj<0.05&abs(log2FoldChange)>0.5&gene%in%hto_signature,na.rm=T),by="lineage"]

res_hto_linf<-unique(res_hto_lin[,.(lineage,n.genes,n.degs,pct.degs)])

?chisq.test
chisq.test(res_hto_linf$n.degs,p=res_hto_linf$n.genes, rescale.p = T) 
chisq.test(res_hto_linf$n.degs,p=res_hto_linf$n.genes/sum(res_hto_linf$n.genes)) # p-value < 2.2e-16


res_hto_lin[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5&gene%in%hto_signature] #796/1518

ggplot(res_hto_lin[padj<0.05&abs(log2FoldChange)>0.5])+geom_bar(aes(x=lineage,fill=lineage))+theme_minimal()
ggsave(fp(out,"supp-n_degs_by_lineage_3ctrl_stimulation.pdf"))

ggplot(res_hto_lin[padj<0.05&abs(log2FoldChange)>0.5&gene%in%hto_signature])+geom_bar(aes(x=lineage,fill=lineage))+theme_minimal()
ggsave(fp(out,"supp-n_degs_of_hto_stim_by_lineage_3ctrl_stimulation.pdf"))

#heatmap 
library(pheatmap)
res_hto_all<-rbind(res_hto_lin,fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")[,lineage:="all_cbps"],fill=T)

res_hto_mat<-dcast(res_hto_all[gene%in%hto_signature&lineage%in%c("all_cbps","LT-HSC","HSC","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas")],gene~lineage,value.var ="log2FoldChange")
head(res_hto_mat)
res_hto_mat<-as.matrix(data.frame(res_hto_mat,row.names = "gene"))
head(res_hto_mat)
res_hto_matf<-res_hto_mat[!rowSums(is.na(res_hto_mat))==1,]#removing degs just in all cbps
res_hto_matf[is.na(res_hto_matf)]<-0


pdf(fp(out,"supp-heatmap_activation_signature_by_lineage_3ctrls_no_scaling.pdf"))
break_cols<-c(-30:0/10,1:30*1.25/10)
pheatmap(res_hto_matf,show_rownames = F,na_col = "grey",
         color = colorRampPalette(c("darkblue","white", "red"))(length(break_cols)+1),
         breaks =break_cols )
dev.off()

res_hto_mat_scaled<-t(scale(t(res_hto_matf),center = F))
pdf(fp(out,"supp-heatmap_activation_signature_by_lineage_3ctrls_scaled.pdf"))
pheatmap(res_hto_mat_scaled[,],cluster_cols = T,show_rownames = F,na_col = "grey")
dev.off()

#supp - all compa by lineages 
res_lin_all<-Reduce(rbind,list(fread("outputs/08-HTO_signature/pseudobulk_DESeq2_ctrl_hto/res_pseudobulkDESeq2_by_lineage.csv.gz")[,compa:="ctrl_hto"],
                           fread("outputs/08-HTO_signature/pseudobulk_DESeq2_lga_hto/res_pseudobulkDESeq2_by_lineage.csv.gz")[,compa:="lga_hto"],
                           fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[,compa:="lga_vs_ctrl_hto"]
))
res_lin_all<-res_lin_all[lineage%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")]

res_lin_all$lineage<-factor(res_lin_all$lineage,levels = c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid"))
#barplot degs by compa
ggplot(res_lin_all[padj<0.05&abs(log2FoldChange)>0.5])+geom_bar(aes(x=lineage,fill=lineage))+
  facet_wrap("compa")+scale_x_discrete(guide = guide_axis(angle = 45))

ggsave(fp(out,"suppFig5-barplot_ndegs_by_lineage_and_compa.pdf"))

#barplot degs_sign by compa
ggplot(res_lin_all[padj<0.05&abs(log2FoldChange)>0.5&gene %in% hto_signature])+geom_bar(aes(x=lineage,fill=lineage))+
  facet_wrap("compa")+scale_x_discrete(guide = guide_axis(angle = 45))
ggsave(fp(out,"suppFig5-barplot_ndegs_signature_by_lineage_and_compa.pdf"))

#volcanoplot degs by compa
res_lin_all[padj<0.05&abs(log2FoldChange)>0.5,deg_sig:="deg"]
res_lin_all[padj<0.05&abs(log2FoldChange)>0.5&gene %in% hto_signature,deg_sig:="deg_hto"]
res_lin_all[is.na(deg_sig),deg_sig:="other_gene"]
#ctrl
ggplot(res_lin_all[compa=="ctrl_hto"],aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point(size=0.6)+
  facet_wrap("lineage")+
  scale_color_manual(values = c("red","orange","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"suppfig5A-volcano_degs_hto_sign_ctrl_stim_vs_ctrl_steady_state.pdf"))

#lga
ggplot(res_lin_all[compa=="lga_hto"],aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point(size=0.6)+
  facet_wrap("lineage")+
  scale_color_manual(values = c("red","orange","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"suppfig5B-volcano_degs_hto_sign_lga_stim_vs_lga_steady_state.pdf"))

#lga ctrl
ggplot(res_lin_all[compa=="lga_vs_ctrl_hto"],aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point(size=0.6)+
  facet_wrap("lineage")+
  scale_color_manual(values = c("red","orange","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"suppfig5C-volcano_degs_hto_sign_lga_stim_vs_ctrl_stim.pdf"))


#3A : stimulation sign across group
library(pheatmap)

res_hto<-Reduce(rbind,list(fread("outputs/08-HTO_signature/pseudobulk_DESeq2_ctrl_hto/res_all_cbps_de_analysis.csv")[,compa:="ctrl_hto"],
                           fread("outputs/08-HTO_signature/pseudobulk_DESeq2_lga_hto/res_all_cbps_de_analysis.csv")[,compa:="lga_hto"],
                           fread("outputs/08-HTO_signature/pseudobulk_DESeq2_lga_vs_ctrl_hto/res_all_cbps_de_analysis.csv")[,compa:="lga_vs_ctrl_hto"]
))

res_hto_com<-res_hto[compa%in%c('ctrl_hto','lga_hto')]
res_hto_com<-res_hto_com[gene%in%gene[duplicated(gene)]]
res_hto_com<-res_hto_com[order(gene)]
wilcox.test(rank(res_hto[compa=="ctrl_hto"]$log2FoldChange),rank(res_hto[compa=="lga_hto"]$log2FoldChange))


cor(res_hto_com[compa=="ctrl_hto"]$log2FoldChange,res_hto_com[compa=="lga_hto"]$log2FoldChange,method = "spearman")
cor.test(res_hto_com[compa=="ctrl_hto"]$log2FoldChange,res_hto_com[compa=="lga_hto"]$log2FoldChange,method="spearman")
# p-value < 2.2e-16, rho = 0.70



signat_inter<-intersect(res_hto[compa=="ctrl_hto"&padj<0.05&log2FoldChange>0.5&gene%in%hto_signature]$gene,res_hto[compa=="lga_hto"&padj<0.05&log2FoldChange>0.5&gene%in%hto_signature]$gene)
length(signat_inter) #328 / 372 DEGS LGA, /686 DEGs Ctrl


res_hto_mat<-dcast(res_hto[gene%in%hto_signature],gene~compa,value.var ="log2FoldChange")
head(res_hto_mat)
res_hto_mat<-as.matrix(data.frame(res_hto_mat,row.names = "gene"))
head(res_hto_mat)

pdf(fp(out,"3a-heatmap_antibody_activation_signature_by_group_all_cbps.pdf"),width = 3,height = 6)

pheatmap(res_hto_mat,cluster_cols = F,show_rownames = F,
         color = colorRampPalette(c("darkblue","white", "red"))(42),
         breaks = breakRes<-c(-20:0/10,1:20*2/10))
dev.off()


#same with HSC
res_hto_hsc<-Reduce(rbind,list(fread("outputs/08-HTO_signature/pseudobulk_DESeq2_ctrl_hto/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,compa:="ctrl_hto"],
                           fread("outputs/08-HTO_signature/pseudobulk_DESeq2_lga_hto/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,compa:="lga_hto"],
                           fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][,compa:="lga_vs_ctrl_hto"]
))
res_hto_com<-res_hto_hsc[compa%in%c('ctrl_hto','lga_hto')]
res_hto_com<-res_hto_com[gene%in%gene[duplicated(gene)]]
res_hto_com<-res_hto_com[order(gene)]
wilcox.test(rank(res_hto[compa=="ctrl_hto"]$log2FoldChange),rank(res_hto[compa=="lga_hto"]$log2FoldChange))


cor(res_hto_com[compa=="ctrl_hto"]$log2FoldChange,res_hto_com[compa=="lga_hto"]$log2FoldChange,method = "spearman")
cor.test(res_hto_com[compa=="ctrl_hto"]$log2FoldChange,res_hto_com[compa=="lga_hto"]$log2FoldChange,method="spearman")
# p-value < 2.2e-16, rho = 0.65


signat_inter<-intersect(res_hto_hsc[compa=="ctrl_hto"&padj<0.05&log2FoldChange>0.5&gene%in%hto_signature]$gene,res_hto_hsc[compa=="lga_hto"&padj<0.05&log2FoldChange>0.5&gene%in%hto_signature]$gene)
length(signat_inter) #154 / 182 DEGS LGA, /490 DEGs Ctrl


res_hto_mat<-dcast(res_hto_hsc[gene%in%hto_signature&!is.na(log2FoldChange)],gene~compa,value.var ="log2FoldChange")
res_hto_mat<-as.matrix(data.frame(res_hto_mat,row.names = "gene"))
res_hto_mat<-res_hto_mat[rowSums(is.na(res_hto_mat))==0,]
head(res_hto_mat)

pdf(fp(out,"3a-heatmap_antibody_activation_signature_by_group_hsc.pdf"),width = 3,height = 6)

pheatmap(res_hto_mat,cluster_cols = F,show_rownames = F,
         color = colorRampPalette(c("darkblue","white", "red"))(42),
         breaks = breakRes<-c(-20:0/10,1:20*2/10))

dev.off()



#3B : LGA vs CTRL - activated condition 
res_lin<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")
res_lin[padj<0.05&abs(log2FoldChange)>0.5,deg_sig:="deg"]

res_lin[padj<0.05&abs(log2FoldChange)>0.5&gene %in% hto_signature,deg_sig:="deg_hto"]
res_lin[is.na(deg_sig),deg_sig:="other_gene"]

table(res_lin[padj<0.05&abs(log2FoldChange)>0.5]$lineage)
     # B cell          DC Erythro-Mas         HSC    Lymphoid    MPP/LMPP     Myeloid      T cell 
     #     64           5          30         373           7         311          42           1 


res_lin[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5&log2FoldChange>0]#88 upreg
#dont KLF3 qui represse KLF1 pour rentrer en differentiation vers Erythroid (https://journals.asm.org/doi/full/10.1128/MCB.00173-12 )
#et CYTOR, qui promeut proliferation

res_lin[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5&log2FoldChange<0]#285 dnreg

res_hsc<-res_lin[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5]
fwrite(res_hsc,fp(out,"table-res_LGA_vs_ctrl_hsc_activated_padj0.05_logFC0.5.csv"))


ggplot(res_lin[lineage=="HSC"],aes(x=log2FoldChange,y=-log10(padj),col=deg_sig))+
  geom_point()+
   geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("red","blue","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(fp(out,"3b-volcano_LGA_vs_Ctrl_HSC_activation_response_signature_highlight.pdf"))

#supp pathways HSC up and dn  
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

res_hsc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg.csv")
res_hsc[p.adjust<0.1]
fwrite(res_hsc[p.adjust<0.1],fp(out,"supp-pathway_enrichment_LGA_HSC_activated_padj0.1.csv"))

res_hsc_up<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg_up.csv")
res_hsc_up[p.adjust<0.1]
fwrite(res_hsc_up[p.adjust<0.1],fp(out,"supp-pathway_enrichment_LGA_HSC_activated_up_padj0.1.csv"))

res_hsc_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg_dn.csv")
res_hsc_dn[p.adjust<0.1]
fwrite(res_hsc_dn[p.adjust<0.1],fp(out,"supp-table_pathway_enrichment_LGA_HSC_activated_dn_padj0.1.csv"))

res_hsc_dn[,size.pathway:=as.numeric(str_extract(BgRatio,"^[0-9]+"))]
res_hsc_dn[,pct.overlap:=Count/size.pathway]


ggplot(res_hsc_dn[p.adjust<0.1])+geom_point(aes(x=Description,y=pct.overlap,size=Count,col=p.adjust))+
  scale_x_discrete(limits=res_hsc_dn[p.adjust<0.1][order(p.adjust)]$Description)+
  scale_color_continuous(trans="log10")+scale_color_gradient(low = "cyan",high = "darkblue")+
  theme(axis.text.x = element_text(angle = 75,hjust = 1))
ggsave(fp(out,"suppfig6-table_pathway_enrichment_LGA_HSC_activated_dn_padj0.1.pdf"),height = 7)


#supp : sc [to update]
res_lin_act_sc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")
res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6&lineage_hmap=="HSC"]

ggplot(res_lin_act_sc[lineage_hmap%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],aes(x=avg_logFC,y=-log10(p_val_adj),col=p_val_adj<0.001&abs(avg_logFC)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(p_val_adj<0.001&
                                        abs(avg_logFC)>0.6&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("lineage_hmap")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("outputs/figures_epi_response/figure2/2D-sc_edger_by_lineage_lga_vs_ctrl_activated.pdf")


res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6,deg_sig:="deg"]
res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6&gene %in% hto_signature,deg_sig:="deg_hto"]
res_lin_act_sc[is.na(deg_sig),deg_sig:="other_gene"]


ggplot(res_lin_act_sc[lineage_hmap=="HSC"],aes(x=avg_logFC,y=-log10(p_val_adj),col=deg_sig))+
  geom_point()+
   geom_label_repel(aes(label=ifelse(p_val_adj<0.001&abs(avg_logFC)>0.6&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("red","blue","grey")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("outputs/figures_epi_response/figure2/2D-volcano_LGA_vs_Ctrl_sc_HSC_activation_response_signature_highlight.pdf")
res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6&lineage_hmap=="HSC"]



#figure 4 : Correl Meth - Expr  ####
source("scripts/utils/new_utils.R")
out<-fp(out0,"figure4")
dir.create(out)
#4A : hypermethylation  HSC activated DEGs 
hsc_act_res_pseudo<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]
hsc_act_res_pseudo[padj<0.05&abs(log2FoldChange)>0.5,deg_sig:="deg"]
hsc_act_res_pseudo[padj<0.05&abs(log2FoldChange)>0.5&gene %in% hto_signature,deg_sig:="deg_hto"]
hsc_act_res_pseudo[is.na(deg_sig),deg_sig:="non_deg"]
hsc_act_res_pseudo$deg_sig<-factor(hsc_act_res_pseudo$deg_sig,levels = c('non_deg',"deg","deg_hto"))

meth_res<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")

meth_res[,n.dmc:=sum(pval<0.001&abs(meth.change)>25),by="gene"]
meth_resg<-unique(meth_res[order(gene,pval)][,.(gene,gene_score_add,n.dmc)],by="gene")

res_int_pseudo<-merge(meth_resg,hsc_act_res_pseudo)
res_int_pseudo[padj<0.05&abs(log2FoldChange)>0.5] #293
res_int_pseudo[is.na(padj),padj:=1]
res_int_pseudo[gene=="ARID5A"]

table(res_int_pseudo$deg_sig)
# non_deg     deg deg_hto 
#    8923     146     147 

#DMCs enrich in DEGs
res_int_pseudo[,n.dmc.gene:=sum(n.dmc)/.N,by="deg_sig"]
ggplot(unique(res_int_pseudo[,.(deg_sig,n.dmc.gene)]))+geom_col(aes(x=deg_sig,y=n.dmc.gene,fill=deg_sig))
ggsave(fp(out,"4A-n.dmc_by_gene_barplot.pdf"))

ggplot(res_int_pseudo)+geom_boxplot(aes(x=deg_sig,y=n.dmc,fill=deg_sig),outlier.shape = NA)+coord_cartesian(ylim = c(0,3))
ggsave(fp(out,"4A-n.dmc_by_gene_boxplot.pdf"))

wilcox.test(res_int_pseudo[deg_sig=="non_deg"]$n.dmc,res_int_pseudo[deg_sig=="deg"]$n.dmc)
#^p=0.68
wilcox.test(res_int_pseudo[deg_sig=="non_deg"]$n.dmc,res_int_pseudo[deg_sig=="deg_hto"]$n.dmc)
#p-value = 0.0007632

#genescore enriched in DEGs


ggplot(res_int_pseudo,aes(x=padj<0.05&abs(log2FoldChange)>0.5,y=gene_score_add))+
  geom_boxplot(aes(fill=padj<0.05&abs(log2FoldChange)>0.5),size=0.5,alpha=0.6,outlier.shape = NA)+
  scale_fill_manual(values = c("grey","red"))+
  coord_cartesian(ylim = c(0,2000))+
  theme_classic()
ggsave(fp(out,"4A-gene_score_boxplot_degs_only.pdf"))

#in degs_sign
ggplot(res_int_pseudo,aes(x=deg_sig,y=gene_score_add))+
  geom_boxplot(aes(fill=deg_sig),size=0.5,alpha=0.6,outlier.shape = NA)+
  scale_fill_manual(values = c("grey","red","blue"))+
  coord_cartesian(ylim = c(0,2000))+
  theme_classic()
ggsave(fp(out,"4A-gene_score_boxplot.pdf"))

wilcox.test(res_int_pseudo[deg_sig=="non_deg"]$gene_score_add,res_int_pseudo[deg_sig=="deg"]$gene_score_add)
#p-value = 0.02839

wilcox.test(res_int_pseudo[deg_sig=="non_deg"]$gene_score_add,res_int_pseudo[deg_sig=="deg_hto"]$gene_score_add)
#p-value = 7.494e-07

wilcox.test(res_int_pseudo[deg_sig=="deg_hto"]$gene_score_add,res_int_pseudo[deg_sig=="deg"]$gene_score_add)
#p-value = 0.004262


#4B : volcano plot correl meth - gene expression

genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC","SESN2")

#degs only
res_int_pseudo[padj<0.05&abs(log2FoldChange)>0.5,deg_meth:="deg"]
#res_int_pseudo[padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>400,deg_meth:="deg_meth"]
res_int_pseudo[is.na(deg_meth),deg_meth:="non_degs"]

ggplot(res_int_pseudo,aes(x=log2FoldChange,y=gene_score_add,col=padj<0.05&abs(log2FoldChange)>0.5))+
  geom_point(aes(alpha=padj<0.05&abs(log2FoldChange)>0.5))+
    geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>400&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  coord_cartesian(ylim = c(0,3000))+
  theme_classic()
ggsave(fp(out,"4B-plot_gene_score_degs_lga_hsc_activated.pdf"))

plot(density(res_int_pseudo$gene_score_add))
abline(v=400)
#degs sign
ggplot(res_int_pseudo,aes(x=log2FoldChange,y=gene_score_add,col=deg_sig))+
  geom_point(aes(alpha=deg_sig))+
    geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>200&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red","blue"))+
  coord_cartesian(ylim = c(0,3000))+
  theme_classic()

ggsave(fp(out,"4B-plot_gene_score_degs_lga_hsc_activated.pdf"))
       
ggplot(res_int_pseudo,aes(x=log2FoldChange,y=log(gene_score_add),col=deg_sig))+
  geom_point(aes(alpha=deg_sig))+
    geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>200&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red","blue"))+
  coord_cartesian(ylim = c(0.00001,8))+
  theme_classic()

ggsave(fp(out,"4B-plot_gene_score_degs_lga_hsc_activated_log.pdf"))

ggplot(res_int_pseudo,aes(x=log2FoldChange,y=rank(gene_score_add),col=deg_sig))+
  geom_point(aes(alpha=deg_sig))+
    geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>200&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red","blue"))+
  theme_classic()


quants<-quantile(res_int_pseudo$gene_score_add,1:100/100)
res_int_pseudo[,gene_score_bin:=sum(gene_score_add>=quants),by="gene"]
ggplot(res_int_pseudo,aes(x=log2FoldChange,y=gene_score_bin,col=deg_sig))+
  geom_point()+
    geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>200&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red","blue"))+
  theme_classic()

ggsave(fp(out,"4B-plot_gene_score_bins_degs_lga_hsc_activated260821.pdf"))

ggplot(res_int_pseudo,aes(x=log2FoldChange,y=gene_score_bin,col=deg_sig))+
  geom_point(aes(alpha=deg_sig))+
    geom_label_repel(aes(label=ifelse(padj<0.05&abs(log2FoldChange)>0.5&gene_score_add>200&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red","blue"))+
  theme_classic()

ggsave(fp(out,"4B-plot_gene_score_bins_degs_lga_hsc_activated260821.pdf"))


#supp : sc HSC [to update]
meth_res<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
hsc_act_res<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")[lineage_hmap=="HSC"]
res_int<-merge(meth_res[,.(gene,gene_score_add)],hsc_act_res)
res_int[p_val_adj<0.001]


#figure 5 : SCENIC  ####
source("scripts/utils/new_utils.R")
out<-fp(out0,"figure5")
dir.create(out)

cbps<-readRDS("outputs/10-SCENIC/cbps_with_regulons_activity.rds")

#supp : heatmap regulon by lin
regul_lin <- sapply(split(colnames(cbps), cbps@meta.data$lineage_hmap),
                                     function(cells) {
                                       if(length(cells)>1){
                                         return(rowMeans(as.matrix(cbps@assays$TF_AUC@data)[,cells]))
                                       }else if (length(cells)==1){
                                         return(cbps@assays$TF_AUC@data[,cells])

                                           }else return(NA)

                                       })
regul_lin_scaled<-t(scale(t(regul_lin), center = T, scale=T))

regul_lin_dt<-melt(regul_lin_scaled)
colnames(regul_lin_dt)<-c("regulon","lineage","relative_activity")
regul_lin_dt<-data.table(regul_lin_dt)

regul_lin_dt[,top15:=relative_activity>=sort(relative_activity,decreasing = T)[15],by='lineage']

pdf(fp(out,"3A-heatmap_top15_tf_by_lineage.pdf"),width = 20,height = 5)
ComplexHeatmap::Heatmap(t(regul_lin_scaled[unique(regul_lin_dt[top15==T]$regulon),unique(cbps$lineage_hmap[!cbps$differentiated])]), name="Regulon activity")
dev.off()

#supp: Heatmap TF by celltype [to update]

regul_ct <- sapply(split(colnames(cbps), cbps@meta.data$cell_type_hmap),

                                     function(cells) {
                                       if(length(cells)>1){
                                         return(rowMeans(as.matrix(cbps@assays$SCENIC@data)[,cells]))
                                       }else if (length(cells)==1){
                                         return(cbps@assays$SCENIC@data[,cells])
                                           }else return(NA)

                                       })

regul_ct_scaled<-t(scale(t(regul_ct), center = T, scale=T))

regul_ct_dt<-melt(regul_ct_scaled)
colnames(regul_ct_dt)<-c("regulon","cell_type","relative_activity")
regul_ct_dt<-data.table(regul_ct_dt)

regul_ct_dt[,top15:=relative_activity>=sort(relative_activity,decreasing = T)[15],by='cell_type']

pdf(fp(out,"5A-heatmap_top15_tf_by_cell_type.pdf"),width = 8,height = 20)
ComplexHeatmap::Heatmap(regul_ct_scaled[unique(regul_ct_dt[top15==T]$regulon),], name="Regulon activity")
dev.off()

#regulons list
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")
paste(regulons_list[["EGR1"]],collapse = ", ")
regulons_df<-Reduce(rbind,lapply(names(regulons_list), function(tf)data.table(term=tf,gene=regulons_list[[tf]])))
fwrite(regulons_df,fp(out,"table-regulons_list.csv"))

#supp : regulons activation after stimulation
tf_diff_hto<-fread("outputs/10-SCENIC/regulon_activity_HTO_vs_not_by_lineage.csv.gz")

mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tfs_alt<-tf_diff_hto[p_val_adj<0.001&abs(avg_log2FC)>0.05&lineage=="HSC"&!str_detect(regulon,"e$")]$regulon
tfs_alt20<-head(tfs_alt,20)
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[tfs_alt20,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)

ggplot(tf_alt_act[lineage_hmap=="HSC"])+
  geom_boxplot(aes(x=regulon,y=activity,fill=hto))+theme_minimal()
ggsave(fp(out,"supp-boxplot_tf_activity_change_padj0_log2FC0.05_hto_vs_not.pdf"))

ggplot(tf_alt_act[lineage_hmap=="HSC"])+
  geom_boxplot(aes(x=hto,y=activity,fill=hto))+facet_wrap("regulon",nrow = 2)
  scale_fill_manual(values = c("white","brown1"))+theme_minimal()
ggsave(fp(out,"supp-boxplot_tf_activity_change_padj0_top20_log2FC0.05_hto_vs_not.pdf"))

tf_diff_hto[p_val_adj==0,p_val_adj:=10^-299]

ggplot(tf_diff_hto[lineage=="HSC"&!str_detect(regulon,"e$")],aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.001&avg_log2FC>0.05))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(p_val_adj<10^-250&avg_log2FC>0.05,regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  scale_y_continuous(limits = c(0,350))+
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(tf_diff_hto[lineage=="HSC"&!str_detect(regulon,"e$")&p_val_adj<0.001&avg_log2FC>0.05],aes(x=regulon,y=avg_log2FC,fill=-log10(p_val_adj)))+
  geom_col()+scale_x_discrete(limits=tf_diff_hto[lineage=="HSC"&!str_detect(regulon,"e$")&p_val_adj<0.001&avg_log2FC>0.05][order(-avg_log2FC)]$regulon)
  theme_minimal() 

avg_auc<-as.data.frame(AverageExpression(subset(cbps,lineage_hmap=="HSC"),assays = "TF_AUC",group.by = "hto")$TF_AUC)
avg_auc<-avg_auc[!str_detect(rownames(avg_auc),"e$"),]
avg_auc$steady_state<-avg_auc$`FALSE`
avg_auc$stimulated<-avg_auc$`TRUE`

ggplot(avg_auc)+
  geom_point(aes(x=steady_state,y=stimulated,col=rownames(avg_auc)%in%tfs_alt))+
  scale_color_manual(values = c("grey","red"))

#5A : TF activity alteration

res_tf_diff<-fread("outputs/10-SCENIC/regulon_activity_lga_vs_ctrl_HTO_by_lineage.csv.gz")

fwrite(res_tf_diff[p_val_adj<0.001&abs(avg_log2FC)>0.05&lineage=="HSC"&hto==T],fp(out,"supptabl6-res_tf_activity_alteration_lga_stim_hsc_padj0.001_avg_log2FC0.05.csv"))

res_tf_diff[,auc_change:=avg_log2FC]
table(res_tf_diff[p_val_adj<0.001&abs(avg_log2FC)>0.05]$lineage,res_tf_diff[p_val_adj<0.001&abs(avg_log2FC)>0.05]$hto)
  #            FALSE TRUE
  # B cell          0   10
  # Erythro-Mas     2    0
  # HSC             1    7
  # Lymphoid        0   22
  # MPP/LMPP        1    4
  # Myeloid         0    4
plot(density(res_tf_diff$avg_log2FC))
abline(v=-0.05)
res_tf_diff[p_val_adj<0.001&lineage=="HSC"&abs(avg_log2FC)>0.05]


ggplot(res_tf_diff[p_val_adj<0.001&lineage=="HSC"&abs(avg_log2FC)>0.05&hto==T])+
  geom_col(aes(x=regulon,y=auc_change,fill=-log10(p_val_adj)))+
 scale_x_discrete(limits = res_tf_diff[p_val_adj<0.001&lineage=="HSC"&abs(avg_log2FC)>0.05&hto==T][order(p_val_adj)]$regulon)

ggsave(fp(out,"5a-barplot_tf_change_hsc_lga_vs_ctrl.pdf"))

mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tfs_alt<-res_tf_diff[p_val_adj<0.001&lineage=="HSC"&abs(avg_log2FC)>0.05&hto==T]$regulon
tf_alt_act<-data.table(t(as.matrix(cbps@assays$TF_AUC@data[tfs_alt,])),keep.rownames = "cell")
tf_alt_act<-melt(tf_alt_act,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
tf_alt_act<-merge(tf_alt_act,mtd)
ggplot(tf_alt_act[lineage_hmap=="HSC"&hto==T&regulon!="KLF2e"])+
  geom_boxplot(aes(x=regulon,y=activity,fill=group))+theme_bw()
ggsave(fp(out,"5a-boxplot_tf_activity_change_padj0.001_logFC0.05_lga_hsc_activated.pdf"))

ggplot(tf_alt_act[lineage_hmap=="HSC"&regulon!="KLF2e"])+
  geom_boxplot(aes(x=hto,y=activity,fill=group))+facet_wrap("regulon",nrow = 1)
ggsave(fp(out,"5a-boxplot_tf_activity_change_padj0.001_logFC0.05_lga_hsc_activated_vs_not.pdf"))




for(tf in tfs_alt){
  print(tf)
  print(ggplot(tf_alt_act[regulon==tf])+geom_density(aes(x=activity)))
  tf_alt_act[regulon==tf,activ.thr:=as.numeric(readline("threshold: "))]
  
}

tf_alt_act[,activated:=activity>activ.thr]
tf_alt_act[,pct.activ:=sum(activated)/.N,,by=.(sample_hto,lineage_hmap,regulon)]

tf_alt_actl<-unique(tf_alt_act,by=c("sample_hto",'lineage_hmap',"regulon"))
ggplot(tf_alt_actl[lineage_hmap=="HSC"&regulon!="KLF2e"])+
  geom_boxplot(aes(x=hto,y=pct.activ,fill=group))+facet_wrap("regulon")+theme_minimal()
ggsave(fp(out,"5a-boxplot_pct_activ_tf_of_interest_by_sample_hsc.pdf"))

stat_hsc<-tf_alt_actl[lineage_hmap=="HSC"&regulon!="KLF2e"][,pvalue:=wilcox.test(pct.activ[group=="lga"],pct.activ[group=="ctrl"])$p.value,by=.(hto,regulon)]

unique(stat_hsc,by=c("hto","regulon"))[,.(regulon,hto,pvalue)]

#     regulon   hto     pvalue
#  1:  ARID5A FALSE 0.62820513
#  2:    EGR1 FALSE 0.44522145
#  3:    FOSB FALSE 1.00000000
#  4:    KLF4 FALSE 0.62820513
#  5:    KLF2 FALSE 0.83566434
#  6:     JUN FALSE 0.44522145
#  7:  ARID5A  TRUE 0.05927406
#  8:    EGR1  TRUE 0.18115218
#  9:    FOSB  TRUE 0.18115218
# 10:    KLF4  TRUE 0.10789211
# 11:    KLF2  TRUE 0.14185814
# 12:     JUN  TRUE 0.28238428

#5B : Regulons epigenetic programming
res_gsea<-fread("outputs/11-regulons_enrichment_genescore/res_gsea_genescore_regulons_hiconf.csv.gz")
res_gsea[,regulon:=pathway]
res_gsea[,pct.enrich:=size/size.regulon]

gs_enrich_regulons<-res_gsea[padj<0.005&NES>1.6]$regulon #38

res_gseaf<-res_gsea[regulon%in%gs_enrich_regulons]

regulons_ord<-res_gseaf[order(padj)]$pathway
ggplot(res_gseaf)+geom_point(aes(x=regulon,y=pct.enrich,col=-log10(padj),size=size.regulon))+
  scale_x_discrete(limits=regulons_ord)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
ggsave(fp(out,"5B-dotplot_gsea_regulons_genescore_padj0.005_NES1.6.pdf"))



#5C : Regulatory Network epigenetic programming
res_gsea_cl<-readRDS("outputs/11-regulons_enrichment_genescore/res_gsea_genescore_regulon_non_extended_clusterprofiler.rds")
res_gsea_cl_dt<-data.table(as.data.frame(res_gsea_cl))
res_gsea_cl_dt[p.adjust<0.05&NES>1.6]
res_gsea_cl_sub = res_gsea_cl[res_gsea_cl$p.adjust<0.05&res_gsea_cl$NES>1.6, asis=T]

pdf(fp(out,"5c-emapplot_gsea_genescore_regulon_clusterprofiler_padj0.005_NES1.6.pdf"),width=14,height=8)

emapplot(pairwise_termsim(res_gsea_cl_sub,showCategory = 33),showCategory = 33)

dev.off()
#5D : Regulatory network Expression alteration
stem_regulons<-c("ATF3","FOS","JUNB","KLF2","FOSB","JUN","EGR1","JUND","KLF10","KLF4","ARID5A")
res_gsea_cl_stem_net = res_gsea_cl_sub[res_gsea_cl_sub$ID%in%stem_regulons, asis=T]
class(res_gsea_cl_stem_net)

res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]
fold_changes<-res_degs$log2FoldChange
names(fold_changes)<-res_degs$gene

fold_changes_1<-fold_changes
fold_changes_1[fold_changes_1<(-1)]<- (-1)
fold_changes_1[fold_changes_1>(1)]<-1
pdf(fp(out,"5d-cnetplot_gsea_genescore_stemness_regulons_network_FoldChange_genes_highlight.pdf"),width=10,height=6)

cnetplot(res_gsea_cl_stem_net,showCategory = 11,
         foldChange = fold_changes_1,cex_category=0.8,cex_label_category=0.8,
         node_label = "category",
         )+scale_color_gradient2(limits=c(-max(abs(fold_changes_1)),max(abs(fold_changes_1))),
                                 high = scales::muted("red"),low =scales::muted("blue"))
dev.off()



fold_changes<-res_degs[padj<0.05]$log2FoldChange
names(fold_changes)<-res_degs[padj<0.05]$gene

fold_changes_bin<-fold_changes
fold_changes_bin[fold_changes_bin<0]<- (-1)
fold_changes_bin[fold_changes_bin>0]<-1

pdf(fp(out,"5d-cnetplot_gsea_genescore_stemness_regulons_network_DEGs_padj0.05_highlight.pdf"),width=10,height=6)

cnetplot(res_gsea_cl_stem_net,showCategory = 11,
         foldChange = fold_changes_bin,cex_category=0.8,cex_label_category=0.8,
         node_label = "category",
         )
dev.off()

#test enrichment of downregulated in this network
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")
genes_stem_reg<-Reduce(union,regulons_list[stem_regulons])

res_degs[gene%in%genes_stem_reg]#456
res_degs[padj<0.05&log2FoldChange<0&gene%in%genes_stem_reg] #53
53/456 #12%
over_repr_test_simple(set1 = res_degs[padj<0.05&log2FoldChange<0]$gene,
                      set2=genes_stem_reg,size_universe = nrow(res_degs))

#p= 1.956816e-19

fold_changes<-res_degs[padj<0.05&abs(log2FoldChange)>0.5]$log2FoldChange
names(fold_changes)<-res_degs[padj<0.05&abs(log2FoldChange)>0.5]$gene

fold_changes_bin<-fold_changes
fold_changes_bin[fold_changes_bin<0]<- (-1)
fold_changes_bin[fold_changes_bin>0]<-1

pdf(fp(out,"5d-cnetplot_gsea_genescore_stemness_regulons_network_DEGs_padj0.05_lFC0.5_highlight.pdf"),width=10,height=6)

cnetplot(res_gsea_cl_stem_net,showCategory = 11,
         foldChange = fold_changes_bin,cex_category=0.8,cex_label_category=0.8,
         node_label = "category",
         )
dev.off()

#with label
pdf(fp(out,"5d-cnetplot_gsea_genescore_stemness_regulons_network_FoldChange_label_genes.pdf"),width=14,height=8)
cnetplot(res_gsea_cl_stem_net,showCategory = 11,
         foldChange = fold_changes_bin,cex_category=0.8,cex_label_category=0.8,cex_label_gene=0.5,
         node_label = "gene",
         )+scale_color_gradient2(limits=c(-max(abs(fold_changes_1)),max(abs(fold_changes_1))),
                                 high = scales::muted("red"),low =scales::muted("blue"))
dev.off()


#enrichment regulons - degs
res_gsea_degs<-fread("outputs/11-regulons_enrichment_genescore/res_gsea_regulons_degs.csv")
res_gsea_degs[padj<0.05] #46
res_gsea_degs[padj<0.05&NES<0] #16
#     pathway         pval         padj   log2err         ES       NES size
#  1:  ARID5A 2.166824e-06 1.430104e-05 0.6272567 -0.6541719 -2.434705   24
#  2:   BRCA1 3.564308e-05 1.857192e-04 0.5573322 -0.4309827 -2.055753   64
#  3:    CTCF 1.426337e-03 4.381913e-03 0.4550599 -0.6347891 -1.961397   13
#  4:    E2F1 4.465968e-07 3.158078e-06 0.6749629 -0.3520359 -2.007411  182
#  5:    EGR1 3.416092e-07 2.601486e-06 0.6749629 -0.6717122 -2.517119   26
#  6:     FOS 1.789072e-09 2.530259e-08 0.7881868 -0.4216395 -2.342546  140
#  7:    FOSB 9.421371e-11 1.865431e-09 0.8390889 -0.6252559 -2.798122   50
#  8:     JUN 2.120224e-14 1.049511e-12 0.9759947 -0.5340468 -2.866158  115
#  9:    JUNB 5.370605e-06 3.127587e-05 0.6105269 -0.3345703 -1.887180  164
# 10:   KLF10 6.433794e-09 5.307880e-08 0.7614608 -0.4772567 -2.413191   89
# 11:    KLF2 3.298653e-11 8.164167e-10 0.8513391 -0.5040806 -2.602989   98
# 12:    KLF4 3.695465e-09 3.658511e-08 0.7749390 -0.6298719 -2.707370   42
# 13:    KLF6 5.285320e-03 1.341658e-02 0.4070179 -0.3702973 -1.644896   48
# 14:   NFKB2 2.918157e-03 8.024931e-03 0.4317077 -0.2818358 -1.564189  138
# 15:    RELB 3.386403e-03 8.822472e-03 0.4317077 -0.3801516 -1.688670   48
# 16:    YBX1 8.743607e-04 3.091490e-03 0.4772708 -0.3553222 -1.764431   79

res_gsea_degs[padj<0.05&NES>0] #30

fwrite(res_gsea_degs[padj<0.05],fp(out,"res_gsea_regulons_degs_padj0.05.csv"))

#merge res_gsea_degs, res_gsea_gs

res_gsea_gs<-fread("outputs/11-regulons_enrichment_genescore/res_gsea_genescore_regulons_hiconf.csv.gz",
                   select = c(1:3,6,9),
                     col.names = c("regulon","pval.gs","padj.gs","NES.gs","size.regulon"))
res_gsea_degs<-fread("outputs/11-regulons_enrichment_genescore/res_gsea_regulons_degs.csv",
                     select = c(1:3,6,9),
                     col.names = c("regulon","pval.degs","padj.degs","NES.degs","size.regulon"))

res_merge<-merge(res_gsea_gs,res_gsea_degs,all.x=T)

res_merge
ggplot(res_merge,aes(x=-log10(padj.gs),y=-log10(padj.degs),col=padj.gs<0.01&padj.degs<0.01))+
  geom_point(aes(size=size.regulon))+
    geom_label_repel(aes(label=ifelse(padj.degs<0.01&padj.gs<0.01&regul,regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()

gs_enrich_regulons<-res_merge[padj.gs<0.01&NES.gs>1.6]$regulon #33

ggplot(res_merge,aes(y=NES.gs,x=NES.degs,col=padj.gs<0.01&padj.degs<0.01))+
  geom_point(aes(size=size.regulon))+
    geom_label_repel(aes(label=ifelse(padj.gs<0.01&padj.degs<0.01,regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = "bottom")


ggplot(res_merge,aes(y=NES.gs,x=NES.degs,col=padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6))+
  geom_point(aes(size=size.regulon))+
    geom_label_repel(aes(label=ifelse(padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6,regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = "bottom")

ggsave(fp(out,"5B-dotplot_NES_gsea_regulons_genescore_degs.pdf"))



ggplot(res_merge,aes(x=NES.gs,y=NES.degs,col=padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6))+
  geom_point(aes(size=size.regulon))+
    geom_label_repel(aes(label=ifelse(padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6,regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = "bottom")

ggsave(fp(out,"5B-dotplot_NES_gsea_regulons_genescore_degs2.pdf"))

ggplot(res_merge,aes(y=NES.gs,x=NES.degs,col=padj.gs<0.01&NES.gs>1.6&padj.degs<0.01&abs(NES.degs)>1.6))+
  geom_point(aes(size=size.regulon))+
    geom_label_repel(aes(label=ifelse(regulon%in% c("ARID5A",'EGR1','FOSB',"KLF4","KLF2","JUN"),regulon,"")),
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red"))+
  theme_minimal()+theme(legend.position = "bottom")



res_mergef<-res_merge[regulon%in%gs_enrich_regulons]

regulons_ord<-res_mergef[order(NES.degs)]$regulon
ggplot(res_mergef)+geom_point(aes(x=regulon,y=-log10(padj.degs),col=NES.degs,size=size.regulon))+
  geom_hline(yintercept = 2)+
  scale_x_discrete(limits=regulons_ord)+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_color_gradient2(,low = "darkblue",high ='darkred' )
ggsave(fp(out,"5B-dotplot_gsea_regulons_genescore_padj0.01_NES1.6_degs_merged.pdf"))




#figure 7 : lineage bias  ####
source("scripts/utils/new_utils.R")
out<-fp(out0,"figure7")
dir.create(out)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

#7A : lineage bias
mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd$lineage_hmap<-factor(mtd$lineage_hmap,
                         levels = c("LT-HSC",
                                    "HSC",
                                    "MPP/LMPP",
                                    "Lymphoid",
                                    "B cell",
                                    "T cell",
                                    "Erythro-Mas",
                                    "Mk/Er",
                                    "Myeloid",
                                    "DC"))
#at group level
ggplot(mtd)+geom_bar(aes(x=group,fill=lineage_hmap),position ="fill" )+facet_wrap("hto")
ggsave(fp(out,"7a-lineage_distrib_by_group.pdf"))
chisq.test(table(mtdf[hto==T]$lineage_hmap,mtdf[hto==T]$group)) 
# Pearson's Chi-squared test
# 
# data:  table(mtdf[hto == T]$lineage_hmap, mtdf[hto == T]$group)
# X-squared = 429.88, df = 9, p-value < 2.2e-16

chisq.test(table(mtdf[hto==F]$lineage_hmap,mtdf[hto==F]$group)) 
# 	Pearson's Chi-squared test
# 
# data:  table(mtd[hto == F]$group, mtd[hto == F]$lineage_hmap)
# X-squared = 340.98, df = 9, p-value < 2.2e-16

#at sample lvl
mtd[,n.cells:=.N,"sample_hto"]
mtd[,n.cells.lin:=.N,c("sample_hto","lineage_hmap")]
mtd[,pct.lin:=.N/n.cells,c("sample_hto","lineage_hmap")]
mtsl<-unique(mtd[lineage_hmap%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],by=c("sample_hto","lineage_hmap"))
mtsl$lineage_hmap<-factor(mtsl$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid"))

#rm outlyers
ggplot(mtsl)+geom_boxplot(aes(x=hto,y=pct.lin,fill=group))+facet_wrap("lineage_hmap")
  # scale_x_discrete(guide = guide_axis(angle = 45))+
  # scale_y_continuous(expand = c(0,0))

mtsl[,is.outlier:=pct.lin%in%boxplot.stats(pct.lin)$out,by=.(group,hto,lineage_hmap)]
mtsl[(is.outlier)]$sample_hto 
mtsl[,is.outlier:=pct.lin%in%boxplot.stats(pct.lin)$out,by=.(group,hto,lineage_hmap)]
mtsl[,outlier:=sample_hto%in%sample_hto[(is.outlier)]]
unique(mtsl[(outlier)]$sample_hto) #"lgaF552FALSE"  "ctrlF523TRUE"  "ctrlM537FALSE" "ctrlM530TRUE"
mtslf<-mtsl[(!outlier)]
ggplot(mtslf)+
  geom_boxplot(aes(x=lineage_hmap,y=pct.lin,fill=group))+facet_wrap("hto")
ggsave(fp(out,"7a-boxplot_lineage_proport_in_samples_by_group.pdf"))


mtslf[,pval:=wilcox.test(pct.lin[group=="lga"],pct.lin[group=="ctrl"])$p.value,by=c("hto","lineage_hmap")]

mthl<-unique(mtslf[,.(hto,lineage_hmap,pval)])
mthl[,padj:=p.adjust(pval),by="hto"]
mthl
#       hto lineage_hmap       pval       padj
#  1: FALSE     MPP/LMPP 0.66233766 1.00000000
#  2: FALSE          HSC 1.00000000 1.00000000
#  3: FALSE  Erythro-Mas 0.42857143 1.00000000
#  4: FALSE      Myeloid 0.79220779 1.00000000
#  5: FALSE     Lymphoid 0.79220779 1.00000000
#  6: FALSE       LT-HSC 0.28571429 1.00000000
#  7:  TRUE     MPP/LMPP 0.13203463 0.66017316
#  8:  TRUE          HSC 0.01515152 0.09090909
#  9:  TRUE  Erythro-Mas 1.00000000 1.00000000
# 10:  TRUE     Lymphoid 0.69913420 1.00000000
# 11:  TRUE      Myeloid 0.24025974 0.96103896
# 12:  TRUE       LT-HSC 0.73015873 1.00000000



#change of proportion HSC ?
mts<-data.table(sample=unique(mtsl$sample_hto),
                n.hscs=as.vector(table(mtd[lineage_hmap=="HSC"]$sample_hto)[unique(mtsl$sample_hto)]),
                n.cells=as.vector(table(mtd$sample_hto)[unique(mtsl$sample_hto)]),
                group=unique(mtsl,by="sample_hto")$group,
                hto=unique(mtsl,by="sample_hto")$hto,
                batch=unique(mtsl,by="sample_hto")$batch,
                sex=unique(mtsl,by="sample_hto")$sex)

glm.hsc <- stats::glm(n.hscs~n.cells+group+batch+sex,family=poisson(),data = mts[(hto)])

glm.hsc.null <- stats::glm(n.hscs~group, family="poisson", data = mts[(hto)])

glm.hsc.nocell <- stats::glm(n.hscs~group+batch+sex, family="poisson", data = mts[(hto)])

performance::compare_performance(glm.hsc.null, glm.hsc, glm.hsc.nocell)

performance::check_model(glm.hsc)
performance::check_model(stats::glm(n.hscs~n.cells+group+batch+sex,family=gaussian(),data = mts[(hto)]))
anova(glm.hsc.null, glm.hsc)

summary(glm.hsc)
# Call:
# stats::glm(formula = n.hscs ~ n.cells + group + batch + sex, 
#     family = poisson(), data = mts[hto == T])
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -18.6887   -5.4115   -0.3179    3.2327   12.7806  
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  5.342e+00  5.345e-02  99.931  < 2e-16 ***
# n.cells      5.150e-04  4.654e-05  11.067  < 2e-16 ***
# grouplga    -1.528e-01  4.907e-02  -3.113  0.00185 ** 
# batchcbp3   -1.593e+00  1.539e-01 -10.355  < 2e-16 ***
# batchcbp4   -2.317e-02  5.005e-02  -0.463  0.64342    
# batchcbp8   -2.957e-01  5.627e-02  -5.255 1.48e-07 ***
# sexM        -1.114e-02  5.545e-02  -0.201  0.84072    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
#     Null deviance: 1704.38  on 13  degrees of freedom
# Residual deviance:  934.65  on  7  degrees of freedom
# AIC: 1049.3
# 
# Number of Fisher Scoring iterations: 5



#change of proportion MPP ?
mts<-data.table(sample=unique(mtsl$sample_hto),
                n.mpps=as.vector(table(mtd[lineage_hmap=="MPP/LMPP"]$sample_hto)[unique(mtsl$sample_hto)]),
                n.cells=as.vector(table(mtd$sample_hto)[unique(mtsl$sample_hto)]),
                group=unique(mtsl,by="sample_hto")$group,
                hto=unique(mtsl,by="sample_hto")$hto,
                batch=unique(mtsl,by="sample_hto")$batch,
                sex=unique(mtsl,by="sample_hto")$sex)

glm.mpp<-stats::glm(n.mpps~n.cells+group+batch+sex,family=poisson(),data = mts[hto==T])

summary(glm.mpp)
# Call:
# stats::glm(formula = n.mpps ~ n.cells + group + batch + sex, 
#     family = poisson(), data = mts[hto == T])
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -9.3855  -1.0268   0.1843   2.0875   6.2796  
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.778e+00  5.011e-02  95.356  < 2e-16 ***
# n.cells      1.410e-03  4.196e-05  33.609  < 2e-16 ***
# grouplga    -4.083e-01  4.188e-02  -9.749  < 2e-16 ***
# batchcbp3   -1.165e+00  1.501e-01  -7.766 8.13e-15 ***
# batchcbp4    3.926e-01  4.649e-02   8.444  < 2e-16 ***
# batchcbp8    8.612e-01  4.277e-02  20.137  < 2e-16 ***
# sexM        -6.275e-01  5.099e-02 -12.305  < 2e-16 ***

mts<-data.table(sample=unique(mtsl$sample_hto),
                n.cells_lin=as.vector(table(mtd[lineage_hmap=="Lymphoid"]$sample_hto)[unique(mtsl$sample_hto)]),
                n.cells=as.vector(table(mtd$sample_hto)[unique(mtsl$sample_hto)]),
                group=unique(mtsl,by="sample_hto")$group,
                hto=unique(mtsl,by="sample_hto")$hto,
                batch=unique(mtsl,by="sample_hto")$batch,
                sex=unique(mtsl,by="sample_hto")$sex)

glm.lin<-stats::glm(n.cells_lin~n.cells+group+batch+sex,family=poisson(),data = mts[hto==T])

summary(glm.lin)
# Call:
# stats::glm(formula = n.cells_lin ~ n.cells + group + batch + 
#     sex, family = poisson(), data = mts[hto == T])
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -4.5792  -1.6967  -0.1383   0.3321   5.7411  
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  3.3013612  0.1229754  26.846  < 2e-16 ***
# n.cells      0.0014065  0.0001519   9.261  < 2e-16 ***
# grouplga    -0.6181826  0.1599613  -3.865 0.000111 ***
# batchcbp3   -0.7466297  0.2649369  -2.818 0.004830 ** 
# batchcbp4    0.0499724  0.1358303   0.368 0.712945    
# batchcbp8   -0.9718284  0.1565491  -6.208 5.37e-10 ***
# sexM        -0.5136635  0.1622625  -3.166 0.001547 ** 

#change of cells proportions ?

glm.lin_hto<-stats::glm(n.cells.lin~n.cells+lineage_hmap*group+batch+sex,family=poisson(),data = mtsl[hto==T])
 summary(glm.lin_hto)

# Call:
# stats::glm(formula = n.cells.lin ~ n.cells + lineage_hmap * group + 
#     batch + sex, family = poisson(), data = mtsl[hto == T])
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -25.1204   -2.6923   -0.5314    1.3302   15.3360  
# 
# Coefficients:
#                                    Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                       1.621e-01  2.518e-01   0.644   0.5197    
# n.cells                           1.090e-03  2.884e-05  37.781  < 2e-16 ***
# lineage_hmapHSC                   4.739e+00  2.510e-01  18.881  < 2e-16 ***
# lineage_hmapMPP/LMPP              4.968e+00  2.508e-01  19.811  < 2e-16 ***
# lineage_hmapErythro-Mas           2.933e+00  2.558e-01  11.465  < 2e-16 ***
# lineage_hmapMyeloid               2.242e+00  2.614e-01   8.577  < 2e-16 ***
# lineage_hmapLymphoid              2.992e+00  2.555e-01  11.714  < 2e-16 ***
# grouplga                          4.609e-03  3.293e-01   0.014   0.9888    
# batchcbp3                        -8.902e-01  7.955e-02 -11.190  < 2e-16 ***
# batchcbp4                         1.977e-01  3.086e-02   6.408 1.47e-10 ***
# batchcbp8                         2.743e-01  3.090e-02   8.877  < 2e-16 ***
# sexM                             -3.805e-01  3.439e-02 -11.062  < 2e-16 ***
# lineage_hmapHSC:grouplga         -5.710e-01  3.302e-01  -1.730   0.0837 .  
# lineage_hmapMPP/LMPP:grouplga    -2.014e-01  3.297e-01  -0.611   0.5412    
# lineage_hmapErythro-Mas:grouplga -1.261e-01  3.361e-01  -0.375   0.7076    
# lineage_hmapMyeloid:grouplga      1.136e-01  3.422e-01   0.332   0.7398    
# lineage_hmapLymphoid:grouplga    -2.277e-01  3.360e-01  -0.678   0.4980    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
#     Null deviance: 19907  on 80  degrees of freedom
# Residual deviance:  2302  on 64  degrees of freedom
# AIC: 2799
# 
# Number of Fisher Scoring iterations: 5
 
 
glm.lin<-stats::glm(n.cells.lin~n.cells+lineage_hmap*group*hto+batch+sex,family=poisson(),data = mtsl)
 summary(glm.lin)
 
#  Call:
# stats::glm(formula = n.cells.lin ~ n.cells + lineage_hmap * group * 
#     hto + batch + sex, family = poisson(), data = mtsl)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -27.061   -4.371   -1.138    2.654   24.554  
# 
# Coefficients: (1 not defined because of singularities)
#                                            Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                              -3.700e-01  1.270e-01  -2.913  0.00358 ** 
# n.cells                                   6.254e-04  1.244e-05  50.256  < 2e-16 ***
# lineage_hmapHSC                           3.415e+00  9.910e-02  34.456  < 2e-16 ***
# lineage_hmapMPP/LMPP                      4.497e+00  9.811e-02  45.836  < 2e-16 ***
# lineage_hmapErythro-Mas                   2.974e+00  9.993e-02  29.767  < 2e-16 ***
# lineage_hmapMyeloid                       2.081e+00  1.032e-01  20.165  < 2e-16 ***
# lineage_hmapLymphoid                      2.508e+00  1.013e-01  24.760  < 2e-16 ***
# grouplga                                 -1.892e-01  1.444e-01  -1.310  0.19030    
# htoTRUE                                   6.841e-01  2.789e-01   2.453  0.01416 *  
# batchcbp0_lga                             1.631e+00  6.267e-02  26.031  < 2e-16 ***
# batchcbp2                                 9.624e-02  2.356e-02   4.085 4.42e-05 ***
# batchcbp3                                -9.467e-01  7.774e-02 -12.177  < 2e-16 ***
# batchcbp4                                -2.852e-02  2.508e-02  -1.137  0.25546    
# batchcbp6a                                1.440e+00  5.546e-02  25.967  < 2e-16 ***
# batchcbp6b                                1.631e+00  7.010e-02  23.265  < 2e-16 ***
# batchcbp6c                                1.621e+00  6.341e-02  25.556  < 2e-16 ***
# batchcbp7a                                1.610e+00  6.154e-02  26.167  < 2e-16 ***
# batchcbp7b                                1.474e+00  5.532e-02  26.640  < 2e-16 ***
# batchcbp7c                                1.298e+00  5.379e-02  24.139  < 2e-16 ***
# batchcbp8                                        NA         NA      NA       NA    
# sexM                                      7.526e-02  1.217e-02   6.182 6.32e-10 ***
# lineage_hmapHSC:grouplga                  2.482e-01  1.455e-01   1.706  0.08806 .  
# lineage_hmapMPP/LMPP:grouplga            -3.574e-02  1.444e-01  -0.248  0.80449    
# lineage_hmapErythro-Mas:grouplga         -3.076e-01  1.475e-01  -2.085  0.03706 *  
# lineage_hmapMyeloid:grouplga             -5.657e-02  1.518e-01  -0.373  0.70932    
# lineage_hmapLymphoid:grouplga            -4.499e-01  1.503e-01  -2.993  0.00276 ** 
# lineage_hmapHSC:htoTRUE                   1.320e+00  2.698e-01   4.893 9.92e-07 ***
# lineage_hmapMPP/LMPP:htoTRUE              4.674e-01  2.693e-01   1.736  0.08259 .  
# lineage_hmapErythro-Mas:htoTRUE          -4.524e-02  2.746e-01  -0.165  0.86915    
# lineage_hmapMyeloid:htoTRUE               1.581e-01  2.811e-01   0.562  0.57378    
# lineage_hmapLymphoid:htoTRUE              4.811e-01  2.748e-01   1.751  0.08000 .  
# grouplga:htoTRUE                          5.211e-01  3.589e-01   1.452  0.14645    
# lineage_hmapHSC:grouplga:htoTRUE         -7.784e-01  3.608e-01  -2.157  0.03097 *  
# lineage_hmapMPP/LMPP:grouplga:htoTRUE    -1.248e-01  3.599e-01  -0.347  0.72868    
# lineage_hmapErythro-Mas:grouplga:htoTRUE  2.224e-01  3.671e-01   0.606  0.54461    
# lineage_hmapMyeloid:grouplga:htoTRUE      2.110e-01  3.743e-01   0.564  0.57284    
# lineage_hmapLymphoid:grouplga:htoTRUE     2.631e-01  3.681e-01   0.715  0.47477    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
#     Null deviance: 80741.9  on 156  degrees of freedom
# Residual deviance:  8662.3  on 121  degrees of freedom
# AIC: 9724.1
# 
# Number of Fisher Scoring iterations: 5

#7B : CFU bias

#7C : colony expansion alteration


#table
#number of cells/ech

# genescore, gene

#stimulation signature, pval, FC

#DEGs LGA HSC activ



 
 #figure 6 : Pseudotime ####
source("scripts/utils/new_utils.R")
library(Seurat)
out<-fp(out0,"figure6")
dir.create(out)
#6A : Pseudotime 
pseudo_mtd<-fread("outputs/12A-Pseudotime_integrated/metadata_pseudotime_ComputRoot.csv")
cbps<-AddMetaData(cbps,
                  metadata = data.frame(pseudo_mtd[,.(bc,pseudotime)],row.names = "bc"),
                  col.name = "pseudotime")

Idents(cbps)<-"lineage_hmap"
FeaturePlot(cbps,"pseudotime",label=T,reduction = "ref.umap",cols = c("aquamarine1","darkorange2"))
ggsave(fp(out,"umap_pseudotime.pdf"))

mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd$lineage_hmap<-factor(mtd$lineage,levels = c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid"))

ggplot(mtd)+geom_boxplot(aes(y=pseudotime,x=lineage_hmap),alpha=0.7)+theme_minimal()
ggsave(fp(out,"boxplot_pseudotime_by_lineage.pdf"))

#cor
summary(lm(pseudotime~lineage_hmap,data=mtd))

# Call:
# lm(formula = pseudotime ~ lineage_hmap, data = mtd)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -23.933  -2.542  -0.037   2.215  27.842 
# 
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               1.4181     0.2576   5.505 3.71e-08 ***
# lineage_hmapHSC           6.9209     0.2602  26.603  < 2e-16 ***
# lineage_hmapMPP/LMPP     11.9950     0.2588  46.350  < 2e-16 ***
# lineage_hmapErythro-Mas  15.1831     0.2640  57.520  < 2e-16 ***
# lineage_hmapMyeloid      18.0431     0.2708  66.619  < 2e-16 ***
# lineage_hmapLymphoid     23.0087     0.2671  86.154  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.856 on 45499 degrees of freedom
#   (2492 observations deleted due to missingness)
# Multiple R-squared:  0.5451,	Adjusted R-squared:  0.545 
# F-statistic: 1.09e+04 on 5 and 45499 DF,  p-value: < 2.2e-16

#6B Pseudotime alteration

#distribution 
#rm outliers
mtdf<-mtd[sample_hto%in%unique(mtslf$sample_hto)]

ggplot(mtdf)+geom_boxplot(aes(x=pseudotime,y=group,fill=group),alpha=0.7)+facet_wrap('hto',nrow = 2)+theme_minimal()

wilcox.test(mtdf[group=="ctrl"&hto==T]$pseudotime,mtdf[group=="lga"&hto==T]$pseudotime) #pvalue < 2.2e-16

wilcox.test(mtdf[group=="ctrl"&hto==F]$pseudotime,mtdf[group=="lga"&hto==F]$pseudotime) #pvalue < 2.2e-16

ggsave(fp(out,"boxplot_pseudotime_by_group_hto.pdf"))

# pseudotime influencé par group:hto ?
ggplot(mtdf)+geom_density(aes(x=pseudotime,fill=group),alpha=0.5,alpha=0.7)+facet_wrap('hto')+theme_minimal()
ggsave(fp(out,"density_pseudotime_by_group_hto.pdf"))

#choose model ~ data
#mtd
ggplot(mtdf)+
  geom_density(aes(x=pseudotime)) # binomial
ggplot(mtdf[(hto)])+
  geom_density(aes(x=pseudotime,col = group,group=sample)) #

ggplot(mtd[(hto)])+
  geom_density(aes(x=pseudotime,col = group))


mtdf[,n_cells.sample:=.N,by="sample_hto"]
mtdf[,n_cells.scaled:=scale(n_cells.sample)]
library(lmerTest)

lm_null<-lmer(formula = pseudotime~ 1 + group + (1|sample_hto), data = mtdf[!is.na(pseudotime)&(hto)])
performance::check_model(lm_null)
lm_lin<-lmer(formula = pseudotime ~ 1 +group+ group:lineage_hmap + (1|sample_hto), data = mtdf[!is.na(pseudotime)&(hto)])

lm_ls<-lmer(formula = pseudotime ~ 1 +group+ group:lineage_hmap+n_cells.sample +lineage_hmap+ (1|sample_hto), data = mtdf[!is.na(pseudotime)&(hto)])

lm_lsb<-lmer(formula = pseudotime ~ 1 +group+ group:lineage_hmap+n_cells.sample +batch+ (1|sample_hto), data = mtdf[!is.na(pseudotime)&(hto)])
lm_lsbs<-lmer(formula = pseudotime ~ 1 +group+ group:lineage_hmap+lineage_hmap+n_cells.sample +batch+ sex+(1|sample_hto), data = mtdf[!is.na(pseudotime)&(hto)])

performance::compare_performance(lm_null, lm_lin,lm_ls,lm_lsb,lm_lsbs) #lm_lin best

summary(lm_lin)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pseudotime ~ 1 + group + group:lineage_hmap + (1 | sample_hto)
#    Data: mtdf[!is.na(pseudotime) & (hto)]
# 
# REML criterion at convergence: 59915
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -5.1702 -0.6330 -0.0157  0.5318  6.1511 
# 
# Random effects:
#  Groups     Name        Variance Std.Dev.
#  sample_hto (Intercept)  0.4478  0.6692  
#  Residual               14.3201  3.7842  
# Number of obs: 10888, groups:  sample_hto, 12
# 
# Fixed effects:
#                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                           6.6663     1.1268  2212.3544   5.916 3.81e-09 ***
# grouplga                             -5.5674     1.4248  1517.7642  -3.907 9.74e-05 ***
# groupctrl:lineage_hmapHSC             2.0080     1.0971 10868.5204   1.830   0.0672 .  
# grouplga:lineage_hmapHSC              8.0747     0.8325 10869.7961   9.700  < 2e-16 ***
# groupctrl:lineage_hmapMPP/LMPP        6.5944     1.0971 10867.4584   6.011 1.91e-09 ***
# grouplga:lineage_hmapMPP/LMPP        12.8196     0.8294 10868.4432  15.456  < 2e-16 ***
# groupctrl:lineage_hmapErythro-Mas    10.1700     1.1176 10867.9724   9.100  < 2e-16 ***
# grouplga:lineage_hmapErythro-Mas     15.1684     0.8444 10869.4526  17.963  < 2e-16 ***
# groupctrl:lineage_hmapMyeloid        12.6886     1.1546 10866.9057  10.989  < 2e-16 ***
# grouplga:lineage_hmapMyeloid         18.0219     0.8541 10868.3455  21.102  < 2e-16 ***
# groupctrl:lineage_hmapLymphoid       17.4038     1.1150 10868.6054  15.609  < 2e-16 ***
# grouplga:lineage_hmapLymphoid        24.0667     0.8452 10869.2256  28.474  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lm_b_lin<-lmer(formula = pseudotime ~ 1 +group+ group:lineage_hmap + (1|sample_hto), data = mtdf[!is.na(pseudotime)&!(hto)])
summary(lm_b_lin)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: pseudotime ~ 1 + group + group:lineage_hmap + (1 | sample_hto)
#    Data: mtdf[!is.na(pseudotime) & !(hto)]
# 
# REML criterion at convergence: 154948.2
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -6.3029 -0.6327 -0.0141  0.5659  5.6513 
# 
# Random effects:
#  Groups     Name        Variance Std.Dev.
#  sample_hto (Intercept)  0.07964 0.2822  
#  Residual               14.57641 3.8179  
# Number of obs: 28077, groups:  sample_hto, 11
# 
# Fixed effects:
#                                     Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)                           1.0601     0.4000   992.0184   2.651  0.00816 ** 
# grouplga                             -0.3533     0.7188  2061.0833  -0.492  0.62310    
# groupctrl:lineage_hmapHSC             6.8697     0.3883 28062.1933  17.691  < 2e-16 ***
# grouplga:lineage_hmapHSC              7.1566     0.5876 28063.4479  12.180  < 2e-16 ***
# groupctrl:lineage_hmapMPP/LMPP       12.3462     0.3847 28064.8946  32.092  < 2e-16 ***
# grouplga:lineage_hmapMPP/LMPP        12.5135     0.5852 28064.9969  21.385  < 2e-16 ***
# groupctrl:lineage_hmapErythro-Mas    15.3401     0.3947 28064.2646  38.863  < 2e-16 ***
# grouplga:lineage_hmapErythro-Mas     16.1351     0.5954 28064.9787  27.101  < 2e-16 ***
# groupctrl:lineage_hmapMyeloid        18.6161     0.4068 28064.9726  45.765  < 2e-16 ***
# grouplga:lineage_hmapMyeloid         18.4699     0.6033 28059.7529  30.615  < 2e-16 ***
# groupctrl:lineage_hmapLymphoid       20.8289     0.4078 28064.4742  51.075  < 2e-16 ***
# grouplga:lineage_hmapLymphoid        24.0670     0.6050 28064.6660  39.779  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#figures : percentage of cells below each pseudotime
mtdf
n_bins<-30
pseudotime_thr<-1:n_bins*(max(mtdf$pseudotime,na.rm = T)/n_bins)
pseudo_bins<-data.table(bin=0:(n_bins-1),
                        pseudotime_thr=pseudotime_thr)
mtdf[,bin:=sum(pseudotime>=pseudotime_thr),by="bc"]
summary(mtdf$bin)

pseudo_bins_mtd<-merge(pseudo_bins,mtdf)

pseudo_bins_mtd[,n.cells:=.N,by="sample_hto"]
pseudo_bins_mtd[,n.cells.bin:=.N,by=.(sample_hto,bin)]
pseudo_bins_<-unique(pseudo_bins_mtd,by=c("sample_hto","bin"))
for(gr in c("ctrl","lga")){
  for(h in c(FALSE,TRUE)){
    samples<-unique(pseudo_bins_[group==gr&hto==h]$sample)
    for(i in 1:max(pseudo_bins_$bin)){
      for(s in setdiff(samples,pseudo_bins_[bin==i&group==gr&hto==h]$sample)){
        pseudo_bins_<-rbind(pseudo_bins_,pseudo_bins_[bin==i-1&group==gr&hto==h&sample==s][,bin:=bin+1][,n.cells.bin:=0])
        
      }
      
    }
    
  }
  
}

pseudo_bins_[,pct.cells.bin:=n.cells.bin/n.cells,by=.(sample_hto,bin)]

pseudo_bins_[,n.cells.below.thr:=NA]
pseudo_bins_[,n.cells.below.thr:=as.numeric(n.cells.below.thr)]

for(binn in 0:n_bins){
  pseudo_bins_[bin<=binn,n.cells.below.thr:=ifelse(is.na(n.cells.below.thr),sum(n.cells.bin),n.cells.below.thr),by="sample_hto"]

  }

pseudo_bins_[,pct.cells.below.thr:=n.cells.below.thr/n.cells]
ggplot(pseudo_bins_,aes(x=factor(bin),y=pct.cells.below.thr,fill=group))+
  geom_boxplot()+
facet_wrap("hto")
ggsave(fp(out,"boxplot_cumulated_pseudotime_by_group.pdf"))

#visuellement le pique de l'influence semble être à bin = 7 wich corresponf to HSC cells
#so, is there less cells with pseudotim <7 in stimulated LGA compared to stimulated control ?

wilcox.test(pseudo_bins_[hto==T&bin==7&group=="lga"]$pct.cells.below.thr,pseudo_bins_[hto==T&bin==7&group=="ctrl"]$pct.cells.below.thr)
#p-value = 0.04113

wilcox.test(pseudo_bins_[hto==F&bin==7&group=="lga"]$pct.cells.below.thr,pseudo_bins_[hto==F&bin==7&group=="ctrl"]$pct.cells.below.thr)
# p-value = 0.9307

#for each pseudotime
pseudo_bins_[,pval:=wilcox.test(pct.cells.below.thr[group=="lga"],pct.cells.below.thr[group=="ctrl"])$p.value,by=.(hto,bin)]

unique(pseudo_bins_[,.(bin,hto,pval)])[order(hto,bin)]


