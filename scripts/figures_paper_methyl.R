#figures paper methyl
#see other test on figures_paper_methyl_annexe
source("../methyl/scripts/utils/new_utils.R")
library(limma)
library(pheatmap)
out<-"../methyl/outputs/figures_paper_meth"
dir.create(out)


#validation cohorts
meth<-fread("../methyl/datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
mtd<-fread("../methyl/datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")

var_to_model<-c("Group_Sex","Mat.Age","latino","Group_Complexity_Fac")

mtdf<-mtd[rowSums(is.na(mtd[,.SD,.SDcols=var_to_model]))==0][batch==2]
formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac 

design<-model.matrix(formule,data = data.frame(mtdf,row.names = "sample"))
fit <- lmFit(data.frame(meth,row.names = "locisID")[,mtdf$sample], design)

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

res<-Reduce(rbind,lapply(colnames(cont.matrix), function(comp)data.table(topTable(fit2,coef = comp,n = Inf),keep.rownames = "cpg_id")[,compa:=comp]))
fwrite(res,fp(out,"res_limma_cohort2.tsv.gz"),sep="\t")
res<-fread(fp(out,"res_limma_cohort2.tsv.gz"),sep="\t")

table(res[adj.P.Val<0.2&abs(logFC)>30][,hyper_meth:=logFC>0][,.(hyper_meth,compa)])
p<-ggplot(res[compa%in%c("C.L","CF.LF",'CM.LM',"LM.LF")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=adj.P.Val<0.2&abs(logFC)>30))+
  facet_wrap("compa")+
  scale_color_manual(values = c("grey","red"))
ggsave(fp(out,"volcano_plot_cohort2.png"),plot=p,height=5,width=7)

#comp to cohort 1 
mtdf<-mtd[rowSums(is.na(mtd[,.SD,.SDcols=var_to_model]))==0][batch==1]
formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac 

design<-model.matrix(formule,data = data.frame(mtdf,row.names = "sample"))
fit <- lmFit(data.frame(meth,row.names = "locisID")[,mtdf$sample], design)

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

#merge res cohort 1 and 2
res<-rbind(res[,cohort:=2],
           Reduce(rbind,lapply(colnames(cont.matrix),
                               function(comp)data.table(topTable(fit2,coef = comp,n = Inf),
                                                        keep.rownames = "cpg_id")[,compa:=comp]))[,cohort:=1])
fwrite(res,fp(out,"res_limma_cohorts.tsv.gz"),sep="\t")
#cl
ggplot(res[compa%in%c("C.L")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<10^-4&abs(logFC)>30),size=)+
  facet_wrap("cohort")+
  scale_color_manual(values = c("grey","red"))

#cflf
ggplot(res[compa%in%c("CF.LF")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<10^-4&abs(logFC)>30))+
  facet_wrap("cohort")+
  scale_color_manual(values = c("grey","red"))

#cmlm
ggplot(res[compa%in%c("CM.LM")])+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<10^-4&abs(logFC)>30))+
  facet_wrap("cohort")+
  scale_color_manual(values = c("grey","red"))

table(res[P.Value<10^-4&abs(logFC)>30][,hyper_meth:=logFC>0][,.(hyper_meth,compa,cohort)])




#res meth CL
res_cl<-fread("../methyl/outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")

ggplot(res_cl)+
  geom_point(aes(x=meth.change,y=-log10(pval),col=pval<10^-3&abs(meth.change)>30),size=0.1)+
  scale_color_manual(values = c("grey","red"))


res_cl[pval<10^-3&meth.change>30]
res_cl[,padj:=p.adjust(pval,method = 'BH')]
res_cl[padj<0.05]
resg_cl<-unique(res_cl[order(gene,pval)],by="gene")
resg_cl[GeneScore>60&pval<10^-3]

plot(density(resg_cl$GeneScore))
abline(v=60)
ggplot(resg_cl)+
  geom_point(aes(x=GeneScore,y=-log10(pval),col=pval<10^-3&abs(GeneScore)>60),size=0.1)+
  scale_color_manual(values = c("grey","red"))

res<-fread("../methyl/outputs/model14_without_iugr/2020-10-08_all_res_with_perm.csv")



#OR 

renv::install("bioc::clusterProfiler")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
res_kegg<-enrichKEGG(bitr(resg_cl[GeneScore>60&pval<10^-3]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(res_kegg)


emapplot(pairwise_termsim(res_kegg))
saveRDS(res_kegg,"../methyl/res_kegg_or_cl_temp.rds")

#Integration basal
source("../methyl/scripts/utils/new_utils.R")
library(Seurat)

out<-"../methyl/outputs/figures_paper_meth"

  #Lineage bias
cbps<-readRDS("../singlecell/outputs/02-hematopo_datasets_integration/cbps0-8_clean/cbps0-8_clean.rds")
mtd<-data.table(cbps@meta.data,keep.rownames ="bc")


cbps_cl<-subset(cbps,hto==F&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L")
mtd_cl<-data.table(cbps_cl@meta.data,keep.rownames ="bc")
mtd_cl[,n.sample:=.N,by=c("sample","orig.ident")]
mtd_cl[,pct.ct:=.N/n.sample,by=c("sample","orig.ident","cell_type")]
mtd_cl[,pct.lin:=.N/n.sample,by=c("sample","orig.ident","lineage1")]

lineages<-c("HSC","MPP","Erythroid","Lymphoid","Myeloid","LMPP")


ggplot(unique(mtd_cl[lineage1%in%lineages][,.SD,.SDcols=unique(colnames(mtd_cl))],by=c("sample","orig.ident","lineage1")))+
  geom_boxplot(aes(x=group,y=pct.lin,fill=group))+facet_wrap("lineage1",scales = "free")


  #DEGs
degs_cl<-fread("outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")

ggplot(degs_cl,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6,gene,"")),
                   max.overlaps=3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#by lineage
degs_lin<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_by_lin/res_de_analysis_all_genes.csv")
ggplot(degs_lin,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6,gene,"")),
                   max.overlaps=3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  facet_wrap("lineage")+
  theme_minimal() +
  theme(legend.position = "bottom")

#single cell degs

degs_cl_sc<-fread("outputs/08-DEGs_LGA_no_stress/sc_edger_deseq2_all_cbps/res_de_analysis_all_genes.csv")

degs_cl_sc[p_val_adj<0.01&abs(avg_logFC)>0.6]

ggplot(degs_cl_sc,aes(x=avg_logFC,y=-log10(p_val),col=p_val_adj<0.1&abs(avg_logFC)>0.6))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")



#correl scdegs with methyl

res_cl<-fread("../methyl/outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")

res_cl[gene=="HES1"]

res_cl_merge<-merge(degs_cl_sc,res_cl)
res_cl_merge<-unique(res_cl_merge,by="gene")
ggplot(res_cl_merge) + geom_point(aes(x=-log10(p_val_adj)*avg_logFC,y=GeneScore))
ggplot(res_cl_merge) + geom_point(aes(x=avg_logFC,y=GeneScore))

#Over repre test GO BP and KEGG pathway

res_cl_merge[,gs_scaled:=scale(GeneScore,center = F)]
ggplot(res_cl_merge)+geom_density(aes(x=gs_scaled))
ggplot(res_cl_merge)+geom_density(aes(x=avg_logFC))

res_cl_merge[,correl:=ifelse(gs_scaled>1&avg_logFC>1,"Hypermet-Upreg",
                             ifelse(gs_scaled>1&avg_logFC<(-1),"Hypermet-Downreg",
                              ifelse(gs_scaled<(-1)&avg_logFC<(-1),"Hypomet-Downreg",
                               ifelse(gs_scaled<(-1)&avg_logFC>1,"Hypomet-Upreg","no change"))))]
ggplot(res_cl_merge) + geom_point(aes(x=avg_logFC,y=gs_scaled,col=correl))

#enrichKEGG
library(clusterProfiler)

#downreg
res_kegg<-enrichKEGG(bitr(res_cl_merge[correl=="Hypermet-Downreg"]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.1)

dotplot(res_kegg,showCategory =20)

#upreg
res_kegg_up<-enrichKEGG(bitr(res_cl_merge[correl=="Hypermet-Upreg"]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.1)
as.data.frame(res_kegg_up) #work
dotplot(res_kegg_up,showCategory =20)


#INTEGR WITH STRESS
  #Lineage bias
cbps_clh<-subset(cbps,hto==T&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L")
mtd_clh<-data.table(cbps_clh@meta.data,keep.rownames ="bc")
mtd_clh[,n.sample:=.N,by=c("sample","orig.ident")]
mtd_clh[,pct.ct:=.N/n.sample,by=c("sample","orig.ident","cell_type")]
mtd_clh[,pct.lin:=.N/n.sample,by=c("sample","orig.ident","lineage1")]
ggplot(unique(mtd_clh[lineage1%in%lineages][,.SD,.SDcols=unique(colnames(mtd_clh))],by=c("sample","orig.ident","lineage1")))+
  geom_boxplot(aes(x=group,y=pct.lin,fill=group))+facet_wrap("lineage1",scales = "free")

#correl scdegs with methyl
res_cl<-fread("../methyl/outputs/model14_without_iugr/2020-09-16_res_C.L_with_GeneScore_and_permut.csv")
degs_clh_sc<-fread("../singlecell/outputs/07-DEGs_LGA_stress/sc_edger_deseq2_all_cbps/res_de_analysis_all_genes.csv")
res_clh_merge<-merge(degs_clh_sc,unique(res_cl[order(gene,pval)],by="gene"))

ggplot(res_clh_merge) + geom_point(aes(x=-log10(p_val_adj)*avg_logFC,y=GeneScore))
ggplot(res_clh_merge) + geom_point(aes(x=avg_logFC,y=GeneScore))


#enrichKEGG with clusterprofiler

res_keggh<-enrichKEGG(bitr(res_clh_merge[correl=="Hypermet-Downreg"]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.1)
as.data.frame(res_keggh) #work
dotplot(res_keggh,showCategory =20)

emapplot(pairwise_termsim(res_keggh))

res_keggh_up<-enrichKEGG(bitr(res_clh_merge[correl=="Hypermet-Upreg"]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.1)
as.data.frame(res_keggh_up) #work
dotplot(res_keggh_up,showCategory =20)

#wnt and growth hormone genes
genes_wnt<-tr(as.data.frame(res_keggh)["hsa04310","geneID"],tradEntrezInSymbol = T)
genes_growth<-tr(as.data.frame(res_keggh)["hsa04935","geneID"],tradEntrezInSymbol = T)

res_keggh_dt<-data.table(as.data.frame(res_keggh))
res_keggh_dt[,gene:=paste(tr(geneID,tradEntrezInSymbol = T),collapse = "|"),by="ID"]
res_keggh_dt[str_detect(gene,"HES1")]

#HSC specific
DefaultAssay(cbps_clh)<-"SCT"
Idents(cbps_clh)<-"lineage"
FeaturePlot(cbps_clh,features =genes_of_interest[1:6] , max.cutoff = "q95")
FeaturePlot(cbps_clh,features =genes_of_interest[7:length(genes_of_interest)], max.cutoff = "q95")
FeaturePlot(cbps_clh,features ="SOCS3",label=T ,max.cutoff = "q95")
FeaturePlot(cbps_clh,features ="SOCS1",label=T ,max.cutoff = "q95")
FeaturePlot(cbps_clh,features ="JUN",label=T )
FeaturePlot(cbps_clh,features ="FOS",label=T )
FeaturePlot(cbps_clh,features ="HES1",label=T ,max.cutoff = "q95")

#degs by lineage
#pseudo bulk

#in basal
degs_cl_ps_lin1<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_by_lineage1/res_de_analysis_all_genes.csv")

ggplot(degs_cl_ps_lin1[padj<0.1&abs(log2FoldChange)>0.6])+geom_bar(aes(x=lineage,fill=lineage))

ggplot(degs_cl_ps_lin1)+
  geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.05&abs(log2FoldChange)>0.6),size=0.2)+
               facet_wrap("lineage")+
               scale_color_manual(values = c("grey","red"))

ggplot(degs_cl_ps_lin1,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point(size=0.2)+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%union(genes_wnt,genes_growth),gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))+
  facet_wrap("lineage")


genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",                        "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

p1<-ggplot(degs_clh_ps_lin1[lineage%in%c("HSC")],aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))


p2<-ggplot(degs_clh_ps_lin1[lineage%in%c("MPP")],aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point(si)+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))

p1+p2+plot_layout(guides = "collect")

#in hto
degs_clh_ps_lin<-fread("../singlecell/outputs/07-DEGs_LGA_stress/pseudobulk_deseq2_by_lin/res_de_analysis_all_genes.csv")
ggplot(degs_clh_ps_lin[padj<0.1&abs(log2FoldChange)>0.6])+geom_bar(aes(x=lineage,fill=lineage))

ggplot(degs_clh_ps_lin)+
  geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.05&abs(log2FoldChange)>0.6),size=0.2)+
               facet_wrap("lineage")+
               scale_color_manual(values = c("grey","red"))

ggplot(degs_clh_ps_lin,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point(size=0.2)+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%union(genes_wnt,genes_growth),gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))+
  facet_wrap("lineage")


genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",                        "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

p1<-ggplot(degs_clh_ps_lin[lineage%in%c("HSC")],aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))


p2<-ggplot(degs_clh_ps_lin[lineage%in%c("MPP")],aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))

p1+p2+plot_layout(guides = "collect")

#degs by lineage1
#pseudo bulk
degs_clh_ps_lin1<-fread("../singlecell/outputs/07-DEGs_LGA_stress/pseudobulk_deseq2_by_lineage1/res_de_analysis_all_genes.csv")

table(degs_clh_ps_lin1[padj<0.05]$lineage)

ggplot(degs_clh_ps_lin1[padj<0.1&abs(log2FoldChange)>0.6])+geom_bar(aes(x=lineage,fill=lineage))

ggplot(degs_clh_ps_lin1)+
  geom_point(aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.05&abs(log2FoldChange)>0.6),size=0.2)+
               facet_wrap("lineage")+
               scale_color_manual(values = c("grey","red"))

ggplot(degs_clh_ps_lin1,aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point(size=0.2)+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%union(genes_wnt,genes_growth),gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))+
  facet_wrap("lineage")


genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",                        "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

p1<-ggplot(degs_clh_ps_lin1[lineage%in%c("HSC")],aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))


p2<-ggplot(degs_clh_ps_lin1[lineage%in%c("MPP")],aes(x=log2FoldChange,y=-log10(pvalue),col=padj<0.1&abs(log2FoldChange)>0.6))+
  geom_point(si)+ 
  geom_label_repel(aes(label = ifelse(padj<0.1&
                                        abs(log2FoldChange)>0.6&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps=3000,
                   segment.color = 'grey50')+
               scale_color_manual(values = c("grey","red"))

p1+p2+plot_layout(guides = "collect")
#scDEGs


#SCENIC

cbps<-readRDS("outputs/06-SCENIC/cbps0-8_clean/all_cbps_with_regulons_activity.rds")
tf_lin<-fread("outputs/06-SCENIC/cbps0-8_clean/tf_markers_lineage.csv.gz")
tf_lin[lineage=="Erythroid"]
DefaultAssay(cbps)<-"TF_AUC"
FeaturePlot(cbps,"GATA1",label = T)

tf_lin[lineage=="Myeloid"]$regulon
FeaturePlot(cbps,"CEBPA",label = T)

tf_lin[lineage=="Lymphoid"]
tf_lin[lineage=="Lymphoid"&avg_log2FC>0]
FeaturePlot(cbps,"IRF8",label = T,min.cutoff = 0.2,max.cutoff = 0.4)

tf_lin[lineage=="HSC"&avg_log2FC>0]
tf_lin[lineage=="HSC"&avg_log2FC>0]$regulon
FeaturePlot(cbps,c("EGR1","FOSB","JUN","STAT3"),label = T,min.cutoff = 0.1,max.cutoff = 0.6,pt.size = 0.1
            )

tf_lin[lineage=="MPP"&avg_log2FC>0]
tf_lin[lineage=="MPP"&avg_log2FC>0]$regulon
FeaturePlot(cbps,c("MEIS1e"),label = T,min.cutoff = 0.1,max.cutoff = 0.5,pt.size = 0.2)


#tf_bias_b<-#[to do]
#tf_bias_b[regulon%in%c("STAT3","JUN","MEIS1e","SPI1","EBF")&p_val_adj<0.001]

tf_bias_h<-fread("outputs/06-SCENIC/cbps0-8_clean/regulon_activity_after_stress_lga_vs_ctrl_by_lineage.csv.gz")
tf_bias_h[regulon%in%c("STAT3","JUN","MEIS1e","SPI1","EBF")&p_val_adj<0.001]
