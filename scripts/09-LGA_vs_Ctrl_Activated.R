library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("scripts/utils/new_utils.R")
out<-"outputs/09-LGA_vs_Ctrl_Activated"
dir.create(out)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps_h<-subset(cbps,hto==T)


#Distribution change
cbps_h$lineage_hmap<-factor(cbps_h$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))


mtd<-data.table(cbps_h@meta.data,keep.rownames="bc")
mtd[,n.sample:=.N,by="sample_hto"]
mtd[,pct.lin:=.N/n.sample,by=c("sample_hto","lineage_hmap")]

ggplot(unique(mtd,by=c("sample","lineage_hmap")))+
  geom_boxplot(aes(x=lineage_hmap,y=pct.lin,fill=group))


#expr change in lineage
#pseudobulk
Idents(cbps_h)<-"lineage_hmap"
res_lin<-Reduce(rbind,lapply(levels(cbps_h),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_h,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                                 design = ~ group+orig.ident+sex)
  
  dds <- DESeq(dds)
  
  
  mod_mat <- model.matrix(design(dds), colData(dds))
  lga <- colMeans(mod_mat[dds$group == "lga", ])
  ctrl <- colMeans(mod_mat[dds$group == "ctrl", ])


  res <- results(dds,contrast = lga-ctrl,alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

fwrite(res_lin,fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"))

#volcano by lin
res_lin<-fread(fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"))

genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC","PIK3R1","MT-ND3")

ggplot(res_lin[lineage%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],aes(x=log2FoldChange,y=-log10(padj),col=padj<0.11&abs(log2FoldChange)>0.6))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(padj<0.05&
                                        abs(log2FoldChange)>0.6&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("lineage")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")
#ggsave("outputs/figures_epi_response/figure2/2D-pseudo_bulk_deseq2_by_lineage_lga_vs_ctrl_activated.pdf")

#sc 
#run 09A
res_lin_act_sc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")
table(res_lin_act_sc[p_val_adj<0.001&abs(avg_logFC)>0.6]$lineage_hmap)

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


#pathways/biological process up and dn in HSC ####

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
res_dt<-res_lin[lineage=="HSC"]

#gseGO
res_dt[,deg_score:=-log10(padj)*log2FoldChange]
#res_dt[,deg_score_rk:=rank(deg_score)]#du plus surexpr au plus downreg

res_dt1<-merge(res_dt,
              data.table(bitr(res_dt$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F))[,gene:=SYMBOL])[!is.na(ENTREZID)]
genelist<-res_dt1[order(-deg_score)]$deg_score
names(genelist)<-res_dt1[order(-deg_score)]$ENTREZID

res_gsea_go<- gseGO(geneList     =genelist , 
                    ont="BP",
                    exponent = 1,
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
saveRDS(res_gsea_go,fp(out,"res_gsea_go_bp.rds"))
dotplot(res_gsea_go,showCategory=50)

gsea_go_dt<-data.table(as.data.frame(res_gsea_go))
gsea_go_dt[,sens:=ifelse(NES>0,"up","dn")]
gsea_go_dt[,n.enriched:=length(tr(core_enrichment)),"ID"]
fwrite(gsea_go_dt,fp(out,"res_gsea_go_bp.csv.gz"))

gsea_go_dt[p.adjust<0.05][order(p.adjust)][1:50]


ggplot(gsea_go_dt[p.adjust<0.05],aes(x=Description,y=-log10(p.adjust)))+
  geom_point(aes(size=n.enriched,col=NES))+
  facet_grid(~sens,scales = "free_x",space = "free_x")+
  scale_x_discrete(guide = guide_axis(angle=80))+
  scale_y_continuous(limits = c(0,))+
  scale_color_gradient2(low = "blue",mid="grey",high = "red")
ggplot(gsea_go_dt[p.adjust<0.05])+geom_point(aes(x=Description))

#with rank degscore
res_dt[log2FoldChange>0,deg_score:=rank(log2FoldChange*-log10(padj))]
res_dt[log2FoldChange<=0,deg_score:=-rank(abs(log2FoldChange)*-log10(padj))]


res_dt1<-merge(res_dt,
              data.table(bitr(res_dt$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F))[,gene:=SYMBOL])[!is.na(ENTREZID)]
genelist<-res_dt1[order(-deg_score)]$deg_score
names(genelist)<-res_dt1[order(-deg_score)]$ENTREZID

res_gsea_go<- gseGO(geneList     =genelist , 
                    ont="BP",
                    exponent = 1,
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
gsea_go_dt<-data.table(as.data.frame(res_gsea_go))
gsea_go_dt[,sens:=ifelse(NES>0,"up","dn")]
gsea_go_dt[,n.enriched:=length(tr(core_enrichment)),"ID"]
gsea_go_dt[p.adjust<0.05]
gsea_go_dt[p.adjust<0.05&NES<0]$Description

saveRDS(res_gsea_go,fp(out,"res_gsea_go_bp_rank.rds"))
fwrite(gsea_go_dt,fp(out,"res_gsea_go_bp_rank.csv.gz"))


#enrich
# -kegg
possible_genes<-rownames(readRDS('outputs/06-integr_singlecell_cbps/cbps_light.rds'))

res_kegg<-enrichKEGG(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1,
                 universe = bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
res_kegg_dt<-data.table(as.data.frame(res_kegg))
res_kegg_dt[p.adjust<0.1]#4 : TNF, MAPK, Foxo
saveRDS(res_kegg,fp(out,"res_kegg.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg.csv"))

#up
res_kegg_up<-enrichKEGG(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1,
                 universe = bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
res_kegg_dt<-data.table(as.data.frame(res_kegg_up))
res_kegg_dt[p.adjust<0.1]#1 : Phosphatidylinositol signaling system 
saveRDS(res_kegg_up,fp(out,"res_kegg_up.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg_up.csv"))

#dn
res_kegg_dn<-enrichKEGG(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
res_kegg_dt<-data.table(as.data.frame(res_kegg_dn))
res_kegg_dt[p.adjust<0.1]#12: TNF, MAPK
saveRDS(res_kegg_dn,fp(out,"res_kegg_dn.rds"))
fwrite(res_kegg_dt,fp(out,"res_kegg_dn.csv"))


#gseKegg
fcs<-res_dt$log2FoldChange*-log10(res_dt$pvalue)
names(fcs)<-possible_genes
fcs<-sort(fcs,decreasing = T)
head(fcs)

genes.df<-bitr(names(fcs),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
head(genes.df)

fcs<-fcs[genes.df$SYMBOL]
names(fcs)<-genes.df$ENTREZID
head(fcs)
fcs_rank<-fcs
fcs_rank[fcs>0]<-rank(fcs[fcs>0])
fcs_rank[fcs<=0]<-(-rank(abs(fcs[fcs<=0])))
head(fcs_rank)
tail(fcs_rank)
res_gsek<-gseKEGG(geneList     = fcs,
                        exponent=0,
                                    eps=0,
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 1,
                        verbose = FALSE)

saveRDS(res_gsek,fp(out,"res_hto_signature_gsea_kegg.rds"))
res_kegg_dt<-data.table(as.data.frame(res_gsek))
res_kegg_dt[p.adjust<0.2]#no
fwrite(res_kegg_dt,fp(out,"res_hto_signature_kegg.csv"))



res_gsek<-gseGO(geneList     = fcs_rank,
                        exponent=1,
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 0.1,
                        verbose = FALSE)

saveRDS(res_gsek,fp(out,"res_hto_signature_gsea_kegg.rds"))
res_kegg_dt<-data.table(as.data.frame(res_gsek))
res_kegg_dt#yes !


# -go

#molecular function
res_go_mf<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_mf))
res_go_dt[p.adjust<0.1]#14:HSP, methyltrafe
saveRDS(res_go_mf,fp(out,"res_go_mf.rds"))
fwrite(res_go_dt,fp(out,"res_go_mf.csv.gz"))

#biological process
res_go_bp<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_bp))
res_go_dt[p.adjust<0.1]#47
saveRDS(res_go_bp,fp(out,"res_go_bp.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp.csv.gz"))

#up
res_go_mf_up<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_mf_up))
res_go_dt[p.adjust<0.1]#0
saveRDS(res_go_mf_up,fp(out,"res_go_mf_up.rds"))
fwrite(res_go_dt,fp(out,"res_go_mf_up.csv.gz"))

res_go_bp_up<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_bp_up))
res_go_dt[p.adjust<0.1]#0
saveRDS(res_go_bp_up,fp(out,"res_go_bp_up.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp_up.csv.gz"))

#dn
res_go_mf_dn<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_mf_dn))
res_go_dt[p.adjust<0.1]#13
saveRDS(res_go_mf_dn,fp(out,"res_go_mf_dn.rds"))
fwrite(res_go_dt,fp(out,"res_go_mf_dn.csv.gz"))

res_go_bp_dn<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_bp_dn))
res_go_dt[p.adjust<0.1]#80
saveRDS(res_go_bp_dn,fp(out,"res_go_bp_dn.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp_dn.csv.gz"))


res_gsego_bp_dn<-gseGO(geneList     = fcs_rank, 
                       ont="BP",
                        minGSSize    = 50,
                        pvalueCutoff = 0.05,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)

res_go_dt<-data.table(as.data.frame(res_gsego_bp_dn))
res_go_dt[p.adjust<0.1]#



