#figures these
out<-"outputs/28-figures_thesis_complementary_works"
dir.create(out)
library(Seurat)
source("scripts/utils/new_utils.R")

#I)	Is there a stress response in our samples ?###
# A) RNA velocity analysis####
#fig1A:velocity stream
#ok
#fig1B: genes most actively expressed 
#ok
#(bigger RNA velocity and variance across HSPCs) 
#=genes related to control of activation like FOS, DUSP1, FOSB, ZFP36, HES1, and NFKBIA 

#fig1C:regulons enrichment in the main contributors of this transcriptional dynamics model
#=EGR1 regulon being the most enriched, while KLF2 and JUN are in the top5 

#fig1D: Diff bias
#LGA HSCs differentiate more  toward more differentiated cells compared to CTRL HSC 
#but not others lin


# B)	 rapidly processed cells.###
#   1)	HTO protocol effect####
#fig2A: strong sample wide gene expression difference between the 2 protocols 
res_hto<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")
res_hto[padj<0.05&abs(log2FoldChange)>0.5]
res_hto[padj<0.05&log2FoldChange>0.5] #1075/1518
res_hto_gobp_up<-fread("outputs/08-HTO_signature/res_hto_signature_go_bp_up.csv")
res_hto_gobp_up[order(pvalue)][1:20]

#fig2B: KEGG/GO enrichemnt in HTO upreg
#fig2C: subpop DE analysis
#fig2D: KEG/GO enrichment
#fig2E: regulons enrichement in DEGS
#fig2F: cells prepared with HTO are more committed in the mitosis process (G2/M phase) 
#than non HTO-prepared cells

#   2) LGA vs CTRL####
#fig3A: any significant gene expression change across lineage between LGA and CTRL samples 
#fig3B: no significant subpopulation proportion change between LGA and CTRL


#II)	Can we validate the EGR1/ KLF2/KLF4 regulatory network ?###  
# A)	siRNAs validation####
#fig4A: siKLF2  reduce expression of KLF2 through RT-qPCR
#fig4B: UMAP + cell type markers expression
sicbl<-readRDS("outputs/25-siRNAs_KLF2_cbl/sicbls2.rds")
DimPlot(sicbl)

key_genes_lin<-c("ID1","DUSP2", #LT-HSC
           "AVP","FOS", #HSC
           "MLLT3","CDK6", #MPP
           "CD99", #LMPP
           "LTB", #CLP
           "TNFAIP3","CD7", #T cell
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "PLEK","HBD", #Mk/Er
           "MPO","CEBPA","CTSG","AZU1") #GMP
DotPlot(sicbl,features = key_genes_lin,cols = c("white","black"))+coord_flip()

#fig4C: 4% of cells was annotated as HSC while 24% of CD34+ cells  was HSC in our previous unincubated samples
mtdsi<-data.table(sicbl@meta.data)
mtdsi[,pct.ct:=.N/nrow(mtdsi)*100,by="cell_type2"]

#sicbl$lineage_hmap<-ifelse(sicbl$cell_type2%in%c("MPP","LMPP", "MPP"),"MPP/LMPP",)
p1<-ggplot(unique(mtdsi,by="cell_type2"))+geom_col(aes(x=cell_type2,y=pct.ct,fill=cell_type2))+
  labs(,y="Proportion (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45))

mtdp<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdp<-mtdp[(hto)]
mtdp[,pct.ct:=.N/nrow(mtdp)*100,by="lineage_hmap"]
mtdps<-unique(mtdp[lineage_hmap%in%c('HSC',"MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")],by="lineage_hmap")
mtdps[,lineage_hmap:=factor(lineage_hmap,levels=c('HSC',"MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")),by="lineage_hmap"]

p2<-ggplot(mtdps)+
  geom_col(aes(x=lineage_hmap,y=pct.ct,fill=lineage_hmap))+
  labs(y="Proportion (%)")+
  scale_x_discrete(guide = guide_axis(angle = 45),limits=c('HSC',"MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid"))
p1+p2

#fig4D: slight reduction of HSC and MPP cells in siKLF2 compared to siCTRL
mtdsi<-data.table(sicbls@meta.data)
mtdsi[,pct.hscmpp:=sum(cell_type2%in%c("HSC","MPP"))/.N,by="sirna"]

ggplot(unique(mtdsi,by="sirna"))+geom_col(aes(x=sirna,y=pct.hscmpp,fill=sirna))+
  scale_fill_manual(values = c("black","grey"))

#chiSquared test
chisq.test(c(sum(mtdsi$sirna=="KLF2"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="KLF2"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   p = c(sum(mtdsi$sirna=="CTRL"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="CTRL"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   rescale.p = T)$p.value #0.00035

chisq.test(c(sum(mtdsi$sirna=="KLF2"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="KLF2"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   y = c(sum(mtdsi$sirna=="CTRL"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="CTRL"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   )$p.value #0.00035



#fig4E: KLF2 regulons are specifically active in HSC 
sicbl<-readRDS("outputs/25-siRNAs_KLF2_cbl/sicbls2_with_regulons_activity.rds")
FeaturePlot(sicbl,"KLF2",min.cutoff = "q10",max.cutoff = "q95",label=T,repel=F)


#fig4F: DE regulons KLF2 between siKLF2 and siCTRL 
#few DEGs p-value threshold (n=4), but 15/89 (16%) genes at nominal pvalue 

#fig4G: 4 regulons significantly enriched in downregulated genes (adjusted p-value <0.1),
#with KLF2 regulon being the second most enriched regulon behind STAT3 

# B)	Single cell multimodal validation ####
cbls<-readRDS("outputs/27-multiomics_cbl_stim/cbls_stim_final.rds")
cbls #4354
table(cbls$condition) #2001    2353

table(cbls$group) 
  # CTRL  Doublet      LGA Negative 
  #   1835       70     1667      753 

#fig5A: UMAP + cell type markers expression

#fig5B: net decrease of HSC in the stimulated conditions compared to unstimulated conditions 
#accompagnied by a net increase of MPP/LMPP cells as well as DC progenitors
mtds<-data.table(cbls@meta.data,keep.rownames = "bc")
mtds[,n.sample:=.N,by="condition"]

mtds[,pct.ct:=.N/n.sample,by=c("condition","ident.2")]

lvls<-c("HSC","HSC-Mk","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas","DC Prog")
mtdsf<-mtds[ident.2%in%lvls]

mtdsf[,lineage:=factor(ident.2,levels = lvls)]
ggplot(unique(mtdsf,by=c("condition","lineage")))+
  geom_col(aes(x=lineage,y=pct.ct,fill=condition),position = position_dodge(width=0.9),width = 0.8,)+
  scale_fill_manual(values = c("black","grey"))+
  theme_minimal()


#fig5C: stimulated to unstimulated HSC, including 239 upregulated and 342 downregulated genes
DefaultAssay(cbls)<-"SCT"
Idents(cbls)
res<-data.table(FindMarkers(cbls,group.by="condition",ident.1="stim_2h",logfc.threshold = 0,subset.ident = "HSC"),keep.rownames = "gene")

res<-fread("outputs/27-multiomics_cbl_stim/res_hsc_stim_vs_no_stim.csv.gz")
down_lga<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][padj<0.05&log2FoldChange<(-0.5)]
up_lga<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"][padj<0.05&log2FoldChange>0.5]

res[,lga_degs_cat:=ifelse(p_val_adj<0.05&abs(avg_log2FC)>0.25&gene%in%c(down_lga$gene,up_lga$gene),ifelse(gene%in%down_lga$gene,"down_in_lga","up_in_lga"),"")]

ggplot(res,aes(x=avg_log2FC,y=-log10(p_val),col=p_val_adj<0.05&abs(avg_log2FC)>0.25))+
  geom_point()+
  geom_label_repel(aes(label=lga_degs_cat),max.overlaps = 2000)+
  scale_color_manual(values = c("grey","red"))+theme_minimal()+theme(legend.position = "bottom")


ggplot(res,aes(x=avg_log2FC,y=-log10(p_val),col=p_val_adj<0.05&abs(avg_log2FC)>0.25))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(p_val_adj<0.05&abs(avg_log2FC)>0.25&gene%in%down_lga$gene,gene,"")),max.overlaps = 2000)+
  scale_color_manual(values = c("grey","red"))+theme_minimal()+theme(legend.position = "bottom")


#fig5D: upregulated genes are strongly enriched for regulons of JUN, FOSB, ARID5A, FOS, and JUNB,
#being the top5 most enriched regulons;
res<-fread("outputs/27-multiomics_cbl_stim/res_regulons_enrichment_in_upregulated_genes_stim_hsc2h.csv.gz")


ggplot(res[padj<0.05][order(padj)])+
  geom_point(aes(x=term,y=-log10(pval),size=n.overlap,col=fold.enrichment))+
  coord_cartesian(ylim = c(0,10))+
  scale_color_gradient2(low = "grey",mid = "orange",high = "red",midpoint = 2)+
  scale_x_discrete(limits=c(res[padj<0.05][order(-pval)]$term),guide = guide_axis(angle=60))+theme_minimal()

#while EGR1 and KLF2 regulons are also found significantly enriched
#fig5F:ZFP36L2 from this network having a lot of peaks (26) associated with it suggesting tight regulation
DefaultAssay(cbls)<-"peaks"
links<-fread("outputs/27-multiomics_cbl_stim/res_peak_genes_all_genes_regulons_linkage.csv")
links[,score:=as.numeric(score)]
links[gene=="ZFP36L2"]
links<-links[as.numeric(score)]
library(GenomicRanges)
linksgr<-makeGRangesFromDataFrame(links,keep.extra.columns = T)
Links(cbls)<-linksgr
p1 <- CoveragePlot(
  object = cbls,
  region = "ZFP36L2",
  features = "ZFP36L2",
  expression.assay = "SCT",
  idents = lvls,
  extend.upstream = 100000,
  extend.downstream = 10000
)
p1
p2 <- CoveragePlot(
  object = cbls,
  region = "LYZ",
  features = "LYZ",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)
#fig5G: 46% of genes of EGR1 regulon have an EGR1 motif on a peak associated to their expression 
#(p-value < 0.05, over-representation test), while 45% for KLF2 (p-value < 0.01, over-representation test)
#and only 38% for KLF4 (ns)

res_enr<-fread("outputs/27-multiomics_cbl_stim/res_enrichment_cres_egrkl_in_regulons.csv")
ggplot(res_enr[term%in%c('EGR1',"KLF2","KLF4")])+
  geom_point(aes(x=term,size=fold.enrichment,y=pct.term.overlap*100,col=pval<0.05))+
  facet_wrap("query")+scale_color_manual(values =  c("grey","red"))+theme_bw()+
  labs(x="Regulon",y="% of the regulon")


#fig5H : DMCs associated CREs are strongly enriched for EGR1/KLF2/KLF4 TF motifs with 91% of 
#these CREs having at least one of this TF motif (p-value < 0.0001
cpgs_in_cres<-fread("outputs/27-multiomics_cbl_stim/cpgs_in_CREs.csv.gz")


res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

cres_meth<-merge(cpgs_in_cres,res_meth,by="cpg_id")


res_motifpeak_genes<-fread("outputs/27-multiomics_cbl_stim/res_peaks_genes_links_EGR1_KLF2_KLF4_motifs_anno.csv")
res_motifpeak_genesf<-res_motifpeak_genes[!is.na(motif)]
tfs<-c("EGR1","KLF2","KLF4")
genes_egrncres<-unique(res_motifpeak_genesf[motif.name%in%tfs][,.(gene,motif.name,peak)])
cres_cpgs_egrn<-merge(cpgs_in_cres,genes_egrncres,by="peak",,allow.cartesian=TRUE)

cres_meth_egrn<-merge(cres_cpgs_egrn,res_meth,by="cpg_id")

#enrichment  for egrn CREs in DMCs asso CREs ?
univ<-unique(cres_meth$peak)
OR(set1 = unique(cres_meth[P.Value<0.001&abs(logFC)>25]$peak),
   set2= unique(cres_meth_egrn$peak),size_universe = length(univ)) #p = 5.710426e-20

intersect(cres_meth_egrn$peak,unique(cres_meth[P.Value<0.001&abs(logFC)>25]$peak))
#1532/2311(66%) peak asso with cpgs have EGR1/klf2/klf4. 198/216(91%) peaks with DMCs have EGR1/klf2/klf4
renv::install("bioc::ggVennDiagram")
library("ggVennDiagram")
x<-list(`CREs containing DMCs` = unique(cres_meth[P.Value<0.001&abs(logFC)>25]$peak),
   `CREs containing EGR1, KLF2, or KLF4 motif`= unique(cres_meth_egrn$peak))
ggVennDiagram(x, label_alpha = 0,label_color = "white" )

