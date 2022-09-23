out<-'outputs/26-thesis_additional_data'
dir.create(out)
source("scripts/utils/new_utils.R")


#0) Quiescence bias ?
mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd$Phase<-factor(mtd$Phase,levels = c("G1","S","G2M"))

#G2M score
ggplot(mtd)+geom_boxplot(aes(x=hto,y=G2M.Score,fill=group),outlier.shape = NA)+
  coord_cartesian(ylim = c(-0.2,0.12))


#pct Phase HSC/MPP
mtdf<-subset(mtd,lineage_hmap%in%c("HSC","MPP/LMPP"))
mtdf[,n.sample:=.N,.(sample_hto)]
mtdf[,pct.phase:=.N/n.sample,by=.(sample_hto,Phase)]
mts<-unique(mtdf,by=c("sample_hto","Phase"))

ggplot(mts)+geom_boxplot(aes(x=Phase,y=pct.phase*100,fill=group))+labs( x="Phase", y = "HSC/MPP (%)",fill="group")+
  facet_wrap("hto")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))


#bof
#I)Validate differentiation bias with other approach

#RNA velocity####
mtdff<-fread(fp(out,"pseudo_bias_rna_velo_based_lineages.csv.gz"))

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

#cell level
ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))+
  scale_fill_manual(values = c("white","grey"))+
  theme_minimal()
ggsave(fp(out,"boxplot_pseudobias_lga_ctrl.pdf"))

#test
mtdff[,p_val:=wilcox.test(pseudo_bias[group=="lga"],pseudo_bias[group=="ctrl"])$p.value,by="lineage_hmap"]
unique(mtdff[,.(lineage_hmap,p_val)])
#    lineage_hmap        p_val
# 1:     MPP/LMPP 2.433119e-01
# 2:          HSC 9.583314e-12 *
# 3:  Erythro-Mas 7.495362e-06
# 4:     Lymphoid 5.262100e-03
# 5:      Myeloid 8.889535e-08
# 6:       LT-HSC 7.779419e-02
# 7:         <NA> 4.746781e-01


#sample level
mtdffl<-mtdff[lineage_hmap%in%lins]
mtdffl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdffl[,avg_pseudobias:=mean(pseudo_bias),by=.(lineage_hmap,sample)]

mts<-unique(mtdffl,by=c("lineage_hmap","sample"))
ggplot(mts[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))+
  scale_fill_manual(values = c("white","grey"))+
  theme_minimal()
ggsave(fp(out,"boxplot_pseudobias_sample_lvl_lga_ctrl.pdf"))

#test
mts[,p_val:=wilcox.test(pseudo_bias[group=="lga"],pseudo_bias[group=="ctrl"])$p.value,by="lineage_hmap"]
unique(mts[,.(lineage_hmap,p_val)])
#    lineage_hmap      p_val
# 1:     MPP/LMPP 0.28238428
# 2:          HSC 0.05927406 .
# 3:  Erythro-Mas 0.49084249
# 4:     Lymphoid 0.85181485
# 5:      Myeloid 0.57276057
# 6:       LT-HSC 0.31515152



#HTO based stimulation ####
#from 08-HTO_signature get figures:
#- cellcycle_activation

#degs hto vs not by samples
#all lin
library(Seurat)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
cbpsf<-subset(cbps,sample%in%c("ctrlM555","ctrlM518","ctrlM537"))
table(cbpsf$hto)
cbpsf
DefaultAssay(cbpsf)<-"SCT"
#see cst2

#II)Validate KLF2, EGR1, KLF4 regulatory network upon stimulation
#siRNAs####
sicbl<-readRDS("outputs/25-siRNAs_KLF2_cbl/sicbls.rds")
sicbl #  3649 cells

DimPlot(sicbl, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) 
ggsave(fp(out,"umap_sicbls.pdf"))

#degs in MPP/LMPP ?
Idents(sicbl)<-"predicted.lineage"
res<-FindMarkers(sicbl,group.by = "sirna",ident.1 = "KLF2",ident.2 = "CTRL",subset.ident = "MPP/LMPP",logfc.threshold = 0)
res<-data.table(res,keep.rownames="gene")
fwrite(res,fp(out,"res_wilcox_sirna_KLF2_vs_ctrl_MppLmpp.csv"))

res[p_val_adj<0.001&abs(avg_log2FC)>0.25]
ggplot(res,aes(x=avg_log2FC,y=-log10(p_val),col=p_val_adj<0.001&abs(avg_log2FC)>0.25))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-4.5,4.5))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4A-volcano_HSC.pdf"),height = 6.6)

#res gost
resgo<-fread("outputs/25-siRNAs_KLF2_cbl/res_gost_degs_siklf2_vs_ctrl_mpplmpp.csv.gz")
resgo20<-head(resgo[p_value<0.01][str_detect(source,"GO")],20)

ggplot(resgo20)+geom_col(aes(x=-log10(p_value),y=term_name,fill=precision))

resgo20small<-head(resgo[p_value<0.01][str_detect(source,"GO")][term_size<500],20)


ggplot(resgo20small)+geom_col(aes(x=-log10(p_value),y=term_name,fill=precision))


#diff within MPP/LMPP ?
ggplot(sicbl@meta.data)+geom_bar(aes(x=sirna,fill=predicted.cell_type),position = "dodge")
ggplot(sicbl@meta.data)+geom_bar(aes(x=sirna,fill=predicted.cell_type),position = "fill")

mtd<-data.table(sicbl@meta.data,keep.rownames = "bc")
mtd[,n.cells:=.N,by="sirna"]

mtd[,pct.ct:=.N/n.cells,by=c("sirna","predicted.cell_type")]

ct_order<-c("MPP","MPP-Ery","EMP","EMP-cycle","ErP-cycle","Mk/Er","LMPP","GMP","GMP-cycle","DC","CLP","T cell")
mtd[,cell_type:=factor(predicted.cell_type,levels = ct_order)]

mts<-unique(mtd[cell_type%in%ct_order],by=c("sirna","predicted.cell_type"))
ggplot(mts)+geom_col(aes(x=cell_type,y=pct.ct,fill=sirna),position = "dodge")

ct_order2<-c("MPP","MPP-Ery","LMPP")

mtsf<-unique(mtd[cell_type%in%ct_order2],by=c("sirna","predicted.cell_type"))
ggplot(mtsf)+geom_col(aes(x=cell_type,y=pct.ct,fill=sirna),position = "dodge",width = 0.6)

#cell cycle regulation ?

mtd<-data.table(sicbl@meta.data,keep.rownames = "bc")
mtd$Phase<-factor(mtd$Phase,levels = c("G1","S","G2M"))

#G2M score
ggplot(mtd)+geom_boxplot(aes(x=sirna,y=G2M.Score),outlier.shape = NA)+
  coord_cartesian(ylim = c(-0.2,0.12))


#pct Phase
p1<-ggplot(mts)+geom_col(aes(x=Phase,y=pct.phase*100,fill=sirna),position = "dodge")
p1+labs( x="Phase", y = "HSPC (%)",fill="siRNA")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

#degs MPP, to reduce cell type bias
Idents(sicbl)<-"predicted.cell_type"
res<-FindMarkers(sicbl,group.by = "sirna",ident.1 = "KLF2",ident.2 = "CTRL",subset.ident = "MPP",logfc.threshold = 0)
res<-data.table(res,keep.rownames="gene")
fwrite(res,fp(out,"res_wilcox_sirna_KLF2_vs_ctrl_Mpp.csv.gz"))

res[p_val_adj<0.001&abs(avg_log2FC)>0.25] #101
ggplot(res,aes(x=avg_log2FC,y=-log10(p_val),col=p_val_adj<0.001&abs(avg_log2FC)>0.25))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-4.5,4.5))+
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(fp(out,"4A-volcano_MPP.pdf"),height = 6.6)

res[p_val_adj<0.001&avg_log2FC>0.25] #101

#res go
resgo<-fread("outputs/25-siRNAs_KLF2_cbl/res_gost_degs_siklf2_vs_ctrl_mpp.csv.gz")
resgo20<-head(resgo[p_value<0.01][str_detect(source,"GO")],20)

ggplot(resgo20)+geom_col(aes(x=-log10(p_value),y=term_name,fill=precision))

resgo20small<-head(resgo[p_value<0.01][str_detect(source,"GO")][term_size<500],20)


ggplot(resgo20small)+geom_col(aes(x=-log10(p_value),y=term_name,fill=precision))


#overlap KLF2 regulons ?####
regulons<-fread("outputs/10-SCENIC/regulons_extended.csv")
res<-fread(fp(out,"res_wilcox_sirna_KLF2_vs_ctrl_Mpp.csv.gz"))

intersect(regulons[tf=="KLF2"]$gene,
          res[p_val_adj<0.001]$gene)

intersect(regulons[tf=="KLF2"]$gene,
          res[p_val_adj<0.05]$gene)

#3/265
OR(set1 = regulons[tf=="KLF2"]$gene,
   set2=res[p_val_adj<0.05&abs(avg_log2FC)>0.25]$gene,
   size_universe = nrow(res))


res<-fread(fp(out,"res_wilcox_sirna_KLF2_vs_ctrl_MppLmpp.csv"))

intersect(regulons[tf=="KLF2e"]$gene,
          res[p_val_adj<0.001]$gene)
# 
#  [1] "ARF6"     "ARRDC2"   "CALM1"    "DDX3X"    "DNAJA1"   "HLA-C"    "HNRNPUL1"
#  [8] "HSPA1A"   "ITM2B"    "MYADM"    "PLCG2"    "PPP1R10"  "RAB21"    "RBM39"   
# [15] "SRSF5"    "SRSF7"    "UBC"      "ZFP36L2"  "EIF4A2"  
#19/109 19/2012

OR(set1 = regulons[tf=="KLF2"]$gene,
   set2=res[p_val_adj<0.001&abs(avg_log2FC)>0.25]$gene,
   size_universe = nrow(res))

#not signif.. but expected because not in HSC


#
