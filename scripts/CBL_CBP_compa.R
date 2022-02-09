#CBL == CBPs ?
library("Seurat")
source("scripts/utils/new_utils.R")
out<-"outputs/CBL_CBP_compa"
dir.create(out)
#basalement

#ct proport
cbls<-readRDS("../singlecell/outputs/cbl_hto_stim/cbls_integrated_and_hmap.rds")
DimPlot(cbls,group.by = "lineage_hmap",label=T,reduction = "ref.umap")

mtd_cbls<-data.table(cbls@meta.data,keep.rownames = "bc")
mtd_cbls[,hto:=orig.ident!="cbl_ctrl"]
mtd_cbps<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtds<-rbind(mtd_cbps[,source:="CBPs"],mtd_cbls[,source:="CBLs"],fill=T)
ggplot(mtds[hto==F])+geom_bar(aes(x=source,fill=lineage_hmap),position = 'fill')+facet_wrap("hto")
#++ myeloid in cbls basal or other sous pop ?
DimPlot(subset(cbls,orig.ident=="cbl_ctrl"),group.by = "lineage_hmap",label=T)
FeaturePlot(cbls,group.by = "lineage_hmap",label=T,reduction = "umap")
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps.rds")
cbps_all<-merge(cbps,cbls,merge.dr = c("ref.pca","ref.umap"))

cbps_all$source<-ifelse(str_detect(cbps_all$orig.ident,"cbl_"),"CBLs","CBPs")
cbps_all$hto[str_detect(cbps_all$orig.ident,"cbl_")]<-cbps_all[,str_detect(cbps_all$orig.ident,"cbl_")]$orig.ident!="cbl_ctrl"

VlnPlot(subset(cbps_all,hto==F),"MPO",group.by = "lineage_hmap",split.by = "source")
DimPlot(subset(cbls,orig.ident=="cbl_ctrl"),group.by = "lineage_hmap",label=T,reduction = "ref.umap")

#=> no MPO expr in CBLs => not same sous pop
DimPlot(subset(cbls,orig.ident=="cbl_ctrl"),group.by = "lineage_hmap",label=T,reduction = "integrated.umap")
cbls_myelo<-subset(cbls,lineage_hmap=="Myeloid")
DimPlot(subset(cbls_myelo,orig.ident=="cbl_ctrl"),group.by = "cell_type_hmap",label=T,reduction = "integrated.umap")
FeaturePlot(subset(cbls_myelo,orig.ident=="cbl_ctrl"),features = "MPO",reduction = "integrated.umap")

#degs within ct
cbps_b<-subset(cbps_all,hto==F)
#need rm ct with less than 10 cells by condition
table(cbps_b$lineage_hmap,cbps_b$source)

  #             CBLs  CBPs
  # 18             53   163
  # B cell        266   671
  # DC           2055   120
  # Erythro-Mas  1554  3677
  # HSC          4226  7445
  # LT-HSC        783   195
  # Lymphoid      618  2181
  # Mk/Er          57    65
  # MPP/LMPP     6995 19006
  # Myeloid       739  1680
  # T cell        201   272
#ok
Idents(cbps_b)<-"lineage_hmap"
degs<-Reduce(rbind,lapply(levels(cbps_b),function(ct)data.table(FindMarkers(cbps_b,
                                                                                subset.ident=ct,
                                                                                group.by = "source",
                                                                                ident.1 = "CBLs",
                                                                                ident.2 ="CBPs",
                                                                                test.use = "wilcox",
                                                                                logfc.threshold = 0.25,
                                                                                min.pct = 0.1,
                                                                                assay="SCT"
                                                                                ),
                                                                    keep.rownames = "gene")[,cell_type:=ct]))


fwrite(degs,fp(out,"degs_cbls_vs_cbps_basal_by_ct_sct_wilcox.csv.gz"))
degs[p_val_adj<0.05&cell_type=="HSC"]
ggplot(degs[p_val_adj<0.05&abs(avg_log2FC)>0.5])+geom_bar(aes(x=cell_type))

ggplot(degs)+geom_point(aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.05&abs(avg_log2FC)>0.5))+facet_wrap("cell_type")

#rep au stress HTOstim###
#cbls hto proto vs ctrl at 37°C
cbls$hto<-cbls$orig.ident!="cbl_ctrl"
Idents(cbls)<-"lineage_hmap"
table(degs_hto[])
degs_hto<-Reduce(rbind,lapply(levels(cbls),function(ct)data.table(FindMarkers(cbls,
                                                                                subset.ident=ct,
                                                                                group.by = "hto",
                                                                                ident.1 = "TRUE",
                                                                                ident.2 ="FALSE",
                                                                                test.use = "wilcox",
                                                                                logfc.threshold = 0.25,
                                                                                min.pct = 0.1,
                                                                                assay="SCT"
                                                                                ),
                                                                    keep.rownames = "gene")[,cell_type:=ct]))



fwrite(degs_hto,fp(out,"degs_cbls_hto_proto_vs_37degre_wilcox_sct.csv.gz"))
table(degs_hto[p_val_adj<0.05&(avg_log2FC)>0.5]$cell_type)
   #      18      B cell          DC Erythro-Mas         HSC      LT-HSC    Lymphoid       Mk/Er 
   #        6         194         300         107         202         283         343          12 
   # MPP/LMPP     Myeloid      T cell 
   #      240         312         163 

genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC","PIK3R1","MT-ND3")


ggplot(degs_hto[cell_type%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.05&abs(avg_log2FC)>0.25))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(p_val_adj<0.05&
                                        abs(avg_log2FC)>0.5&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("cell_type")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")


#cbls hto proto with Ab vs hto proto without Ab
cbls$hto<-cbls$orig.ident!="cbl_ctrl"
Idents(cbls)<-"lineage_hmap"
degs_hto<-Reduce(rbind,lapply(levels(cbls),function(ct)data.table(FindMarkers(cbls,
                                                                                subset.ident=ct,
                                                                                group.by = "orig.ident",
                                                                                ident.1 = "cbl_hto",
                                                                                ident.2 ="cbl_stmd",
                                                                                test.use = "wilcox",
                                                                                logfc.threshold = 0.25,
                                                                                min.pct = 0.1,
                                                                                assay="SCT"
                                                                                ),
                                                                    keep.rownames = "gene")[,cell_type:=ct]))


fwrite(degs_hto,fp(out,"degs_cbls_hto_proto_with_Ab_vs_without.csv.gz"))

#
table(degs_hto[p_val_adj<0.05&abs(avg_log2FC)>0.5]$cell_type)
     # B cell          DC Erythro-Mas         HSC      LT-HSC    Lymphoid    MPP/LMPP     Myeloid 
     #      8          15           8          11          17          10          11          10 
     # T cell 
     #      5 

genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC","PIK3R1","MT-ND3")


ggplot(degs_hto[cell_type%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")],aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.05&abs(avg_log2FC)>0.25))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(p_val_adj<0.05&
                                        abs(avg_log2FC)>0.5&gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 5000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  facet_wrap("cell_type")+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#cbls hto proto vs cbls directly put in 10x after frozen
cbls_cold<-subset(cbls,orig.ident=="cbl_stmd")
cbls_ctrl<-readRDS("../singlecell/outputs/")

#cbps vs cbls stim
cbps_hto<-subset(cbps_all,hto==T)
#ct prop
ggplot(mtds[hto==T])+geom_bar(aes(x=source,fill=lineage_hmap),position = 'fill')+facet_wrap("hto")
#++ DCs in cbls basal or other sous pop ?
DimPlot(subset(cbps_hto,source=="CBLs"),group.by = "lineage_hmap",label=T,reduction = "ref.umap")
FeaturePlot(cbps_hto,group.by = "lineage_hmap",label=T,reduction = "umap")
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps.rds")
cbps_all<-merge(cbps,cbls)
cbps_all$source<-ifelse(str_detect(cbps_all$orig.ident,"cbl_"),"CBLs","CBPs")
cbps_all$hto[is.na(cbps_all$hto)]<-cbps_all[,is.na(cbps_all$hto)]$orig.ident=="cbl_ctrl"

VlnPlot(subset(cbps_all,hto==F),"MPO",group.by = "lineage_hmap",split.by = "source")
DimPlot(subset(cbls,orig.ident=="cbl_ctrl"),group.by = "lineage_hmap",label=T,reduction = "ref.umap")

#=> based on MPO markers, same sous pop

#degs within ct
cbps_b<-subset(cbps_all,hto==F)
#need rm ct with less than 10 cells by condition
table(cbps_b$lineage_hmap,cbps_b$source)

  #             CBLs  CBPs
  # 18             53   163
  # B cell        266   671
  # DC           2055   120
  # Erythro-Mas  1554  3677
  # HSC          4226  7445
  # LT-HSC        783   195
  # Lymphoid      618  2181
  # Mk/Er          57    65
  # MPP/LMPP     6995 19006
  # Myeloid       739  1680
  # T cell        201   272
#ok
Idents(cbps_b)<-"lineage_hmap"
degs<-Reduce(rbind,lapply(levels(cbps_b),function(ct)data.table(FindMarkers(cbps_b,
                                                                                subset.ident=ct,
                                                                                group.by = "source",
                                                                                ident.1 = "CBLs",
                                                                                ident.2 ="CBPs",
                                                                                test.use = "wilcox",
                                                                                logfc.threshold = 0.25,
                                                                                min.pct = 0.1,
                                                                                assay="SCT"
                                                                                ),
                                                                    keep.rownames = "gene")[,cell_type:=ct]))


fwrite(degs,fp(out,"degs_cbls_vs_cbps_basal_by_ct_sct_wilcox.csv.gz"))
degs[p_val_adj<0.05]
