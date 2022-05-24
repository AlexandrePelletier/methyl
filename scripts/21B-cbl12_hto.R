
out<-"outputs/21-CREs_inference_with_sc_multiomics_cbls"

source("scripts/utils/new_utils.R")
library(Seurat)


set.seed(1234)


# load the RNA and HTO data

counts <- Read10X_h5("~/RUN/Run_698_single_cell/Output/cellranger_count/cbl12/outs/raw_feature_bc_matrix.h5")
#counts$`Antibody Capture`
cbl12<-readRDS("outputs/21-CREs_inference_with_sc_multiomics_cbls/cbl12.rds")

hto_mat<-t(as.matrix(counts$`Antibody Capture`))
rownames(hto_mat)<-"HTO-1"
colnames(hto_mat)
cbl12[["HTO"]]<-CreateAssayObject(counts = hto_mat[,colnames(cbl12),drop=FALSE])
DefaultAssay(cbl12)<-"HTO"
FeaturePlot(cbl12,"HTO-1",max.cutoff = "q95")
summary(cbl12@assays$HTO@counts[1,])
 # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
 #    0.00    15.00    20.00    55.13    25.00 84860.00 
cbl12@assays$HTO@counts[1,cbl12@assays$HTO@counts[1,]>100]<-100
summary(cbl12@assays$HTO@counts[1,])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.00   15.00   20.00   21.55   25.00  100.00

VlnPlot(cbl12,"HTO-1",group.by = "orig.ident",pt.size = 0,slot = "counts")

VlnPlot(cbl12,"HTO-1",group.by = "orig.ident",pt.size = 0,slot = "counts")+
  geom_hline(yintercept=33 )



sum(cbl12@assays$HTO@counts[1,]>33) #565/6500 (~8%) cells with High HTO count, expected 14%

VlnPlot(cbl12,"HTO-1",group.by = "orig.ident",pt.size = 0,slot = "counts")+
  geom_hline(yintercept=30 )
sum(cbl12@assays$HTO@counts[1,]>30) #824/6500 (~13%) cells with High HTO count, expected 14%
cbl12$hto<-cbl12@assays$HTO@counts[1,]>30

VlnPlot(cbl12,c("nCount_RNA","nCount_ATAC"),group.by = "hto",pt.size = 0)

DimPlot(cbl12,reduction = "umap",group.by = "lineage_hmap")
DimPlot(cbl12,reduction = "ref.umap",group.by = "lineage_hmap")

saveRDS(cbl12,"outputs/21-CREs_inference_with_sc_multiomics_cbls/cbl12.rds")

#find LT-hsc cells ?
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
DefaultAssay(hmap)<-"integrated"

hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

DefaultAssay(cbl12)<-"SCT"

anchors <- FindTransferAnchors(
    reference = hmap,
    query = cbl12,
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )

cbl12 <- MapQuery(
    anchorset = anchors, 
    query = cbl12,
    reference = hmap,
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )


DimPlot(cbl12, reduction = "ref.umap", group.by =  "predicted.cell_type", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

#nop
