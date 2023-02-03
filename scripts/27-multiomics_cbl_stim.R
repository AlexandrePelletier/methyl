out<-"outputs/27-multiomics_cbl_stim"
dir.create(out)
source("scripts/utils/new_utils.R")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

##cblA####

# load the RNA and ATAC data
counts <- Read10X_h5("~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-a/outs/filtered_feature_bc_matrix.h5")
fragpath <- "~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-a/outs/atac_fragments.tsv.gz"
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
#genome(annotation) <- "hg38"

head(annotation)

# create a Seurat object containing the RNA adata
cblA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
cblA[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
cblA
# An object of class Seurat 
# 149849 features across 2944 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
#  1 other assay present: ATAC

#QC
DefaultAssay(cblA) <- "ATAC"

cblA <- NucleosomeSignal(cblA)
cblA <- TSSEnrichment(cblA)

VlnPlot(
  object = cblA,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
cblA <- subset(
  x = cblA,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1
)

cblA # 2001  cells

#peak calling 
# call peaks using MACS2
#reticulate::py_install(packages ="MACS2")

peaks <- CallPeaks(cblA, macs2.path = "renv/python/virtualenvs/renv-python-3.7/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(cblA),
  features = peaks,
  cells = colnames(cblA)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
cblA[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(cblA,fp(out,"cblA.rds"))

#Gene expression data processing
DefaultAssay(cblA) <- "RNA"
cblA <- SCTransform(cblA)
cblA <- RunPCA(cblA)
cblA <- CellCycleScoring(cblA,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)
cblA$CC.Difference <- cblA$S.Score - cblA$G2M.Score
cblA <- RunUMAP(cblA,dims = 1:30)
DimPlot(cblA,group.by="Phase")
cblA[["percent.mt"]] <- PercentageFeatureSet(object = cblA, assay = "RNA",pattern = "^MT-")
FeaturePlot(cblA,"percent.mt")
cblA <- SCTransform(cblA,vars.to.regress = c("percent.mt"))
cblA <- RunPCA(cblA)
cblA <- RunUMAP(cblA,dims = 1:30)

FeaturePlot(cblA,"percent.mt")

#DNA accessibility data processing
DefaultAssay(cblA) <- "peaks"
cblA <- FindTopFeatures(cblA, min.cutoff = 5)
cblA <- RunTFIDF(cblA)
cblA <- RunSVD(cblA)

# load CBPs reference (hematomap)
hmap <- readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
DefaultAssay(hmap) <- "integrated"

DefaultAssay(cblA) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = hmap,
  query = cblA,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:30
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = hmap$cell_type,
  weight.reduction = cblA[['pca']],
  dims = 1:30
)

cblA <- AddMetaData(
  object = cblA,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(cblA) <- "predicted.id"

# # set a reasonable order for cell types to be displayed when plotting
# levels(cblA) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
#                   "CD8 Naive", "dnT",
#                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
#                  "NK Proliferating", "gdT",
#                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
#                  "CD14 Mono", "CD16 Mono",
#                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

#Joint UMAP visualization
#Using the weighted nearest neighbor methods in Seurat v4, we can compute a joint neighbor graph that represent both the gene expression and DNA accessibility measurements.

# build a joint neighbor graph using both assays
cblA <- FindMultiModalNeighbors(
  object = cblA,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
cblA <- RunUMAP(
  object = cblA,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(cblA, label = TRUE, reduction = "umap") + NoLegend() #integration pas terrible

cblA$cell_type_hmap<-Idents(cblA)
cblA[["lineage_hmap"]]<-sapply(as.character(cblA@meta.data$cell_type_hmap), function(ct){
  if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
  else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
  else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
  else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
  else if(ct%in%c("CLP","proB"))return("Lymphoid")
  else return(ct)
  
})
DimPlot(cblA,group.by = "lineage_hmap", label = TRUE, reduction = "umap") + NoLegend()


FeaturePlot(cblA,c("AVP","MPO","CD99","GATA1"),reduction = "umap",max.cutoff = "q95")

FeaturePlot(cblA,c("AVP","SOCS3","JUNB","EGR1"),reduction = "umap",max.cutoff = "q95")

cblA <- RunUMAP(cblA,dims = 1:30,reduction.name = "umap_sct")
DimPlot(cblA, label = TRUE, reduction = "umap_sct") +
  NoLegend() # better with the joint umap 

DimPlot(cblA, label = TRUE, reduction = "umap") + NoLegend() #what are the cluster of cells unannoted in the left ?
cblA<-FindClusters(cblA,resolution = 0.6,graph.name = "wsnn")
DimPlot(cblA, label = TRUE, reduction = "umap") + NoLegend() #what are the cluster 2?

markers<-FindMarkers(cblA,ident.1 = 2)
head(markers,20) #RP related, so protein translation related
#because cells have been put in culture a night.
#need to rm this effect 
DimHeatmap(cblA,cells = 500,dims = 1:10)
#rm PC1
cblA <- FindMultiModalNeighbors(
  object = cblA,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(2:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
# build a joint UMAP visualization
cblA <- RunUMAP(
  object = cblA,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(cblA, label = TRUE, reduction = "umap",group.by = "cell_type_hmap") +
  NoLegend() #greater resemblance but still big RP effect

# #rm rp influence from sct
# rb.genes <- rownames(cblA)[grep("^RP[SL]",rownames(cblA))]
# cblA<-PercentageFeatureSet(cblA,pattern ="^RP[SL]",col.name = "percent.rb" )
# FeaturePlot(cblA,"percent.rb",max.cutoff = "q95")
# cblA<-
# percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
# alldata <- AddMetaData(alldata, percent.ribo, col.name = "percent.ribo")
# 
# DimPlot(cblA, label = TRUE, reduction = "umap",group.by = "lineage_hmap") 
# cblA <- SCTransform(cblA,vars.to.regress = c("percent.mt","percent.rb"))
# cblA <- RunPCA(cblA)
# cblA <- RunUMAP(cblA,dims = 1:30)
# 
# FeaturePlot(cblA,"percent.rb")
# 
# # transfer cell type labels from reference to query
# transfer_anchors <- FindTransferAnchors(
#   reference = hmap,
#   query = cblA,
#   normalization.method = "SCT",
#   reference.reduction = "pca",
#   recompute.residuals = FALSE,
#   dims = 1:30
# )
# 
# predictions <- TransferData(
#   anchorset = transfer_anchors, 
#   refdata = hmap$cell_type,
#   weight.reduction = cblA[['pca']],
#   dims = 1:30
# )
# 
# cblA <- AddMetaData(
#   object = cblA,
#   metadata = predictions
# )
# 
# # set the cell identities to the cell type predictions
# Idents(cblA) <- "predicted.id"
# 
# # # set a reasonable order for cell types to be displayed when plotting
# # levels(cblA) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
# #                   "CD8 Naive", "dnT",
# #                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
# #                  "NK Proliferating", "gdT",
# #                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
# #                  "CD14 Mono", "CD16 Mono",
# #                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")
# 
# #Joint UMAP visualization
# #Using the weighted nearest neighbor methods in Seurat v4, we can compute a joint neighbor graph that represent both the gene expression and DNA accessibility measurements.
# 
# # build a joint neighbor graph using both assays
# cblA <- FindMultiModalNeighbors(
#   object = cblA,
#   reduction.list = list("pca", "lsi"), 
#   dims.list = list(1:40, 2:30),
#   modality.weight.name = "RNA.weight",
#   verbose = TRUE
# )
# 
# # build a joint UMAP visualization
# cblA <- RunUMAP(
#   object = cblA,
#   nn.name = "weighted.nn",
#   assay = "RNA",
#   verbose = TRUE
# )
# 
# DimPlot(cblA, label = TRUE, reduction = "umap") + NoLegend() #integration pas terrible
# 
# cblA$cell_type_hmap<-Idents(cblA)
# cblA[["lineage_hmap"]]<-sapply(as.character(cblA@meta.data$cell_type_hmap), function(ct){
#   if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
#   else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
#   else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
#   else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
#   else if(ct%in%c("CLP","proB"))return("Lymphoid")
#   else return(ct)
#   
# })
# DimPlot(cblA,group.by = "lineage_hmap", label = TRUE, reduction = "umap") +
#   NoLegend() #loss biological variation (Myeloid mix with MPP now)
# 
# 
# FeaturePlot(cblA,c("AVP","MPO","CD99","GATA1"),reduction = "umap",max.cutoff = "q95",order = T)
# 
# FeaturePlot(cblA,c("AVP","SOCS3","JUNB","EGR1"),reduction = "umap",max.cutoff = "q95")

#ccl SCT percent mt better


saveRDS(cblA,fp(out,"cblA.rds"))



#make same for cblB and C and integrate with cbl12 (harmony), then map SCT or integrated matrix
#cblB####
# load the RNA and ATAC data
counts <- Read10X_h5("~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-b/outs/filtered_feature_bc_matrix.h5")
fragpath <- "~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-b/outs/atac_fragments.tsv.gz"


# create a Seurat object containing the RNA adata
cblB <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
cblB[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
cblB
# An object of class Seurat 
# 153901 features across 4003 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
#  1 other assay present: ATAC
#QC
DefaultAssay(cblB) <- "ATAC"

cblB <- NucleosomeSignal(cblB)
cblB <- TSSEnrichment(cblB)

VlnPlot(
  object = cblB,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
cblB <- subset(
  x = cblB,
  subset = nCount_ATAC < 50000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1
)

cblB # 2353  cells

VlnPlot(
  object = cblB,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

peaks <- CallPeaks(cblB, macs2.path = "renv/python/virtualenvs/renv-python-3.7/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(cblB),
  features = peaks,
  cells = colnames(cblB)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
cblB[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(cblB,fp(out,"cblB.rds"))

#Gene expression data processing
DefaultAssay(cblB) <- "RNA"
cblB <- SCTransform(cblB)
cblB <- RunPCA(cblB)
cblB <- CellCycleScoring(cblB,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)
cblB$CC.Difference <- cblB$S.Score - cblB$G2M.Score
cblB <- RunUMAP(cblB,dims = 1:30)
DimPlot(cblB,group.by="Phase")

cblB[["percent.mt"]] <- PercentageFeatureSet(object = cblB, assay = "RNA",pattern = "^MT-")
FeaturePlot(cblB,"percent.mt")
cblB<-PercentageFeatureSet(cblB,pattern ="^RP[SL]",col.name = "percent.rb" )
FeaturePlot(cblB,"percent.rb")

cblB <- SCTransform(cblB,vars.to.regress = c("percent.mt"))
cblB <- RunPCA(cblB)
cblB <- RunUMAP(cblB,dims = 1:30)

FeaturePlot(cblB,"percent.mt")

#DNA accessibility data processing
DefaultAssay(cblB) <- "peaks"
cblB <- FindTopFeatures(cblB, min.cutoff = 5)
cblB <- RunTFIDF(cblB)
cblB <- RunSVD(cblB)

# load CBPs reference (hematomap)
hmap <- readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
DefaultAssay(hmap) <- "integrated"

DefaultAssay(cblB) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = hmap,
  query = cblB,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:30
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = hmap$cell_type,
  weight.reduction = cblB[['pca']],
  dims = 1:30
)

cblB <- AddMetaData(
  object = cblB,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(cblB) <- "predicted.id"

cblB$cell_type_hmap<-Idents(cblB)
cblB[["lineage_hmap"]]<-sapply(as.character(cblB@meta.data$cell_type_hmap), function(ct){
  if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
  else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
  else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
  else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
  else if(ct%in%c("CLP","proB"))return("Lymphoid")
  else return(ct)
  
})
DimHeatmap(cblB,cells = 500,dims = 1:10)
#rm PC1
cblB <- FindMultiModalNeighbors(
  object = cblB,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(2:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
# build a joint UMAP visualization
cblB <- RunUMAP(
  object = cblB,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(cblB, label = TRUE, reduction = "umap",group.by = "cell_type_hmap") +
  NoLegend() #greater resemblance but still big RP effect
saveRDS(cblB,fp(out,"cblB.rds"))

#cblC####
# load the RNA and ATAC data
counts <- Read10X_h5("~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-c/outs/filtered_feature_bc_matrix.h5")
fragpath <- "~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-c/outs/atac_fragments.tsv.gz"



# create a Seurat object containing the RNA adata
cblC <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
cblC[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
cblC
# An object of class Seurat 
# 129370 features across 6431 samples within 2 assays 
# Active assay: RNA (36601 features, 0 variable features)
#  1 other assay present: ATAC
#QC
DefaultAssay(cblC) <- "ATAC"

cblC <- NucleosomeSignal(cblC)
cblC <- TSSEnrichment(cblC)

VlnPlot(
  object = cblC,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
VlnPlot(
  object = cblC,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC"),
  ncol = 2,
  pt.size = 0,
  log = T
)
# filter out low quality cells
cblC <- subset(
  x = cblC,
  subset = nCount_ATAC < 50000 &
    nCount_RNA < 5000 &
    nCount_ATAC > 100 &
    nCount_RNA > 100 &
    nucleosome_signal < 1 &
    TSS.enrichment > 1
)

cblC # 2599  cells

VlnPlot(
  object = cblC,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

peaks <- CallPeaks(cblC, macs2.path = "renv/python/virtualenvs/renv-python-3.7/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(cblC),
  features = peaks,
  cells = colnames(cblC)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
cblC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(cblC,fp(out,"cblC.rds"))

#Gene expression data processing
DefaultAssay(cblC) <- "RNA"
cblC <- SCTransform(cblC)
cblC <- RunPCA(cblC)
cblC <- CellCycleScoring(cblC,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)
cblC$CC.Difference <- cblC$S.Score - cblC$G2M.Score
cblC <- RunUMAP(cblC,dims = 1:30)
DimPlot(cblC,group.by="Phase")

cblC[["percent.mt"]] <- PercentageFeatureSet(object = cblC, assay = "RNA",pattern = "^MT-")
FeaturePlot(cblC,"percent.mt")
cblC<-PercentageFeatureSet(cblC,pattern ="^RP[SL]",col.name = "percent.rb" )
FeaturePlot(cblC,"percent.rb")

cblC <- SCTransform(cblC,vars.to.regress = c("percent.mt"))
cblC <- RunPCA(cblC)
cblC <- RunUMAP(cblC,dims = 1:30)

FeaturePlot(cblC,"percent.mt")

#DNA accessibility data processing
DefaultAssay(cblC) <- "peaks"
cblC <- FindTopFeatures(cblC, min.cutoff = 5)
cblC <- RunTFIDF(cblC)
cblC <- RunSVD(cblC)

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = hmap,
  query = cblC,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:30
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = hmap$cell_type,
  weight.reduction = cblC[['pca']],
  dims = 1:30
)

cblC <- AddMetaData(
  object = cblC,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(cblC) <- "predicted.id"

cblC$cell_type_hmap<-Idents(cblC)
cblC[["lineage_hmap"]]<-sapply(as.character(cblC@meta.data$cell_type_hmap), function(ct){
  if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
  else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
  else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
  else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
  else if(ct%in%c("CLP","proB"))return("Lymphoid")
  else return(ct)
  
})
DimHeatmap(cblC,cells = 500,dims = 1:10)
#rm PC1
cblC <- FindMultiModalNeighbors(
  object = cblC,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(2:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
# build a joint UMAP visualization
cblC <- RunUMAP(
  object = cblC,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(cblC, label = TRUE, reduction = "umap",group.by = "cell_type_hmap") +
  NoLegend() #greater resemblance but still big RP effect
saveRDS(cblC,fp(out,"cblC.rds"))

#INTEGRATION ####
cblA$condition<-"no_stim"

cblB$condition<-"stim_2h"

cblC$condition<-"stim_24h"

cbls_list<-list(cblA,cblB,cblC)

#map on hmap
library(parallel)
DefaultAssay(hmap)<-"integrated"
hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

cbls_list<-mclapply(cbls_list,function(x){
  #message("calculate CC.Difference for",x@project.name)
  if(!"S.Score"%in%colnames(x@meta.data)){
    x<-SCTransform(x,method = "glmGamPoi")
    x <- CellCycleScoring(x,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)


  }
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
  },mc.cores = 6)

cbls_list<-mclapply(cbls_list, SCTransform,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi",mc.cores = 6)


anchors <- list()
for (i in 1:length(cbls_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = hmap,
    query = cbls_list[[i]],
    k.filter = 200,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )
}


for (i in 1:length(cbls_list)) {
  cbls_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = cbls_list[[i]],
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )
}

#merge

cbls_stim <- merge(cbls_list[[1]], cbls_list[2:length(cbls_list)],merge.dr = c("ref.pca","ref.umap"))
DimPlot(cbls_stim, reduction = "ref.umap", group.by =  "predicted.cell_type", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

DimPlot(cbls_stim, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

DimPlot(cbls_stim, reduction = "ref.umap", group.by =  "predicted.id", label = TRUE, repel = TRUE, label.size = 3) 
#don't work well

#just transfer label

for (i in 1:length(cbls_list)) {
  cbls_list[[i]] <- SCTransform(cbls_list[[i]],vars.to.regress = c("percent.mt"))
  cbls_list[[i]] <- RunPCA(cbls_list[[i]])
}

for (i in 1:length(cbls_list)) {
  # transfer cell type labels from reference to query
  transfer_anchors <- FindTransferAnchors(
    reference = hmap,
    query = cbls_list[[i]],
    normalization.method = "SCT",
    reference.reduction = "pca",
    recompute.residuals = FALSE,
    dims = 1:30
  )
  
  predictions <- TransferData(
    anchorset = transfer_anchors, 
    refdata = hmap$cell_type,
    weight.reduction = cbls_list[[i]][['pca']],
    dims = 1:30
  )
  
  cbls_list[[i]] <- AddMetaData(
    object = cbls_list[[i]],
    metadata = predictions
  )
  
  # set the cell identities to the cell type predictions
  Idents(cbls_list[[i]]) <- "predicted.id"
  
    cbls_list[[i]]$cell_type_hmap<-Idents(cbls_list[[i]])
    cbls_list[[i]][["lineage_hmap"]]<-sapply(as.character(cbls_list[[i]]@meta.data$cell_type_hmap), function(ct){
      if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
      else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
      else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
      else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
      else if(ct%in%c("CLP","proB"))return("Lymphoid")
      else return(ct)
  })

}

#without stim 24h because too much different
cbls_stim <- merge(cbls_list[[1]], cbls_list[[2]])
saveRDS(cbls_stim,fp(out,"cbls_stim2.rds"))

# unsupervised clustering
features<-SelectIntegrationFeatures(cbls_list, nfeatures = 3000)

cbls_stim <- RunPCA(cbls_stim, assay = "SCT",features=features)
cbls_stim <- RunUMAP(cbls_stim, reduction = 'pca', dims = 1:30,reduction.name = "umap")
DimPlot(cbls_stim, group.by = c('condition',"lineage_hmap"),reduction = "umap", label = TRUE)

#harmony ?
library(harmony)
cbls_stim<-RunHarmony(cbls_stim,reduction = "pca",group.by.vars = "condition",assay.use = "SCT")
cbls_stim <- RunUMAP(cbls_stim, reduction = 'harmony', dims = 1:30,reduction.name = "humap")
DimPlot(cbls_stim, group.by = c('condition',"lineage_hmap"),reduction = "humap", label = TRUE)

#DimPlot(cbls_stim, group.by = 'cell_type_hmap',reduction = "humap", label = TRUE)
#DimPlot(cbls_stim, group.by = 'lineage_hmap',reduction = "humap", label = TRUE)
table(cbls_stim$condition,cbls_stim$lineage_hmap)
  #            18 B cell   DC Erythro-Mas  HSC LT-HSC Lymphoid Mk/Er MPP/LMPP Myeloid T cell
  # no_stim     3     22   10         218  435     46       58     1     1014     160     34
  # stim_24h   32     22   15         522  537     12      109     5     1039     266     40
  # stim_2h     3     36   23         220  532      0       62     1     1230     212     34
VlnPlot(cbls_stim, features=c("EGR1","MPO","LTB","GATA2"), group.by = 'lineage_hmap')

#harmony lsi
DefaultAssay(cbls_stim) <- "peaks"
cbls_stim <- FindTopFeatures(cbls_stim, min.cutoff = 5)
cbls_stim <- RunTFIDF(cbls_stim)
cbls_stim <- RunSVD(cbls_stim)
cbls_stim<-RunHarmony(cbls_stim,reduction = "lsi",group.by.vars = "condition",assay.use = "peaks",
                      reduction.save = "harmony_lsi", project.dim = FALSE)

#DimHeatmap(cbls_stim,reduction = "harmony",cells=500,dims = 1:6)
cbls_stim <- FindMultiModalNeighbors(
  object = cbls_stim,
  reduction.list = list("harmony", "harmony_lsi"), 
  dims.list = list(1:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)


# build a joint UMAP visualization
cbls_stim <- RunUMAP(
  object = cbls_stim,
  reduction.name = "multi_humap",
  nn.name = "weighted.nn",
  assay = "SCT",
  verbose = TRUE
)
cbls_stim<-FindClusters(cbls_stim,resolution = 0.8,graph.name = "wsnn")
cbls_stim$multi_clusters<-Idents(cbls_stim)
DimPlot(cbls_stim, group.by = c('lineage_hmap','multi_clusters'),reduction = "multi_humap", label = TRUE)

ggplot(cbls_stim@meta.data)+geom_bar(aes(x=multi_clusters,fill=lineage_hmap),position = "fill")
#ggplot(cbls_stim@meta.data)+geom_bar(aes(x=multi_clusters,fill=cell_type_hmap),position = "fill")

DefaultAssay(cbls_stim)<-"SCT"
FeaturePlot(cbls_stim,features = c("EGR1",'SOCS3'),
            reduction = "multi_humap",label=T)
FeaturePlot(cbls_stim,features = c("AVP",'DUSP2'),
            reduction = "multi_humap",label=T)

FeaturePlot(cbls_stim,features = c('MPO',"CEBPA"),reduction = "multi_humap")
FeaturePlot(cbls_stim,features = c('LTB',"CD99"),reduction = "multi_humap",label=T)

FeaturePlot(cbls_stim,features = c('GATA2',"GATA1"),reduction = "multi_humap",label=T)

new.idents<-c(
  "0" = "HSC",
  "1"="MPP/LMPP",
  "2" = "2",
  "3" = "Myeloid",
  "4" = "Lymphoid",
  "5"="5",
  "6"="Erythro-Mas",
  "7"= "7",
  "8"='8',
  "9"="9",
  "10"="Lymphoid",
  "11"="T cell",
  "12"="B cell",
  "13"="13",
  "14"="DC",
  "15"="15",
  '16'='16'
)

cbls_stim<-RenameIdents(cbls_stim,new.idents)
cbls_stim$ident.1<-Idents(cbls_stim)
# reticulate::py_install("leidenalg")
# cbls_stim<-FindClusters(cbls_stim,resolution = 0.8,graph.name = "wsnn",algorithm = 4)
# DimPlot(cbls_stim,reduction = "multi_humap", label = TRUE)

DimPlot(cbls_stim,reduction = "multi_humap", label = TRUE)
FeaturePlot(cbls_stim,features = c('percent.rb'),reduction = "multi_humap")

#annot undefined clusters
markers<-data.table(FindAllMarkers(cbls_stim,assay = "SCT"))
markers[cluster=="7"][1:20]
markers[cluster=="5"][1:20]
markers[cluster=="8"][1:20]
markers[cluster=="9"][1:20]
markers[cluster=="2"][1:20] #DC progen
markers[cluster=="15"][1:20] #MT cells
markers[cluster=="17"][1:20] #MT cells 2


markers_ct<-fread("../singlecell/ref/hematopo_markers.csv")
markersa<-merge(markers,markers_ct,all.x=T)
markersa[cluster=="13"][order(p_val_adj)][1:20] #CD4+ T cell

markersa[cluster=="16"][order(p_val_adj)][1:20] #
markersa[cluster=="16"][order(p_val_adj)][!is.na(cell_type)] #Mac?

markersa[cluster=="9"][order(p_val_adj)][1:20] #CD4+ T cell 2
markersa[cluster=="9"][order(p_val_adj)][!is.na(cell_type)] #

markersa[cluster=="2"][order(p_val_adj)][1:20] #
markersa[cluster=="2"][order(p_val_adj)][!is.na(cell_type)] # DC Prog
FeaturePlot(cbls_stim,"GATA1",reduction = "multi_humap")

markersa[cluster=="8"][order(p_val_adj)][1:20] #
markersa[cluster=="8"][order(p_val_adj)][!is.na(cell_type)] # Mk P


markersa[cluster=="7"][order(p_val_adj)][1:20] #
markersa[cluster=="7"][order(p_val_adj)][!is.na(cell_type)] # MPP-cycle
FeaturePlot(cbls_stim,c("PROM1","ZBTB16","KIT","PLCG2"),reduction = "multi_humap")

markersa[cluster=="5"][order(p_val_adj)][1:40] #HSPC-RP
markersa[cluster=="5"][order(p_val_adj)][!is.na(cell_type)] #
FeaturePlot(cbls_stim,c("KIT","HBA2"),reduction = "multi_humap")

DimPlot(cbls_stim,group.by = "Phase",reduction = "multi_humap")
head(cbls_stim[[]])

new.idents<-c(
  "2" = "DC Prog",
  "5"="HSPC-RP",
  "7"= "MPP-cycle",
  "8"='HSC-Mk',
  "9"="CD4+ T cell",
  "13"="CD4+ T cell 2",
  "15"="MT cells",
  '16'='Mac?',
  "17"="MT cells 2"
)
cbls_stim<-RenameIdents(cbls_stim,new.idents)
cbls_stim$ident.2<-Idents(cbls_stim)
DimPlot(cbls_stim,reduction = "multi_humap", label = TRUE)
FeaturePlot(cbls_stim,c("LTB","VPREB1"),
            reduction = "multi_humap",max.cutoff = "q90")

saveRDS(cbls_stim,fp(out,"cbls_stim.rds"))



#EFFECT of STIMULATION ####
#2h stim vs 2h nostim
cbls<-readRDS("outputs/27-multiomics_cbl_stim/cbls_stim.rds")
DimPlot(cbls,reduction = "multi_humap",label=T)
DimPlot(cbls,reduction = "multi_humap",group.by="Phase",label=T)

#diff of subpop ?####
ggplot(cbls@meta.data)+geom_bar(aes(x=condition,
                                    fill=ident.2),position = "fill")

mtds<-data.table(cbls@meta.data,keep.rownames = "bc")
mtds[,n.sample:=.N,by="condition"]

mtds[,pct.ct:=.N/n.sample,by=c("condition","ident.2")]

levels(cbls)
FeaturePlot(cbls,c("ID1","AVP","SOCS3","EGR1"),min.cutoff = "q10",max.cutoff = "q95",reduction = "multi_humap",order=T,label=T)
FeaturePlot(cbls,c("CDK6"),min.cutoff = "q10",max.cutoff = "q95",reduction = "multi_humap",order=T,label=T)
FeaturePlot(cbls,c("MLLT3"),min.cutoff = "q10",max.cutoff = "q95",reduction = "multi_humap",order=T,label=T)
FeaturePlot(cbls,c("EGR1"),min.cutoff = "q10",max.cutoff = "q95",reduction = "multi_humap",order=T,label=T)
FeaturePlot(cbls,c("percent.mt"),reduction = "multi_humap",label=T)

FeaturePlot(cbls,c("CD8"),min.cutoff = "q10",max.cutoff = "q95",reduction = "multi_humap",order=T,label=T)
res<-FindMarkers(cbls,ident.1 = "T cell",ident.2 = c("CD4+ T cell","CD4+ T cell 2"))
res[res$p_val_adj<0.05,]
lvls<-c("HSC","HSC-Mk","MPP/LMPP","Lymphoid","Myeloid","Erythro-Mas","DC Prog")
mtdsf<-mtds[ident.2%in%lvls]

mtdsf[,lineage:=factor(ident.2,levels = lvls)]
ggplot(unique(mtdsf,by=c("condition","lineage")))+
  geom_col(aes(x=lineage,y=pct.ct,fill=condition),position='dodge')


#global transcription####
unique(cbls$condition)
res<-data.table(FindMarkers(cbls,group.by="condition",ident.1="stim_2h"),keep.rownames = "gene")
res[avg_log2FC>0][1:40]
table(res[p_val_adj<0.001]$avg_log2FC>0)
# FALSE  TRUE 
#   514   281 
res[p_val_adj<0.001&avg_log2FC>0]$gene

#lineage spe transcription####
res<-Reduce(rbind,lapply(lvls,function(ct)data.table(FindMarkers(cbls,
                                group.by="condition",
                                ident.1="stim_2h",
                                subset.ident = ct),keep.rownames = "gene")[,lineage:=ct]))

table(res[p_val_adj<0.05]$lineage)
  # DC Prog Erythro-Mas         HSC      HSC-Mk    Lymphoid    MPP/LMPP     Myeloid 
  #       752         495         581          52         414         670         396 
res[,lineage:=factor(lineage,levels = lvls)]
ggplot(res[p_val_adj<0.05])+geom_bar(aes(x=lineage,fill=avg_log2FC>0),position = "dodge")
res_hsc<-res[lineage=="HSC"]
res_hsc[avg_log2FC>0&p_val_adj<0.05]
# BHLHE40 4.450740e-68  1.6303004 0.808 0.195 1.004443e-63
#  2:     IL2RA 7.826327e-61  2.1002804 0.683 0.102 1.766246e-56
#  3:    IL1RAP 3.002129e-59  1.5955557 0.727 0.132 6.775204e-55
#  4:      CD69 2.903954e-58  2.1083886 0.915 0.491 6.553644e-54
#  5:      BCL2 7.357406e-58  2.0549240 0.738 0.165 1.660419e-53
#  6:    SLC2A3 1.082858e-56  1.1810779 0.568 0.030 2.443795e-52
#  7:      LCP2 1.096051e-50  1.7843891 0.915 0.596 2.473568e-46
#  8:       AHR 4.396296e-50  1.5001468 0.952 0.613 9.921562e-46
#  9:     PRR5L 2.990148e-48  1.2530315 0.598 0.095 6.748166e-44
# 10:     SOCS2 1.002651e-47  1.3655687 0.690 0.190 2.262783e-43
# 11:    ZNF608 1.964223e-45  1.4810815 0.690 0.204 4.432859e-41
# 12:     NFAT5 1.032965e-44  1.1123009 1.000 1.000 2.331196e-40
# 13:  MIR155HG 1.298810e-44  1.5508363 0.701 0.207 2.931154e-40
# 14:    MAP4K3 2.131921e-43  1.3424428 0.963 0.840 4.811318e-39
# 15:     HDAC4 3.525387e-43  1.1821080 0.657 0.180 7.956094e-39
# 16:      PRKX 1.676129e-39  1.1900518 0.565 0.115 3.782688e-35
# 17:     PDE4B 1.951790e-39  1.6047605 0.590 0.137 4.404800e-35
# 18:     CYTIP 5.843758e-39  1.4819810 0.613 0.157 1.318819e-34
# 19:    CASP10 6.408574e-37  0.9809874 0.716 0.282 1.446287e-32
# 20:      STK3 4.192562e-36  1.1153834 0.948 0.723 9.461773e-32
# 21:     DUSP6 5.131293e-36  0.8121466 0.421 0.037 1.158030e-31
# 22:       CFH 5.789217e-36  1.2587129 0.683 0.244 1.306511e-31
# 23:     PTGS2 7.611609e-36  1.4075827 0.502 0.092 1.717788e-31
# 24:      EGR1 6.099288e-34  1.0050089 0.402 0.032 1.376487e-29
# 25:     CLIP1 1.133406e-33  1.0075317 0.834 0.451 2.557871e-29
# 26: MAP4K3-DT 2.982554e-33  1.1064176 0.568 0.150 6.731028e-29
# 27:   DENND4A 4.308448e-33  1.1591597 0.982 0.883 9.723306e-29
# 28:     PTBP3 4.071283e-32  0.9614828 0.908 0.623 9.188071e-28
# 29:      EBF4 4.766243e-32  1.1038126 0.661 0.282 1.075646e-27
# 30:     STAT3 1.612967e-30  0.9104706 0.937 0.771 3.640144e-26
# 31:    FNDC3A 3.818236e-30  1.0061841 0.982 0.898 8.616995e-26
# 32:     POC1B 2.662243e-29  0.9661114 0.749 0.349 6.008151e-25
# 33:     SBNO2 1.336698e-28  0.7556747 0.539 0.160 3.016660e-24
# 34:      TCHH 1.884640e-28  0.5831869 0.273 0.000 4.253254e-24
# 35:     NABP1 7.978733e-28  0.8367724 0.579 0.195 1.800640e-23
# 36:     SOCS3 7.035601e-27  0.5688760 0.273 0.005 1.587794e-22
# 37:     TAFA2 1.146266e-26  1.1265603 0.605 0.219 2.586893e-22
# 38:      IER3 4.926587e-26  0.7855755 0.506 0.147 1.111832e-21
# 39:   LAPTM4B 8.119249e-26  0.8907663 0.797 0.444 1.832352e-21
# 40:     RRAS2 1.962653e-25  0.8495764 0.720 0.369 4.429315e-21
#          gene        p_val avg_log2FC pct.1 pct.2    p_val_adj
table(res_hsc[p_val_adj<0.05]$avg_log2FC>0)
# FALSE  TRUE 
#   202   187 

res_hsc[p_val_adj<0.001&avg_log2FC>0]$gene
# [1] "BHLHE40"      "IL2RA"        "IL1RAP"       "CD69"         "BCL2"        
# [6] "SLC2A3"       "LCP2"         "AHR"          "PRR5L"        "SOCS2"       
# [11] "ZNF608"       "NFAT5"        "MIR155HG"     "MAP4K3"       "HDAC4"       
# [16] "PRKX"         "PDE4B"        "CYTIP"        "CASP10"       "STK3"        
# [21] "DUSP6"        "CFH"          "PTGS2"        "EGR1"         "CLIP1"       
# [26] "MAP4K3-DT"    "DENND4A"      "PTBP3"        "EBF4"         "STAT3"       
# [31] "FNDC3A"       "POC1B"        "SBNO2"        "TCHH"         "NABP1"       
# [36] "SOCS3"        "TAFA2"        "IER3"         "LAPTM4B"      "RRAS2"       
# [41] "IL4R"         "CISH"         "BCL2L1"       "SOD2"         "HUWE1"       
# [46] "DENND3"       "FLOT1"        "YES1"         "KSR1"         "PLSCR1"      
# [51] "BATF"         "LTA"          "MCTP2"        "PRKCE"        "FOSL2"       
# [56] "BTBD7"        "CCND2"        "LATS2"        "PPP1R16B"     "ARNTL2-AS1"  
# [61] "FRMD6"        "KLF6"         "CGAS"         "FAM241A"      "RAPGEF6"     
# [66] "PPP6R3"       "HAPLN3"       "FOSL1"        "PIM1"         "MYO1E"       
# [71] "DDX21"        "XIST"         "HIVEP2"       "EMP1"         "PTPRC"       
# [76] "SOCS1"        "MAP3K20"      "STAT4"        "RBM27"        "IL6ST"       
# [81] "SOCS2-AS1"    "BHLHE40-AS1"  "AC090115.1"   "NAMPT"        "ADAMTSL4-AS1"
# [86] "VMP1"         "DEPP1"        "UBE2W"        "ERI2"         "USP32"       
# [91] "ARAP2"        "TNFRSF10D"    "CDKN1A"       "SPRED1"       "MYC"         
# [96] "SAMSN1"       "RHOH"         "ANKRD12"      "KCNT2"        "LYRM1"       
# [101] "RAB11FIP1"    "FAM83D"       "XBP1"         "PHTF2"        "NLRC5"       
# [106] "NIBAN2"       "CH25H"        "TRPM6"        "ETS2"         "AFDN"        
# [111] "EIF4A1"       "CYTOR"        "APP"          "IL3RA"        "MFHAS1"      
# [116] "PPFIBP1"      "SPRED2"       "BACE2"        "GALNT4"       "CCNL1"       
# [121] "PRNP"         "MIR17HG"      "CTBP2"        "ZSWIM6"       "LRRC8B"      
# [126] "PRDM8"        "NBEAL1"       "MAST4"        "SLC7A1"       "PABPC1"      
# [131] "SLCO4A1"      "FOS"          "CASP8"        "DMD"          "ICAM1"       
# [136] "LDLR"         "HS6ST2"       "ZMYND8"       "CDKL4"        "NR4A3"       
# [141] "FLT1"         "JUN"          "XRN2"         "ZFP36"        "YBX3"        
# [146] "RRBP1"        "NFKB1"        "SCARF1"       "SAMD4A"       "CHRM3"       
# [151] "SLC9A7"       "LTB"          "INTS14"       "CEMIP2"       "MYO9B"       
# [156] "CHSY1"        "POC1B-AS1"    "ZHX2"         "USP15"        "IRF1"        
# [161] "MB21D2"       "LITAF"        "ARIH1"        "MCL1"         "PPARD"       
# [166] "ETV6"         "SMG1"         "MGAT5"        "MMP25-AS1"    "MAP3K14"     
# [171] "SLC7A5"       "SLC1A5"       "CA13"         "RAB27A"       "TAB2"        
# [176] "SKI"          "MAFF"         "TAF4B"        "PITRM1"       "CAB39"       
# [181] "SLC20A1"      "NDFIP2"       "FNDC3B"       "TRAF6"        "GTF2F2"      
# [186] "STON2"        "PAPSS1"  

res_hsc[p_val_adj<0.001&avg_log2FC<0]$gene

fwrite(res,fp(out,'res_lineage_stim_vs_no_stim.csv.gz'))

#enrichment for which regulons ?
regulons<-fread("outputs/10-SCENIC/regulons.csv")

hsc_unstim<-colnames(cbls)[cbls$condition=="no_stim"&cbls$ident.2=="HSC"]
length(hsc_unstim) #401
hsc_stim<-colnames(cbls)[cbls$condition=="stim_2h"&cbls$ident.2=="HSC"]
length(hsc_stim) #271

genes_express_hsc<-rownames(cbls)[rowSums(cbls@assays$RNA[,hsc_unstim]>0)>0.1*length(hsc_unstim)|rowSums(cbls@assays$RNA[,hsc_stim]>0)>0.1*length(hsc_stim)]
length(genes_express_hsc) #7867
res<-OR3(querys = res_hsc[p_val_adj<0.05&avg_log2FC>0]$gene,
         terms_list = split(regulons$gene,regulons$tf),
         background = intersect(genes_express_hsc,unique(regulons$gene)))
res[padj<0.05][order(pval)]
#       term term.size n.query n.overlap pct.query.overlap  precision pct.term.overlap
#  1:    JUN        35      33        13        0.39393939 0.39393939        0.3714286
#  2:   FOSB        19      33        10        0.30303030 0.30303030        0.5263158
#  3: ARID5A        10      33         7        0.21212121 0.21212121        0.7000000
#  4:    FOS        32      33        10        0.30303030 0.30303030        0.3125000
#  5:   JUNB        40      33        11        0.33333333 0.33333333        0.2750000
#  6:   BCL3        26      33         9        0.27272727 0.27272727        0.3461538
#  7:  NFKB1         9      33         5        0.15151515 0.15151515        0.5555556
#  8:   JUND        56      33        11        0.33333333 0.33333333        0.1964286
#  9:   ATF3        26      33         7        0.21212121 0.21212121        0.2692308
# 10:  NFKB2        35      33         8        0.24242424 0.24242424        0.2285714
# 11:   XBP1       104      33        14        0.42424242 0.42424242        0.1346154
# 12: HOXA10        41      33         8        0.24242424 0.24242424        0.1951220
# 13:  CREB5        10      33         4        0.12121212 0.12121212        0.4000000
# 14:   EGR1        10      33         4        0.12121212 0.12121212        0.4000000
# 15:   IRF9        17      33         5        0.15151515 0.15151515        0.2941176
# 16:  STAT3        17      33         5        0.15151515 0.15151515        0.2941176
# 17:  HOXA9        71      33        10        0.30303030 0.30303030        0.1408451
# 18:  KLF10        21      33         5        0.15151515 0.15151515        0.2380952
# 19:   KLF2        24      33         5        0.15151515 0.15151515        0.2083333
# 20:  FOSL2         3      33         2        0.06060606 0.06060606        0.6666667
fwrite(res,fp(out,"res_regulons_enrichment_in_upregulated_genes_stim_hsc2h.csv.gz"))


ggplot(res[padj<0.05][order(padj)])+
  geom_point(aes(x=term,y=-log10(pval),size=n.overlap,col=fold.enrichment))+
  coord_cartesian(ylim = c(0,10))+scale_x_discrete(limits=c(res[padj<0.05][order(pval)]$term))

#GRN validation####
#genes of EGR1/KLF2/KLF4 network
egr1n<-unique(regulons[tf%in%c("EGR1","KLF2","KLF4")]$gene)

#findLinked peaks
DefaultAssay(cbls) <- "peaks"

# first compute the GC content for each peak
cbls <- RegionStats(cbls, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
cbls <- LinkPeaks(
  object = cbls,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = egr1n
)

res_peak_genes<-data.table(as.data.frame(Links(cbls)))
res_peak_genes #267
fwrite(res_peak_genes,fp(out,"res_peak_genes_egr1_klf2_klf4_regulons_linkage.csv"))

#same for all genes in regulons

cbls <- LinkPeaks(
  object = cbls,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = unique(regulons$gene)
)

res_peak_genes<-data.table(as.data.frame(Links(cbls)))
res_peak_genes #4281
res_peak_genes[pvalue<0.05]
fwrite(res_peak_genes,fp(out,"res_peak_genes_all_genes_regulons_linkage.csv"))
res_peak_genes<-fread(fp(out,"res_peak_genes_all_genes_regulons_linkage.csv"))

length(unique(res_peak_genes$gene)) #1487
length(unique(regulons$gene)) #2476
1487/2476
#n links in EGR1N
egr1n<-fread("outputs/16-GRN_final/egr1_KLF2_KLF4_network_tf_target_interactions.csv")
egr1n_genes<-unique(c(egr1n$target,egr1n$tf))
egr1n_genes# 123 genes

unique(res_peak_genes[gene%in%egr1n_genes],by="gene") #78/123 egr1n genes have linked peak (before was 41)

unique(res_peak_genes[gene%in%egr1n_genes]$gene)
#  [1] "ID3"      "RCAN3"    "ZC3H12A"  "PTCH2"    "JUN"      "RAP1A"    "RHOC"     "LMNA"     "IER5"     "ID2"      "ZFP36L2" 
# [12] "PELI1"    "ARID5A"   "DUSP2"    "CXCR4"    "NR4A2"    "TIPARP"   "SKIL"     "HES1"     "CCNG2"    "PRDM8"    "PTGER4"  
# [23] "IL6ST"    "EGR1"     "DUSP1"    "SQSTM1"   "TRIM41"   "IER3"     "HLA-E"    "GNL1"     "HSPA1A"   "CDKN1A"   "CITED2"  
# [34] "EZR"      "HSPB1"    "IFRD1"    "TACC1"    "ZBTB10"   "KLF10"    "NDRG1"    "SPTAN1"   "NFKB2"    "TRIM8"    "CDKN1C"  
# [45] "AHNAK"    "ETS1"     "EPS8"     "TUBA1A"   "GRASP"    "BTG1"     "UBL3"     "TSC22D1"  "ZFP36L1"  "FOS"      "KLF13"   
# [56] "AKAP13"   "SOCS1"    "PLCG2"    "EIF1"     "SOCS3"    "TGIF1"    "NFATC1"   "MIDN"     "GADD45B"  "JUNB"     "IER2"    
# [67] "DNAJB1"   "ARRDC2"   "ZFP36"    "SERTAD3"  "HNRNPUL1" "FOSB"     "MYADM"    "PRNP"     "RBM39"    "MAFF"     "TSC22D3" 
# [78] "IDS"     

#n of peaks by genes
npg<-table(res_peak_genes[gene%in%egr1n_genes]$gene)
npg[order(-as.vector(npg))]
# ZFP36L2   ARRDC2    CXCR4     HES1    SOCS1    SOCS3   TIPARP      ID3    PLCG2    AHNAK     ETS1    TRIM8  ZFP36L1  GADD45B 
#       26       14        9        8        8        8        8        7        7        6        6        6        6        5 
#     IER3     JUNB    NFKB2    PTCH2      FOS    GRASP      ID2     MAFF  TSC22D1    ZFP36   AKAP13   ARID5A   DNAJB1     EGR1 
#        5        5        5        5        4        4        4        4        4        4        3        3        3        3 
#     EIF1    HSPB1     IER2      JUN     LMNA   NFATC1    RCAN3   SPTAN1    CCNG2   CDKN1A   CDKN1C     EPS8     GNL1     IER5 
#        3        3        3        3        3        3        3        3        2        2        2        2        2        2 
#    KLF13    MYADM    NDRG1   PTGER4    RBM39   SQSTM1    TACC1   TRIM41  TSC22D3   ZBTB10     BTG1   CITED2    DUSP1    DUSP2 
#        2        2        2        2        2        2        2        2        2        2        1        1        1        1 
#      EZR     FOSB    HLA-E HNRNPUL1   HSPA1A      IDS    IFRD1    IL6ST    KLF10     MIDN    NR4A2    PELI1    PRDM8     PRNP 
#        1        1        1        1        1        1        1        1        1        1        1        1        1        1 
#    RAP1A     RHOC  SERTAD3     SKIL    TGIF1   TUBA1A     UBL3  ZC3H12A 
#        1        1        1        1        1        1        1        1 



 #CCL : 78/123 of Egr1/KLF2/KLF4 regulons have CREs

#edit 23/11/22 
#DMCs DEGs enrichment in peak down prev atac but with peak-gene link of multiome

#linkPeaks for every closest genes of background used for enrichment (peaks test for DA pct.det>0.05)
#run 27A


stims<-readRDS("outputs/27-multiomics_cbl_stim/cbls_stim_final.rds")


pgl<-fread("outputs/27-multiomics_cbl_stim/res_peak_genes_all_linkage.csv.gz")
pgl[,start.peak:=start(peak)]
pgl[,end.peak:=end(peak)]

#stats
pgl#15k peak gene link

length(unique(pgl$peak))#13k peaks
length(unique(pgl$gene))#7k genes


#overlap this CREs  with peak scATAC-seq

peaks<-fread("outputs/14-DMCs_atac_integr/cbps_atacs_peaks_by_lineage.csv.gz")
peaks_in_cres<-bed_inter(peaks[,.(seqnames,start,end,peaks)][order(seqnames,start)],
          unique(pgl[,.(seqnames,start.peak,end.peak,peak)][order(seqnames,start.peak,end.peak)]),
          select = c(4,1,8,6,7),col.names = c("peak_scatac","chr","peak_scmulti","start.peak","end.peak"))

peaks_in_cres #13k overlap
#N of peak overlap 
nrow(peaks_in_cres)/length(unique(pgl$peak))
unique(peaks_in_cres)
pgla<-merge(pgl[,peak_scmulti:=peak],peaks_in_cres)


#N peak-genes validated by multiome
#across the 13k CREs
peaks_hsc_anno<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")
phaf<-peaks_hsc_anno[query_region%in%pgla$peak_scatac]
phaf#11k
phaf[,well_linked_genes:=gene_name%in%pgla[peak_scatac==query_region]$gene,by="query_region"]
table(unique(phaf,by="query_region")$well_linked_genes)
# FALSE  TRUE 
#  8241  3198 
3198 / 8241 #38%
pgla

#DMCs/DEGs enrichemnt

#merge with res_degs
peaks_hsc_genes<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")

res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")
pgla_e<-merge(pgla,res_degs[lineage=="HSC"],by = "gene")
pgla_e#11143 peaks asso to express genes
length(unique(pgla_e$gene)) #4786 genes

#merge with res_cpgs
peaks_cpgs<-fread("outputs/14-DMCs_atac_integr/cpgs_in_lin_OCRs.csv.gz")
peaks_meth<-merge(peaks_cpgs,fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz"))
peaks_meth[,peak_scatac:=peaks]
peaks_methf<-peaks_meth[peak_scatac%in%pgla$peak_scatac]
length(unique(peaks_methf$peaks)) #7597 peaks
length(unique(peaks_methf$cpg_id)) #52339 CpGs


#term list : peaks DMCs/DEGs
meth_degs_cres_list<-list(CpGs_inCREs=unique(peaks_methf$peaks),
                         DMCs_inCREs=unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                         DEGs_withCREs=unique(pgla_e[padj<0.05&abs(log2FoldChange)>0.5]$peak_scatac),
                         DMCs_inCREs_with_DEGs=intersect(unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(pgla_e[padj<0.05&abs(log2FoldChange)>0.5]$peak_scatac)))
lapply(meth_degs_cres_list, length)
# $CpGs_inCREs
# [1] 7597
# 
# $DMCs_inCREs
# [1] 651
# 
# $DEGs_withCREs
# [1] 327
# 
# $DMCs_inCREs_with_DEGs
# [1] 18
# 

lapply(meth_degs_cres_list, head)

#query list : peaks up dn
peaks_hsc_lga_xy<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_hsc_without_xy.csv.gz")
peaks_hsc_lga_xyf<-peaks_hsc_lga_xy[peak%in%pgla$peak_scatac]
peaks_hsc_lga_xyf#106 /290

da_cres_list<-list(da_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&abs(avg_log2FC)>0.25]$peak),
                down_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC<(-0.25)]$peak),
                up_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC>0.25]$peak))
lapply(da_cres_list, length)
# $da_CREs
# [1] 102
# 
# $down_CREs
# [1] 85
# 
# $up_CREs
# [1] 17

lapply(da_cres_list, head)

#background : peaks in CREs

background<-unique(pgla_e$peak_scatac)
length(background) #10k peaks
head(background)

res_or_da_cres_dmcs_degs<-OR3(da_cres_list,
                    meth_degs_cres_list,
                    background = background,
                    overlap_column = F) 
res_or_da_cres_dmcs_degs

fwrite(res_or_da_cres_dmcs_degs,fp(out,"res_da_cres_enrichment_for_dmcs_degs.csv"))

#2
meth_degs_cres_list<-list(DMCs=unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                         DEGs=unique(pgla_e[padj<0.05&abs(log2FoldChange)>0.5]$peak_scatac),
                         upreg=unique(pgla_e[padj<0.05&log2FoldChange>0.5]$peak_scatac),
                         downreg=unique(pgla_e[padj<0.05&log2FoldChange<(-0.5)]$peak_scatac),
                         DMCs_DEGs=intersect(unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(pgla_e[padj<0.05&abs(log2FoldChange)>0.5]$peak_scatac)),
                         DMCs_upreg=intersect(unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(pgla_e[padj<0.05&log2FoldChange>0.5]$peak_scatac)),
                         DMCs_downreg=intersect(unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(pgla_e[padj<0.05&log2FoldChange<(-0.5)]$peak_scatac)))
lapply(meth_degs_cres_list, length)
# $CpGs
# [1] 7597
# 
# $DMCs
# [1] 651
# 
# $DEGs
# [1] 327
# 
# $upreg
# [1] 96
# 
# $downreg
# [1] 231
# 
# $DMCs_DEGs
# [1] 18
# 
# $DMCs_upreg
# [1] 4
# 
# $DMCs_downreg
# [1] 14


#query list : peaks up dn
peaks_hsc_lga_xy<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_hsc_without_xy.csv.gz")
peaks_hsc_lga_xyf<-peaks_hsc_lga_xy[peak%in%pgla$peak_scatac]
peaks_hsc_lga_xyf#106 /290

da_cres_list<-list(da_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&abs(avg_log2FC)>0.25]$peak),
                down_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC<(-0.25)]$peak),
                up_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC>0.25]$peak))
lapply(da_cres_list, length)
# $da_CREs
# [1] 102
# 
# $down_CREs
# [1] 85
# 
# $up_CREs
# [1] 17

lapply(da_cres_list, head)

#background : peaks in CREs

background<-unique(pgla_e$peak_scatac)
length(background) #10k peaks
head(background)

res_or_da_cres_dmcs_degs<-OR3(da_cres_list,
                    meth_degs_cres_list,
                    background = background,
                    overlap_column = F) 
res_or_da_cres_dmcs_degs

fwrite(res_or_da_cres_dmcs_degs,fp(out,"res_da_cres_enrichment_for_dmcs_degs2.csv"))


#        query.                  term term.size n.query n.overlap pct.query.overlap  precision pct.term.overlap background_size pct.term.background         pval         padj
#  1:   da_CREs           CpGs_inCREs      6034      97        90        0.92783505 0.92783505      0.014915479            9699         0.622125992 3.740578e-12 1.496231e-11
#  2:   da_CREs           DMCs_inCREs       540      97        16        0.16494845 0.16494845      0.029629630            9699         0.055675843 8.094559e-05 1.618912e-04
#  3:   da_CREs         DEGs_withCREs       327      97         7        0.07216495 0.07216495      0.021406728            9699         0.033714816 4.515888e-02 4.515888e-02
#  4:   da_CREs DMCs_inCREs_with_DEGs        18      97         2        0.02061856 0.02061856      0.111111111            9699         0.001855861 1.364740e-02 1.819654e-02
#  5: down_CREs           CpGs_inCREs      6034      82        82        1.00000000 1.00000000      0.013589659            9699         0.622125992 1.016335e-17 4.065340e-17
#  6: down_CREs           DMCs_inCREs       540      82        15        0.18292683 0.18292683      0.027777778            9699         0.055675843 3.914645e-05 7.829290e-05
#  7: down_CREs         DEGs_withCREs       327      82         7        0.08536585 0.08536585      0.021406728            9699         0.033714816 2.036627e-02 2.036627e-02
#  8: down_CREs DMCs_inCREs_with_DEGs        18      82         2        0.02439024 0.02439024      0.111111111            9699         0.001855861 9.895466e-03 1.319395e-02
#  9:   up_CREs           CpGs_inCREs      6034      15         8        0.53333333 0.53333333      0.001325820            9699         0.622125992 8.358412e-01 1.000000e+00
# 10:   up_CREs           DMCs_inCREs       540      15         1        0.06666667 0.06666667      0.001851852            9699         0.055675843 5.768065e-01 1.000000e+00
# 11:   up_CREs         DEGs_withCREs       327      15         0        0.00000000 0.00000000      0.000000000            9699         0.033714816 1.000000e+00 1.000000e+00
# 12:   up_CREs DMCs_inCREs_with_DEGs        18      15         0        0.00000000 0.00000000      0.000000000            9699         0.001855861 1.000000e+00 1.000000e+00
#     fold.enrichment     query
#  1:       1.4913941   da_CREs
#  2:       2.9626575   da_CREs
#  3:       2.1404521   da_CREs
#  4:      11.1099656   da_CREs
#  5:       1.6073914 down_CREs
#  6:       3.2855691 down_CREs
#  7:       2.5319982 down_CREs
#  8:      13.1422764 down_CREs
#  9:       0.8572754   up_CREs
# 10:       1.1974074   up_CREs
# 11:       0.0000000   up_CREs
# 12:       0.0000000   up_CREs

#end edit 23/11/22

#edit 27/11 
#DMCs/ DA CREs enrichment in DEGs
#term list : peaks DMCs/ up and dn
meth_da_cres_list<-list(CpGs=unique(peaks_methf$peaks),
                         DMCs=unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                         da_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&abs(avg_log2FC)>0.25]$peak),
                        down_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC<(-0.25)]$peak),
                         up_CREs=unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC>0.25]$peak),
                         DMCs_down_CREs=intersect(unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC<(-0.25)]$peak)),
                        DMCs_up_CREs=intersect(unique(peaks_methf[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(peaks_hsc_lga_xyf[p_val_adj<0.05&avg_log2FC>0.25]$peak)))
lapply(meth_da_cres_list, length)
# $CpGs
# [1] 7597
# 
# $DMCs
# [1] 651
# 
# $da_CREs
# [1] 102
# 
# $down_CREs
# [1] 85
# 
# $up_CREs
# [1] 17
# 
# $DMCs_down_CREs
# [1] 15
# 
# $DMCs_up_CREs
# [1] 1

#query list : CREs of DEGs
degs_cres_list<-list(DEGs=unique(pgla_e[padj<0.05&abs(log2FoldChange)>0.5]$peak_scatac),
                      upreg=unique(pgla_e[padj<0.05&log2FoldChange>0.5]$peak_scatac),
                      downreg=unique(pgla_e[padj<0.05&log2FoldChange<(-0.5)]$peak_scatac))
lapply(degs_cres_list, length)
# $DEGs
# [1] 327
# 
# $upreg
# [1] 96
# 
# $downreg
# [1] 231


#background : peaks in CREs

background<-unique(pgla_e$peak_scatac)
length(background) #10k peaks
head(background)

res_or_degs_cres_dmcs_da<-OR3(degs_cres_list,
                    meth_da_cres_list,
                    background = background,
                    overlap_column = F) 
res_or_degs_cres_dmcs_da

fwrite(res_or_degs_cres_dmcs_da,fp(out,"res_degs_cres_enrichment_for_dmcs_da_cres.csv"))




#end edit 27/11/22

#TF motif analysis on the peaks####
#Adding motif information to the Seurat object
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
cbls <- AddMotifs(
  object = cbls,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

motifs<-GetMotifIDs(object = cbls,motif.names = c("EGR1","KLF2","KLF4"))

motifs_peaks<-GetMotif(cbls,motifs = motifs,peaks=unique(res_peak_genes$peak))
motifs_peaks

res_motifpeak_genes<-merge(res_peak_genes,motifs_peaks,by="peak",all.x=T)
res_motifpeak_genes
res_motifpeak_genes[,motif.name:=as.character(motif.name)]

fwrite(res_motifpeak_genes,
       fp(out,"res_peaks_genes_links_EGR1_KLF2_KLF4_motifs_anno.csv"))

#EGR1 regulons genes with EGR1 motif on CRE
res_motifpeak_genes[,in_regulon:=gene%in%regulons[tf==motif.name[1]]$gene,by=motif.name]
res_motifpeak_genes[motif.name=="EGR1"&in_regulon==T]
intersect(res_motifpeak_genes[motif.name=="EGR1"]$gene,
          regulons[tf=="EGR1"]$gene)

#12/26, 46%
res_motifpeak_genes<-fread(fp(out,"res_peaks_genes_links_EGR1_KLF2_KLF4_motifs_anno.csv"))

res_motifpeak_genesf<-res_motifpeak_genes[!is.na(motif)]

#enrichement for him and others
DefaultAssay(cbls)<-"SCT"
res_enr<-OR3(querys = lapply(split(res_motifpeak_genesf$gene,res_motifpeak_genesf$motif.name),unique),
         terms_list = split(regulons$gene,regulons$tf),
         background = intersect(rownames(cbls),unique(regulons$gene)))

res_enr[padj<0.05][order(padj)]

res_enr[padj<0.05][query.=="EGR1"][order(-fold.enrichment)] #46%
res_enr[query.=="EGR1"&term=="EGR1"]
#    query. term term.size n.query n.overlap pct.query.overlap  precision
# 1:   EGR1 EGR1        26     714        12        0.01680672 0.01680672
#    pct.term.overlap background_size pct.term.background       pval      padj
# 1:        0.4615385            2475          0.01050505 0.04489891 0.1442207
#    fold.enrichment                                                    genes.overlap
# 1:        1.599871 PTCH2|JUN|TUBA1A|SOCS3|JUNB|IER2|ZFP36|FOSB|DUSP2|MAFF|GNL1|IER3
#    query
# 1:  EGR1


#KLF2
res_enr[padj<0.05][query.=="KLF2"][order(-fold.enrichment)] #44%
res_enr[query.=="KLF2"&term=="KLF2"]
#    query. term term.size n.query n.overlap pct.query.overlap  precision
# 1:   KLF2 KLF2        98     815        44        0.05398773 0.05398773
#    pct.term.overlap background_size pct.term.background        pval       padj
# 1:        0.4489796            2475          0.03959596 0.007815648 0.05523058
#    fold.enrichment
# 1:        1.363466

#KLF4
res_enr[padj<0.05][query.=="KLF4"][order(-fold.enrichment)] #38%
res_enr[query.=="KLF4"&term=="KLF4"]
#    query. term term.size n.query n.overlap pct.query.overlap  precision
# 1:   KLF4 KLF4        42     933        16        0.01714898 0.01714898
#    pct.term.overlap background_size pct.term.background      pval      padj
# 1:        0.3809524            2475           0.0169697 0.5374044 0.8224684
#    fold.enrichment
# 1:        1.010565

res_enr[,regulon:=term]
res_enr[,motif:=query.]
res_enr[,pct.regulon_with_CREsmotifs:=pct.term.overlap]
res_enr[,.(regulon,motif,pct.regulon_with_CREsmotifs,pval)][regulon%in%c("EGR1","KLF2","KLF4")][order(regulon)]

#    regulon motif pct.regulon_with_CREsmotifs        pval
# 1:    EGR1  EGR1                   0.4615385 0.044898906
# 5:    KLF2  KLF2                   0.4489796 0.007815648
# 9:    KLF4  KLF4                   0.3809524 0.537404443

#check pctage egr1 dans cres
unique(res_motifpeak_genesf[motif.name=="EGR1"],by="peak")
length(unique(res_peak_genes$peak))#1005/4015..  25% of CREs have EGR1
res_peak_genes
unique(res_motifpeak_genesf[motif.name=="KLF2"],by="peak")#1232
length(unique(res_motifpeak_genesf$peak))#1232/4015  30% of CREs have KLF2

unique(res_motifpeak_genesf[motif.name=="KLF4"],by="peak")#1552
length(unique(res_motifpeak_genesf$peak))#1552/4015  38% of CREs have KLF4

unique(res_motifpeak_genesf[motif.name%in%tfs],by="peak")#1826/4015 #45% of CREs with motifs

#merge with DMCs####
res_peak_genes[,chr:=seqnames]
res_peak_genes[,start.peak:=start(peak)]
res_peak_genes[,end.peak:=end(peak),by="peak"]

res_cpgs<-fread("outputs/14-DMCs_atac_integr/res_cpgs_hg38.cs.gz")
cpgs_in_cres<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          res_peak_genes[,.(chr,start.peak,end.peak,peak)][order(chr,start.peak,end.peak)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peak","start","end"))

length(unique(cpgs_in_cres$cpg_id))#17146/750k cpgs (instead of 7k for cbl12)
length(unique(cpgs_in_cres$peak)) #2.3k/4.2k peaks

fwrite(cpgs_in_cres,fp(out,"cpgs_in_CREs.csv.gz"))
#% DMCs overlapping CREs
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

cres_meth<-merge(cpgs_in_cres,res_meth,by="cpg_id")
nrow(cres_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]) #6% DMCs are in CREs 

#vs % CpGs overlapping 
nrow(cres_meth)/nrow(res_meth) #2.4% CpGs are in CREs

(nrow(cres_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]))/(nrow(cres_meth)/nrow(res_meth)) #2.47 foldenrichment


OR(set1 = res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id,
   set2= cres_meth$cpg_id,size_universe = nrow(res_meth)) #p = 2.033127e-37


#### CREs of OLD regulons EGRN integration with DMCs
tfs<-c("EGR1","KLF2","KLF4")
cpgs_in_cres_egr1n<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          res_peak_genes[gene%in%regulons[tf%in%tfs]$gene][,.(chr,start.peak,end.peak,peak)][order(chr,start.peak,end.peak)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peak","start","end"))

length(unique(cpgs_in_cres_egr1n$cpg_id))#1312/750k cpgs
length(unique(cpgs_in_cres_egr1n$peak)) #171/257 peaks

fwrite(cpgs_in_cres_egr1n,fp(out,"cpgs_in_CREs_EGR1_KLF2_KLF4_regulons.csv.gz"))

#% DMCs overlapping CREs EGR1_KLF2_KLF4_regulons
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

cres_meth<-merge(cpgs_in_cres,res_meth,by="cpg_id")
cres_meth[P.Value<0.001&abs(logFC)>25] #297
nrow(cres_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]) #2.5% DMCs are in CREs 

cres_meth[P.Value<0.001&abs(logFC)>25] #17 DMCs  in CRE 
cres_meth_genes<-merge(cres_meth,res_peak_genes,by=c("peak","chr"))
cres_meth_genes[P.Value<0.001&abs(logFC)>25]
cres_meth_genes[P.Value<0.001&abs(logFC)>25][,.(chr,peak,pos,P.Value,logFC,gene)]

#       chr                      peak       pos      P.Value    logFC    gene
#  1:  chr1  chr1-181088263-181090448 181089758 1.163360e-05 56.43245    IER5
#  2:  chr1    chr1-23553947-23555213  23554528 6.791300e-04 50.68169     ID3
#  3:  chr1    chr1-44813348-44814070  44813782 7.829432e-05 54.57994   PTCH2
#  4: chr10 chr10-102713406-102715152 102714631 8.426493e-05 51.52761   TRIM8
#  5: chr11     chr11-2883974-2886612   2885092 1.616649e-04 59.46079  CDKN1C
#  6: chr11     chr11-2883974-2886612   2885684 8.917287e-04 41.66439  CDKN1C
#  7: chr11     chr11-2883974-2886612   2885966 7.398161e-05 62.06919  CDKN1C
#  8: chr11   chr11-62545272-62547035  62545815 9.363807e-04 48.48676   AHNAK
#  9: chr11   chr11-62545272-62547035  62545848 7.401320e-04 42.24850   AHNAK
# 10: chr12   chr12-52050386-52051841  52050979 6.832596e-04 51.60495   GRASP
# 11: chr13   chr13-44575715-44578236  44577044 6.739804e-04 38.94290 TSC22D1
# 12: chr14   chr14-68795217-68796424  68795624 6.465818e-04 61.39083 ZFP36L1
# 13: chr15   chr15-31325512-31327254  31326011 5.526391e-04 45.65863   KLF13
# 14: chr15   chr15-84980271-84982396  84980535 4.367939e-04 51.58003  AKAP13
# 15: chr17   chr17-78262975-78263905  78263390 1.892195e-04 61.05853   SOCS3
# 16: chr17   chr17-78358820-78361021  78359375 3.759657e-04 65.01007   SOCS3
# 17: chr17   chr17-78358820-78361021  78359484 3.505579e-04 60.92096   SOCS3


#vs CpGs overlapping 
nrow(cres_meth) #1398 CpGs are in CREs EGRN


#ggplot(cres_meth)+geom_boxplot(aes(x=in_EGRN,y=-log10(P.Value)*logFC)) #nop



#DMCs overlapping new CREs EGR1/KLF2/KLF4 regulatory netw####
cpgs_in_cres<-fread(fp(out,"cpgs_in_CREs.csv.gz"))

genes_egrncres<-unique(res_motifpeak_genesf[motif.name%in%tfs][,.(gene,motif.name,peak)])
cres_cpgs_egrn<-merge(cpgs_in_cres,genes_egrncres,by="peak",,allow.cartesian=TRUE)
cres_meth_egrn<-merge(cres_cpgs_egrn,res_meth,by="cpg_id")

unique(cpgs_in_cres,by="cpg_id") #17k

unique(cres_cpgs_egrn,by="cpg_id") #15k => 15/17, 88% have EGR1 or KLF2/KLF4 motifs

cres_meth[P.Value<0.001&abs(logFC)>25]
unique(cres_meth[P.Value<0.001&abs(logFC)>25],by="peak")
unique(cres_meth_egrn[P.Value<0.001&abs(logFC)>25],by="cpg_id") #243/4815 5% DMCs in CREs EGRN
unique(cres_meth_egrn[P.Value<0.001&abs(logFC)>25],by="peak") #198/216 91% of DMCs asso CREs in CREs EGRN
#243/296 #82%
res_meth[P.Value<0.001&abs(logFC)>25]#4815 DMCs

#vs CpGs overlapping 
nrow(unique(cres_meth_egrn,by="cpg_id")) #15179/750000 2%CpGs are in CREs EGRN
res_meth

#enrichment DMCs in CREs egrn ?
OR(set1 = res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id,
   set2= cres_meth_egrn$cpg_id,size_universe = nrow(res_meth)) #p = 4.883341e-37

#enrichment DMCs in CREs egrn compared to DMCs other CREs ?
univ<-unique(cres_meth$cpg_id)
OR(set1 = intersect(res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id,univ),
   set2= cres_meth_egrn$cpg_id,size_universe = length(univ)) #p = 0.03

#enrichment  for egrn CREs in DMCs asso CREs ?
univ<-unique(cres_meth$peak)
OR(set1 = unique(cres_meth[P.Value<0.001&abs(logFC)>25]$peak),
   set2= unique(cres_meth_egrn$peak),size_universe = length(univ)) #p = 5.710426e-20

intersect(cres_meth_egrn$peak,unique(cres_meth[P.Value<0.001&abs(logFC)>25]$peak))
#1532/2311(66%) peak asso with cpgs have EGR1/klf2/klf4. 198/216(91%) peaks with DMCs have EGR1/klf2/klf4
#increase Meth in CREs EGRN vs CREs ?
cres_meth[,in_EGRN:=cpg_id%in%cres_meth_egrn$cpg_id]
unique(cres_meth[(in_EGRN)],by="cpg_id")
unique(cres_meth,by="cpg_id") #15/17 #88% of CREsCpGs have EGR1/klf2/klf4

cres_meth[,dmc_score:=-log10(P.Value)*logFC]
plot(density(unique(cres_meth,by="cpg_id")$dmc_score))
ggplot(unique(cres_meth,by="cpg_id"))+geom_boxplot(aes(x=in_EGRN,y=-log10(P.Value)*logFC)) #un peu
#ggplot(unique(cres_meth,by="cpg_id"))+geom_bar(aes(x=in_EGRN,fill=P.Value<0.001&abs(logFC)>25),position="fill") #un peu

cres_meth[,FoldChange:=mean(dmc_score[in_EGRN==T])/mean(dmc_score[in_EGRN==F])] #1.39

cres_meth[,avg_log2FC:=log2(FoldChange)] #0.47

cres_meth[,pval:=wilcox.test(dmc_score[in_EGRN==T],
                             dmc_score[in_EGRN==F])$p.value] #p = 3.407891e-30


#CREs egrn of the genes regulons
cpgs_in_cres_egr1n<-fread(fp(out,"cpgs_in_CREs_EGR1_KLF2_KLF4_regulons.csv.gz"))
unique(cpgs_in_cres,by="cpg_id") #17k

unique(cpgs_in_cres_egr1n,by="cpg_id") #1312 
meth_in_cres_egr1n<-merge(cpgs_in_cres_egr1n,res_meth,by="cpg_id")

cres_meth[P.Value<0.001&abs(logFC)>25]
unique(cres_meth[P.Value<0.001&abs(logFC)>25],by="peak")
unique(meth_in_cres_egr1n[P.Value<0.001&abs(logFC)>25],by="cpg_id") #17/4815  DMCs in CREs EGRN regulons
unique(meth_in_cres_egr1n[P.Value<0.001&abs(logFC)>25],by="peak") #13/216 6% of DMCs asso CREs in CREs EGRN regulons
#243/296 #82%
res_meth[P.Value<0.001&abs(logFC)>25]#4815 DMCs
#enrichment  for egrn CREs in DMCs asso CREs ?
univ<-unique(cres_meth$peak)
OR(set1 = unique(cres_meth[P.Value<0.001&abs(logFC)>25]$peak),
   set2= unique(meth_in_cres_egr1n$peak),size_universe = length(univ)) #p =0.89

#EFFECT OF LGA####
##no stim####
unstim<-subset(cbls,condition=="no_stim")
head(colnames(unstim))
#demultiplex based on HTO
hto_mat<-Read10X("~/RUN/Run_733_single-cell/output/count/multi-a/outs/raw_feature_bc_matrix/")$`Antibody Capture`
colnames(hto_mat)<-paste0(colnames(hto_mat),"_1")
cells<-intersect(colnames(unstim),colnames(hto_mat))
length(cells) #1986
ncol(unstim) #2001
unstim<-unstim[,cells]
unstim[["HTO"]]<-CreateAssayObject(as.matrix(hto_mat)[,cells])
rm(hto_mat)

unstim <- NormalizeData(unstim, assay = "HTO", normalization.method = "CLR")
unstim <- HTODemux(unstim, assay = "HTO",positive.quantile = 0.95)
table(unstim$HTO_classification.global)
# Doublet Negative  Singlet 
#       33     1109      844 

table(unstim$hash.ID) #negative = LGA ?
#yes if difference negative / CTRL +  ++ sig
hto_dt<-merge(data.table(unstim@meta.data,keep.rownames = "cell"),melt(data.table(unstim@assays$HTO@data,keep.rownames = "HTO"),id.vars="HTO",variable.name="cell",value.name = "expression"))
ggplot(hto_dt[HTO=="multi-a-CTRL-no-stim"])+geom_density(aes(x=expression,fill=hash.ID),alpha=0.4)

ggplot(hto_dt[HTO=="multi-a-CTRL-no-stim"])+geom_density(aes(x=expression),alpha=0.4)

hto_dt[,median_pos:=median(expression),by=c("HTO","hash.ID")]
med_pos_ctrl<-hto_dt[hash.ID=="multi-a-CTRL-no-stim"&HTO=="multi-a-CTRL-no-stim"]$median_pos[1]
hto_dt[,FC_pos:=expression/med_pos_ctrl]

ggplot(hto_dt[HTO=="multi-a-CTRL-no-stim"])+
  geom_density(aes(x=FC_pos))

#LGA recovery = all with expr HTO CTRL <q75 LGA cells
lga_reco<-hto_dt[HTO=="multi-a-CTRL-no-stim"][expression<quantile(expression[hash.ID=="multi-a-LGA-no-stim"],0.75)]$cell
length(lga_reco)#747
unstim[["group"]]<-sapply(colnames(unstim), function(c){
  x<-unstim@meta.data[c,"hash.ID"]
  return(ifelse(x=="Negative"&c%in%lga_reco,"LGA",str_extract(x,"CTRL|LGA|Doublet|Negative")))
})
table(unstim$group)

    # CTRL  Doublet      LGA Negative 
    #  814       33      756      383 
#save 
saveRDS(unstim, fp(out,"cbls_unstim_hto_anno.rds"))
fwrite(data.table(unstim@meta.data,keep.rownames="bc"), fp(out,"metadata_cbls_unstim_hto_anno.csv.gz"))

#DEGs
#all lin
res<-data.table(FindMarkers(unstim,group.by="group",ident.1="LGA",ident.2 = "CTRL"),keep.rownames = "gene")
res[p_val_adj<0.05] #increase NFKB1, LTB, KDM5D (Lysine Demethylase 5D, demethyling H3K4) decrease ERG,
#38 DEGS, with ++ sex transcrits( increase male in LGA)
fwrite(res,fp(out,"res_hspc_no_stim_lga_vs_ctrl.csv"))
#hsc
table(unstim@meta.data$group,unstim@meta.data$ident.2) #117 CTRL, 194 LGA
Idents(unstim)<-"ident.2"
res_hsc<-data.table(FindMarkers(unstim,group.by="group",ident.1="LGA",ident.2 = "CTRL",subset.ident = 'HSC'),keep.rownames = "gene")
res_hsc[p_val_adj<0.05] #no degs except XIST
fwrite(res_hsc,fp(out,"res_hsc_no_stim_lga_vs_ctrl.csv"))

#DAP hsc
DefaultAssay(unstim)<-"peaks"
metadata <- read.csv(
  file = "~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-a/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
rownames(metadata)<-paste0(rownames(metadata),"_1")
unstim<-AddMetaData(unstim,metadata = metadata[colnames(unstim),])
head(unstim@meta.data)
da_peaks <- FindMarkers(
  object = unstim,
  group.by = "group",
  ident.1 = "LGA", 
  ident.2 ="CTRL",
  subset.ident = "HSC",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments'
)
da_peaks<-data.table(da_peaks,keep.rownames = "peak")
da_peaks[p_val_adj<0.05]
da_peaks[p_val<0.001]

#which genes ?
da_peaks_genes<-merge(da_peaks,res_peak_genes,all.x=T)
da_peaks_genes#only 10/339 link to gene

fwrite(da_peaks_genes[order(p_val)],fp(out,"res_da_peaks_lga_vs_ctrl_hsc_no_stim.csv"))
dapg<-fread(fp(out,"res_da_peaks_lga_vs_ctrl_hsc_no_stim.csv"))
dapg[!is.na(start)]

##stim2h####
stim2<-subset(cbls,condition=="stim_2h")
head(colnames(stim2))
#demultiplex based on HTO
hto_mat<-Read10X("~/RUN/Run_733_single-cell/output/count/multi-b/outs/raw_feature_bc_matrix/")$`Antibody Capture`
colnames(hto_mat)<-paste0(colnames(hto_mat),"_2")
cells<-intersect(colnames(stim2),colnames(hto_mat))
length(cells) #2339
ncol(stim2) #2353
stim2<-stim2[,cells]
stim2[["HTO"]]<-CreateAssayObject(as.matrix(hto_mat)[,cells])
rm(hto_mat)

stim2 <- NormalizeData(stim2, assay = "HTO", normalization.method = "CLR")
stim2 <- HTODemux(stim2, assay = "HTO",positive.quantile = 0.95)
table(stim2$HTO_classification.global)
            # Doublet             Negative  multi-b-LGA-stim-2h multi-b-CTRL-stim-2h 
            #       37                 1363                  911                   28 

table(stim2$hash.ID) #negative = CTRL ?
#yes if difference negative / CTRL +  ++ sig
hto_dt<-merge(data.table(stim2@meta.data,keep.rownames = "cell"),melt(data.table(stim2@assays$HTO@data,keep.rownames = "HTO"),id.vars="HTO",variable.name="cell",value.name = "expression"))
ggplot(hto_dt[HTO=="multi-b-LGA-stim-2h"])+geom_density(aes(x=expression,fill=hash.ID),alpha=0.4)

ggplot(hto_dt[HTO=="multi-b-LGA-stim-2h"])+geom_density(aes(x=expression),alpha=0.4)

hto_dt[,median_pos:=median(expression),by=c("HTO","hash.ID")]
med_pos_lga<-hto_dt[hash.ID=="multi-b-LGA-stim-2h"&HTO=="multi-b-LGA-stim-2h"]$median_pos[1]
hto_dt[,FC_pos:=expression/med_pos_lga]

ggplot(hto_dt[HTO=="multi-b-LGA-stim-2h"])+
  geom_density(aes(x=FC_pos))

#CTRl recovery = all with expr HTO CTRL <q75 LGA cells
ctrl_reco<-hto_dt[HTO=="multi-b-LGA-stim-2h"][expression<quantile(expression[hash.ID=="multi-b-CTRL-stim-2h"],0.75)]$cell
length(ctrl_reco)#1014
stim2[["group"]]<-sapply(colnames(stim2), function(c){
  x<-stim2@meta.data[c,"hash.ID"]
  return(ifelse(x=="Negative"&c%in%ctrl_reco,"CTRL",str_extract(x,"CTRL|LGA|Doublet|Negative")))
})
table(stim2$group)
# 
#       CTRL  Doublet      LGA Negative 
#     1021       37      911      370 
#save 
saveRDS(stim2, fp(out,"cbls_stim2h_hto_anno.rds"))
fwrite(data.table(stim2@meta.data,keep.rownames="bc"), fp(out,"metadata_cbls_stim2h_hto_anno.csv.gz"))
  
#DEGs
#all lin
res<-data.table(FindMarkers(stim2,group.by="group",ident.1="LGA",ident.2 = "CTRL"),keep.rownames = "gene")
res[p_val_adj<0.05] #increase MT gene, NFKB1, IL2RA, LTB, 
#98 DEGS, with ++ sex transcrits( increase male in LGA)
fwrite(res,fp(out,"res_hspc_stim2h_lga_vs_ctrl.csv"))
#hsc
table(stim2@meta.data$group,stim2@meta.data$ident.2) #94 CTRL, 130 LGA
Idents(stim2)<-"ident.2"
res_hsc<-data.table(FindMarkers(stim2,group.by="group",ident.1="LGA",ident.2 = "CTRL",subset.ident = 'HSC'),keep.rownames = "gene")
res_hsc[p_val_adj<0.05] #MT genes upreg
res_hsc[p_val<0.001] #MT genes upreg
#         gene        p_val avg_log2FC pct.1 pct.2    p_val_adj
#  1:     XIST 1.744479e-10 -1.3157759 0.223 0.585 3.936940e-06
#  2:   MT-CYB 1.526677e-08  0.7346400 1.000 1.000 3.445406e-04
#  3:   MT-ND3 3.885416e-08  0.7937510 1.000 1.000 8.768606e-04
#  4:  MT-ATP6 5.650237e-08  0.7281760 1.000 1.000 1.275145e-03
#  5:      UTY 7.147135e-08  0.7114239 0.877 0.447 1.612965e-03
#  6:   MT-CO2 2.764339e-07  0.7866396 1.000 1.000 6.238560e-03
#  7:   MT-ND2 3.193808e-07  0.7363339 1.000 0.989 7.207787e-03
#  8:   MT-CO3 7.346962e-07  0.5781713 1.000 1.000 1.658062e-02
#  9:  ANKRD28 2.393401e-06 -0.7316404 0.892 0.979 5.401427e-02
# 10:    USP9Y 3.594188e-06  0.5234686 0.746 0.383 8.111363e-02
# 11:   VPS13D 4.721657e-06  0.4992228 0.831 0.585 1.065584e-01
# 12:   MT-CO1 1.413581e-05  0.6004498 1.000 1.000 3.190170e-01
# 13:      SPN 1.791290e-05  0.5532827 0.585 0.319 4.042583e-01
# 14:    NEAT1 1.939646e-05 -0.5203589 0.785 0.926 4.377392e-01
# 15:   AGPAT4 2.275412e-05 -0.5774035 0.331 0.606 5.135149e-01
# 16:     ULK4 2.691135e-05 -0.7342443 0.738 0.840 6.073354e-01
# 17:  MT-ND4L 9.017518e-05  0.5851468 0.923 0.745 1.000000e+00
# 18:      FOS 1.052548e-04 -0.9257023 0.177 0.394 1.000000e+00
# 19:   MT-ND1 1.243516e-04  0.5676264 1.000 0.989 1.000000e+00
# 20:     ZEB1 1.333595e-04 -0.4285373 0.869 0.957 1.000000e+00
# 21:  MT-ATP8 1.686675e-04  0.4867637 0.954 0.894 1.000000e+00
# 22: SEPTIN11 2.308526e-04 -0.4677790 0.569 0.745 1.000000e+00
# 23:     GRB2 2.361833e-04  0.4930504 0.677 0.479 1.000000e+00 
# 24:     FOSB 2.615998e-04 -0.6830469 0.177 0.394 1.000000e+00
# 25:    MYADM 2.804907e-04 -0.5477473 0.292 0.511 1.000000e+00
# 26:    PELI1 3.031176e-04 -0.7288248 0.431 0.606 1.000000e+00
# 27:    STAP1 3.212704e-04  0.4483894 0.454 0.245 1.000000e+00
# 28:     ATF7 3.466430e-04 -0.4677790 0.292 0.479 1.000000e+00
# 29:   MT-ND4 3.869307e-04  0.4919673 1.000 1.000 1.000000e+00
# 30:     WHRN 3.871278e-04 -0.2692333 0.062 0.223 1.000000e+00
# 31:   HECTD1 4.014910e-04  0.4003379 0.923 0.734 1.000000e+00
# 32:    ZNF26 5.713949e-04  0.3464296 0.446 0.234 1.000000e+00
# 33:   ELOVL5 6.007681e-04 -0.3806464 0.546 0.745 1.000000e+00
# 34:   EIF1AY 6.220060e-04  0.2635730 0.254 0.074 1.000000e+00
# 35:    BUD13 6.602617e-04 -0.3067165 0.085 0.245 1.000000e+00
# 36:    UTP20 6.604013e-04  0.3878311 0.446 0.234 1.000000e+00
# 37:    TANC1 6.777314e-04 -0.3808433 0.177 0.372 1.000000e+00
# 38:    KDM6A 7.105526e-04 -0.4094803 0.938 0.947 1.000000e+00
# 39:    SIAH1 7.146000e-04 -0.3608638 0.223 0.426 1.000000e+00
# 40:   NFKBIA 8.455195e-04 -0.6187209 0.769 0.851 1.000000e+00
# 41:   ANGPT1 8.691816e-04  0.7117049 0.408 0.213 1.000000e+00
# 42:     PIM3 9.354585e-04 -0.4742340 0.438 0.606 1.000000e+00
# 43:      VIM 9.665343e-04 -0.4401005 0.823 0.872 1.000000e+00
fwrite(res_hsc,fp(out,"res_hsc_stim2h_lga_vs_ctrl.csv"))
res_hsc<-fread(fp(out,"res_hsc_stim2h_lga_vs_ctrl.csv"))
res_hsc[p_val<0.001&gene%in%genes_egrncres$gene]
#        gene        p_val avg_log2FC pct.1 pct.2  p_val_adj
#  1: ANKRD28 2.393401e-06 -0.7316404 0.892 0.979 0.05401427
#  2:     SPN 1.791290e-05  0.5532827 0.585 0.319 0.40425830
#  3:   NEAT1 1.939646e-05 -0.5203589 0.785 0.926 0.43773924
#  4:     FOS 1.052548e-04 -0.9257023 0.177 0.394 1.00000000
#  5:    ZEB1 1.333595e-04 -0.4285373 0.869 0.957 1.00000000
#  6:    FOSB 2.615998e-04 -0.6830469 0.177 0.394 1.00000000
#  7:   KDM6A 7.105526e-04 -0.4094803 0.938 0.947 1.00000000
#  8:  NFKBIA 8.455195e-04 -0.6187209 0.769 0.851 1.00000000
#  9:    PIM3 9.354585e-04 -0.4742340 0.438 0.606 1.00000000
# 10:     VIM 9.665343e-04 -0.4401005 0.823 0.872 1.00000000

#DAP hsc
DefaultAssay(stim2)<-"peaks"
metadata <- read.csv(
  file = "~/RUN/Run_735_single-cell/output/count_multi_gex_atac/multi-b/outs/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)
rownames(metadata)<-paste0(rownames(metadata),"_2")
stim2<-AddMetaData(stim2,metadata = metadata[colnames(stim2),])
head(stim2@meta.data)
da_peaks <- FindMarkers(
  object = stim2,
  group.by = "group",
  ident.1 = "LGA", 
  ident.2 ="CTRL",
  subset.ident = "HSC",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments'
)
da_peaks<-data.table(da_peaks,keep.rownames = "peak")
da_peaks[p_val_adj<0.05]
da_peaks[p_val<0.001]

#which genes ?
da_peaks_genes<-merge(da_peaks,res_peak_genes,all.x=T)
da_peaks_genes[!is.na(gene)]# 44/861 link to gene
#                        peak        p_val avg_log2FC pct.1 pct.2
#  1: chr1-109669671-109669904 4.173117e-03  0.2718667 0.054 0.000
#  2: chr1-196608282-196609512 1.469177e-02 -0.2561756 0.146 0.287
#  3: chr1-198934794-198935909 9.699090e-03  0.2626479 0.215 0.096
#  4:   chr1-24997470-24998092 1.102516e-02 -0.2558516 0.038 0.117
#  5:   chr1-44422196-44422480 6.255123e-05 -0.3422138 0.000 0.096
#  6:    chr10-3757619-3758187 1.595676e-02 -0.3169960 0.015 0.074
#  7:  chr10-72086660-72089123 1.852298e-03 -0.3267952 0.062 0.202
#  8:    chr12-6766098-6767260 5.561893e-03 -0.2751431 0.146 0.298
#  9:  chr12-76878095-76879787 3.367903e-03 -0.2908373 0.192 0.383
# 10:  chr12-89345766-89347105 3.412519e-04  0.3083078 0.231 0.074
# 11:  chr13-23378113-23378545 2.773248e-03  0.3089236 0.131 0.043
# 12:  chr14-35375533-35376076 1.689309e-02  0.2893651 0.085 0.021
# 13:  chr16-10846461-10846998 3.966496e-02 -0.3132519 0.069 0.160
# 14:  chr16-11152943-11153378 3.734860e-02 -0.2942764 0.031 0.096
# 15:  chr16-11152943-11153378 3.734860e-02 -0.2942764 0.031 0.096
# 16:    chr16-4315249-4316616 2.400508e-03 -0.3277656 0.231 0.394
# 17:  chr16-67846304-67847635 1.209741e-05 -0.3235820 0.377 0.628
# 18:  chr16-88455188-88456392 1.531352e-03 -0.2514487 0.431 0.574
# 19:    chr17-1561648-1563082 7.713141e-04 -0.3120462 0.177 0.372
# 20:  chr19-18225621-18226494 3.613614e-03 -0.3516868 0.115 0.245
# 21:  chr19-42596195-42596807 2.241625e-02 -0.2558057 0.023 0.096
# 22:  chr19-45442780-45444572 8.562738e-03 -0.2658373 0.169 0.309
# 23:    chr20-4665733-4666501 2.612731e-02 -0.2850978 0.169 0.298
# 24:  chr22-37538617-37539174 2.138413e-02 -0.2513723 0.069 0.170
# 25:  chr22-50190079-50191772 5.650239e-04 -0.2900150 0.369 0.511
# 26:  chr22-50190079-50191772 5.650239e-04 -0.2900150 0.369 0.511
# 27: chr3-167379782-167380518 3.021301e-02 -0.2558430 0.123 0.245
# 28:   chr3-98778226-98779455 6.791919e-04 -0.2824181 0.008 0.106
# 29:     chr4-4859287-4860197 8.943093e-03 -0.2663502 0.062 0.160
# 30: chr5-142108344-142109587 1.475812e-03 -0.2986766 0.138 0.298
# 31:   chr5-60677549-60678017 4.515451e-03 -0.3114053 0.008 0.074
# 32:   chr5-66821660-66822421 1.097279e-01 -0.2703246 0.054 0.106
# 33: chr6-116392274-116393025 1.162650e-02 -0.3687447 0.031 0.117
# 34: chr6-137873227-137873702 9.861721e-03 -0.3246408 0.023 0.096
# 35: chr6-137939709-137940535 1.313666e-02  0.2645664 0.077 0.011
# 36:   chr6-26025717-26027480 3.154363e-03 -0.2800972 0.146 0.298
# 37:   chr6-26025717-26027480 3.154363e-03 -0.2800972 0.146 0.298
# 38:   chr6-31582099-31582581 2.083149e-02  0.2920471 0.162 0.074
# 39:   chr6-31582099-31582581 2.083149e-02  0.2920471 0.162 0.074
# 40:   chr6-32940717-32941738 1.090455e-02 -0.2828870 0.108 0.223
# 41:   chr6-36682530-36683148 4.619724e-02 -0.2587880 0.077 0.170
# 42:   chr6-42742426-42742819 2.879554e-03  0.2879844 0.123 0.021
# 43:   chr8-67028063-67028956 8.657794e-03 -0.2980655 0.069 0.191
# 44:   chr8-80366775-80367982 2.612219e-03  0.2903173 0.062 0.000
#                         peak        p_val avg_log2FC pct.1 pct.2
#     p_val_adj seqnames     start       end  width strand      score
#  1:         1     chr1 109616104 109669788  53685      * 0.05486597
#  2:         1     chr1 196608897 196651878  42982      * 0.05264962
#  3:         1     chr1 198638671 198935352 296682      * 0.05321227
#  4:         1     chr1  24965121  24997781  32661      * 0.08709927
#  5:         1     chr1  44405194  44422338  17145      * 0.07783265
#  6:         1    chr10   3757903   3785281  27379      * 0.06420560
#  7:         1    chr10  72087892  72089032   1141      * 0.05313932
#  8:         1    chr12   6766679   6867119 100441      * 0.05479687
#  9:         1    chr12  76559809  76878941 319133      * 0.07277046
# 10:         1    chr12  89346436  89353271   6836      * 0.05868408
# 11:         1    chr13  23378329  23411513  33185      * 0.06872973
# 12:         1    chr14  35375805  35404749  28945      * 0.09357242
# 13:         1    chr16  10846730  11256179 409450      * 0.06666609
# 14:         1    chr16  11153161  11256179 103019      * 0.06795666
# 15:         1    chr16  11153161  11636381 483221      * 0.05747239
# 16:         1    chr16   4315933   4425705 109773      * 0.05376017
# 17:         1    chr16  67562407  67846970 284564      * 0.05734457
# 18:         1    chr16  88453317  88455790   2474      * 0.06026599
# 19:         1    chr17   1562365   1628886  66522      * 0.06274285
# 20:         1    chr19  18001132  18226058 224927      * 0.06828800
# 21:         1    chr19  42268537  42596501 327965      * 0.05161902
# 22:         1    chr19  45443676  45467995  24320      * 0.05032910
# 23:         1    chr20   4666117   4686236  20120      * 0.07353515
# 24:         1    chr22  37282027  37538896 256870      * 0.07071337
# 25:         1    chr22  50185915  50190926   5012      * 0.06279094
# 26:         1    chr22  50190926  50307627 116702      * 0.08237722
# 27:         1     chr3 167380150 167653983 273834      * 0.06065837
# 28:         1     chr3  98732236  98778841  46606      * 0.06844178
# 29:         1     chr4   4859742   5019472 159731      * 0.05016340
# 30:         1     chr5 142108505 142108966    462      * 0.05307248
# 31:         1     chr5  60522120  60677783 155664      * 0.05833789
# 32:         1     chr5  66596361  66822041 225681      * 0.10035792
# 33:         1     chr6 116254173 116392650 138478      * 0.08355154
# 34:         1     chr6 137867188 137873465   6278      * 0.10285895
# 35:         1     chr6 137867188 137940122  72935      * 0.06694856
# 36:         1     chr6  25962802  26026599  63798      * 0.05064728
# 37:         1     chr6  26026599  26234933 208335      * 0.06471628
# 38:         1     chr6  31582340  31582522    183      * 0.09535153
# 39:         1     chr6  31582340  31586124   3785      * 0.06267523
# 40:         1     chr6  32668383  32941228 272846      * 0.05102272
# 41:         1     chr6  36594353  36682839  88487      * 0.06034137
# 42:         1     chr6  42742623  42746096   3474      * 0.05592430
# 43:         1     chr8  67028510  67343677 315168      * 0.05069314
# 44:         1     chr8  80367379  80485619 118241      * 0.05783877
#     p_val_adj seqnames     start       end  width strand      score
#         gene   zscore       pvalue   chr start.peak  end.peak
#  1:    AMPD2 3.231849 6.149600e-04  chr1  109669671 109669904
#  2:      CFH 3.093905 9.877044e-04  chr1  196608282 196609512
#  3:    PTPRC 2.221299 1.316535e-02  chr1  198934794 198935909
#  4:    RUNX3 3.961395 3.725651e-05  chr1   24997470  24998092
#  5:   RNF220 2.046988 2.032961e-02  chr1   44422196  44422480
#  6:     KLF6 4.214535 1.251464e-05 chr10    3757619   3758187
#  7:   SPOCK2 2.942674 1.626953e-03 chr10   72086660  72089123
#  8:     TPI1 3.343286 4.139621e-04 chr12    6766098   6767260
#  9:   OSBPL8 2.246333 1.234134e-02 chr12   76878095  76879787
# 10:    DUSP6 2.634789 4.209477e-03 chr12   89345766  89347105
# 11:     SACS 3.035232 1.201755e-03 chr13   23378113  23378545
# 12:   NFKBIA 4.049753 2.563581e-05 chr14   35375533  35376076
# 13:    SOCS1 5.173338 1.149741e-07 chr16   10846461  10846998
# 14:    SOCS1 5.453908 2.463733e-08 chr16   11152943  11153378
# 15:    LITAF 3.584499 1.688629e-04 chr16   11152943  11153378
# 16:    CORO7 2.097989 1.795305e-02 chr16    4315249   4316616
# 17:     CTCF 2.075247 1.898185e-02 chr16   67846304  67847635
# 18:    ZFPM1 2.530140 5.700849e-03 chr16   88455188  88456392
# 19:  SLC43A2 2.232157 1.280229e-02 chr17    1561648   1563082
# 20:   ARRDC2 1.669204 4.753846e-02 chr19   18225621  18226494
# 21:      CIC 2.670169 3.790651e-03 chr19   42596195  42596807
# 22:     FOSB 3.507549 2.261279e-04 chr19   45442780  45444572
# 23:     PRNP 4.380823 5.911603e-06 chr20    4665733   4666501
# 24:    CYTH4 1.782051 3.737047e-02 chr22   37538617  37539174
# 25:    TRABD 2.986937 1.408939e-03 chr22   50190079  50191772
# 26:   PLXNB2 2.866370 2.076046e-03 chr22   50190079  50191772
# 27:    WDR49 2.274836 1.145789e-02  chr3  167379782 167380518
# 28:  ST3GAL6 2.281665 1.125456e-02  chr3   98778226  98779455
# 29:    CYTL1 1.766447 3.866041e-02  chr4    4859287   4860197
# 30:   NDFIP1 2.198716 1.394908e-02  chr5  142108344 142109587
# 31:    PDE4D 3.334035 4.279801e-04  chr5   60677549  60678017
# 32:    MAST4 6.000186 9.854595e-10  chr5   66821660  66822421
# 33:      DSE 4.365604 6.338599e-06  chr6  116392274 116393025
# 34:  TNFAIP3 6.220743 2.474026e-10  chr6  137873227 137873702
# 35:  TNFAIP3 3.486649 2.445563e-04  chr6  137939709 137940535
# 36:   TRIM38 2.247150 1.231522e-02  chr6   26025717  26027480
# 37: HIST1H1D 4.254313 1.048460e-05  chr6   26025717  26027480
# 38:      LTB 6.012863 9.113772e-10  chr6   31582099  31582581
# 39:     LST1 3.935647 4.148648e-05  chr6   31582099  31582581
# 40: HLA-DQB1 3.086569 1.012404e-03  chr6   32940717  32941738
# 41:    SRSF3 3.250356 5.763029e-04  chr6   36682530  36683148
# 42:     TBCC 3.818441 6.714894e-05  chr6   42742426  42742819
# 43:  ARFGEF1 1.887280 2.956133e-02  chr8   67028063  67028956
# 44:   ZBTB10 3.262040 5.530675e-04  chr8   80366775  80367982
#         gene   zscore       pvalue   chr start.peak  end.peak


fwrite(da_peaks_genes[order(p_val)],fp(out,"res_da_peaks_lga_vs_ctrl_hsc_2h_stim.csv"))

dapg<-fread(fp(out,"res_da_peaks_lga_vs_ctrl_hsc_2h_stim.csv"))
dapg[p_val<0.001]
dapg[!is.na(start)][order(p_val)]

#integr with res_degs
res_all<-merge(res_hsc,dapg,by="gene")

res_all
#      gene      p_val.x avg_log2FC.x pct.1.x pct.2.x p_val_adj.x
# 1:   FOSB 0.0002615998   -0.6830469   0.177   0.394           1
# 2:   KLF6 0.0166917825   -0.3977153   0.785   0.809           1
# 3: NFKBIA 0.0008455195   -0.6187209   0.769   0.851           1
# 4: OSBPL8 0.0186909956   -0.2988500   0.838   0.883           1
# 5:  WDR49 0.2044125457   -0.3718099   0.246   0.298           1
#                        peak     p_val.y avg_log2FC.y pct.1.y
# 1:  chr19-45442780-45444572 0.008562738   -0.2658373   0.169
# 2:    chr10-3757619-3758187 0.015956760   -0.3169960   0.015
# 3:  chr14-35375533-35376076 0.016893090    0.2893651   0.085
# 4:  chr12-76878095-76879787 0.003367903   -0.2908373   0.192
# 5: chr3-167379782-167380518 0.030213007   -0.2558430   0.123
#    pct.2.y p_val_adj.y seqnames     start       end  width strand
# 1:   0.309           1    chr19  45443676  45467995  24320      *
# 2:   0.074           1    chr10   3757903   3785281  27379      *
# 3:   0.021           1    chr14  35375805  35404749  28945      *
# 4:   0.383           1    chr12  76559809  76878941 319133      *
# 5:   0.245           1     chr3 167380150 167653983 273834      *
#         score   zscore       pvalue   chr start.peak  end.peak
# 1: 0.05032910 3.507549 2.261279e-04 chr19   45442780  45444572
# 2: 0.06420560 4.214535 1.251464e-05 chr10    3757619   3758187
# 3: 0.09357242 4.049753 2.563581e-05 chr14   35375533  35376076
# 4: 0.07277046 2.246333 1.234134e-02 chr12   76878095  76879787
# 5: 0.06065837 2.274836 1.145789e-02  chr3  167379782 167380518
#with meth ?
cres_meth[,in_egrn_regulons:=peak%in%meth_in_cres_egr1n$peak]
unique(cres_meth,by="cpg_id") #17k

unique(cres_meth[(in_egrn_regulons)],by="cpg_id") #1293/17k
res_all_meth<-merge(res_all,cres_meth,by="peak")
res_all_meth[P.Value<0.01]

#merge all cbls ####
cbls<-readRDS("outputs/27-multiomics_cbl_stim/cbls_stim.rds")
p1<-DimPlot(cbls,label=T,group.by = "ident.2")

p2<-DimPlot(cbls,label=T,reduction = "multi_humap",group.by = "ident.2")

p1+p2

FeaturePlot(cbls,c("LTB","CD99","VPREB1","IGHM"),label=T)
DimPlot(cbls,label=T,group.by = "condition")

#add HTO infos
mtdhto<-rbind(unstim@meta.data,stim2@meta.data)
fwrite(data.table(mtdhto,keep.rownames="bc"),fp(out,"all_metadata_unstim_stim2h.csv.gz"))
mtdhto<-fread(fp(out,"all_metadata_unstim_stim2h.csv.ga"))
mtdhto<-data.frame(mtdhto,row.names = "bc")
cbls$group<-mtdhto[colnames(cbls),"group"]
cbls$atac_peak_region_fragments<-mtdhto[colnames(cbls),"atac_peak_region_fragments"]

#add motifs and save
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
DefaultAssay(cbls)<-"peaks"
cbls <- AddMotifs(
  object = cbls,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
#save
saveRDS(cbls,fp(out,"cbls_stim_final.rds"))
table(cbls$group,cbls$condition)
  #         no_stim stim_2h
  # CTRL         814    1021
  # Doublet       33      37
  # LGA          756     911
  # Negative     383     370

table(cbls$group[cbls$ident.2=="HSC"],cbls$condition[cbls$ident.2=="HSC"])
  #          no_stim stim_2h
  # CTRL         117      94
  # Doublet        4       4
  # LGA          194     130
  # Negative      83      43



#tf-target netw based on regulons scenic + tf-peak-gene multimod ####
cbls<-readRDS("outputs/27-multiomics_cbl_stim/cbls_stim_final.rds")


regulons_ext<-fread("outputs/10-SCENIC/regulons_extended.csv")[str_detect(tf,"e$")]
regulons_ext[,tf:=str_remove(tf,"e$")]
res_peak_genes<-fread("outputs/27-multiomics_cbl_stim/res_peak_genes_all_genes_regulons_linkage.csv")
motifs_peaks<-GetMotif(object = cbls,peaks = unique(res_peak_genes$peak),assay = "peaks")
res_motifpeak_genes<-merge(res_peak_genes,motifs_peaks,by="peak",allow.cartesian=TRUE)
res_motifpeak_genes[,motif.name:=as.character(motif.name)]

regpeak<-unique(res_motifpeak_genes[,.(motif.name,peak,gene)])
regpeak[,tf:=motif.name]
regn<-merge(regulons_ext,regpeak,by=c("tf","gene"))
unique(regn,by=c("tf","gene")) #8694 tf-gene link

#add alterations [p_val<0.05]
  #degs
res_hsc[p_val<0.05,deg:=T]
res_hsc[p_val<0.05,expr_log2FC:=avg_log2FC]

regn<-merge(regn,res_hsc[deg==T,.(gene,deg,expr_log2FC)],all.x=T,by="gene")
  #dap
dapg[p_val<0.05,da:=T]
dapg[p_val<0.05,access_log2FC:=avg_log2FC]
regn<-merge(regn,dapg[da==T,.(peak,da,access_log2FC)],all.x=T,by="peak")

  #dmcs
cpgs_in_cres<-fread("outputs/27-multiomics_cbl_stim/cpgs_in_CREs.csv.gz")

res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

cres_meth<-merge(cpgs_in_cres,res_meth,by="cpg_id")
cres_meth[P.Value<0.05,dmc:=T]
cres_methp<-unique(cres_meth[dmc==T][order(peak,P.Value)],by=c("peak"))
cres_methp[,meth_logFC:=logFC]
regn<-merge(regn,unique(cres_methp[,.(peak,dmc,meth_logFC)]),all.x=T,by="peak")


  #if altered tf gene > all targets altered
regn[tf%in%gene[deg|dmc|da],dtf:=T]
unique(regn[dtf==T]$tf)
#  [1] "FOXP1"   "GABPA"   "ETS1"    "ELF1"    "FOS"     "JUNB"   
#  [7] "ZEB1"    "CTCF"    "IRF2"    "MEIS1"   "USF2"    "RREB1"  
# [13] "MAX"     "STAT3"   "BHLHE40" "KLF2"    "REST"    "RUNX3"  
# [19] "KLF13"   "NFKB2"   "EGR1"    "RELB"    "MEF2D"   "MAFG"   
# [25] "GATA2"   "ETV6"    "KLF10"   "GATA3"   "HLF"     "ERG"    
# [31] "ELK3"    "CEBPG"   "JUN"     "KLF6"    "MEF2C"   "NFIA"   
# [37] "THRB"    "MYB"     "HOXA5"   "NFKB1"   "TGIF1"   "TCF3"   
# [43] "MAFF"    "ETS2"    "MGA"     "MEF2A"   "ETV3"    "POU2F2" 
unique(regn$tf) #48/101 TF altered
#altered gene network

regn_alt<-regn[deg|dmc|da|dtf] #8k/12k
#n of alt by tf-gene
regn_alt[,nalt:=sum(c(deg,da,dmc,dtf),na.rm = T),by=c("tf","gene")]

regn_alt[,nalt_nodmc:=sum(c(deg,da,dtf),na.rm = T),by=c("tf","gene")]

regn_alt[nalt>2&nalt_nodmc>2][order(-nalt)][1:50]
#by tf
regn_alt[,nalt.tf:=sum(nalt),by=c("tf")]
regn_alt[,nalt.tf.norm:=nalt.tf/.N,by=c("tf")]
unique(regn_alt[order(-nalt.tf.norm)],by="tf")
#                          peak     gene      tf motif.name  deg nalt.tf nalt.tf.norm
#   1:  chr1-112465750-112466581     RHOC    KLF2       KLF2   NA 644     6.075472
#   2:      chr1-1058769-1059856     HES4    ZEB1       ZEB1   NA 1771     4.167059
#   3:   chr15-43977093-43977410  TP53BP1     MGA        MGA   NA  16     4.000000
#   4:  chr1-149886343-149888681  PLEKHO1    EGR1       EGR1   NA 1004     3.984127
#   5:    chr1-27111784-27112625    WASF2   MEF2C      MEF2C   NA
#   6:  chr1-181088263-181090448     IER5   KLF10      KLF10   NA
#   7:  chr1-108533621-108534509  FAM102B   RREB1      RREB1   NA
#   8: chr10-102393701-102395972    NFKB2   NFKB1      NFKB1   NA
#   9:  chr1-198633460-198634114    PTPRC   GATA3      GATA3   NA
plot(density(unique(regn_alt,by="tf")$nalt.tf.norm))
#linked alteration if tf neighbor > 3 norm alt.
regn_alt[,tf.alt:=nalt.tf.norm>3]
tfalts<-unique(regn_alt[(tf.alt)]$tf)
tfalts
#  [1] "GABPA" "ELF1"  "ZEB1"  "CTCF"  "RREB1" "KLF2"  "REST"  "NFKB2" "EGR1"  "MEF2D"
# [11] "GATA2" "ETV6"  "KLF10" "GATA3" "MEF2C" "NFKB1" "TCF3"  "MGA" 
regn_alt[,neighb.tf.alt:=gene%in%tfalts]
regn_alt[,linked.tfs.alt:=tf.alt&neighb.tf.alt]
unique(regn_alt[(linked.tfs.alt)],by=c("tf","gene"))
#                         peak  gene    tf motif.name  deg expr_log2FC da
#  1:  chr1-156105983-156106797 MEF2D GABPA      GABPA   NA          NA NA
#  2:  chr1-156398309-156398634 MEF2D MEF2D      MEF2D   NA          NA NA
#  3:  chr1-156786348-156787067 MEF2D  REST       REST   NA          NA NA
#  4:  chr1-156786348-156787067 MEF2D  ZEB1       ZEB1   NA          NA NA
#  5:  chr1-156859645-156861094 MEF2D  CTCF       CTCF   NA          NA NA
#  6: chr10-102393701-102395972 NFKB2  CTCF       CTCF   NA          NA NA
#  7: chr10-102393701-102395972 NFKB2  EGR1       EGR1   NA          NA NA
#  8: chr10-102393701-102395972 NFKB2  KLF2       KLF2   NA          NA NA
#  9: chr10-102393701-102395972 NFKB2 NFKB1      NFKB1   NA          NA NA
# 10: chr10-102393701-102395972 NFKB2 NFKB2      NFKB2   NA          NA NA
# 11: chr10-102616188-102616443 NFKB2  REST       REST   NA          NA NA
# 12:   chr10-31586554-31587437  ZEB1  EGR1       EGR1 TRUE  -0.4285373 NA
# 13:     chr10-8048240-8048996 GATA3  EGR1       EGR1   NA          NA NA
# 14:     chr10-8048240-8048996 GATA3  ELF1       ELF1   NA          NA NA
# 15:     chr10-8048240-8048996 GATA3 GABPA      GABPA   NA          NA NA
# 16:     chr10-8048240-8048996 GATA3 NFKB2      NFKB2   NA          NA NA
# 17:   chr15-41620970-41622187   MGA RREB1      RREB1   NA          NA NA
# 18:   chr16-67226369-67227493  CTCF GABPA      GABPA   NA          NA NA
# 19:   chr16-67226369-67227493  CTCF  REST       REST   NA          NA NA
# 20:   chr16-67561658-67563665  CTCF  CTCF       CTCF   NA          NA NA
# 21:     chr19-1248117-1249915  TCF3  CTCF       CTCF   NA          NA NA
# 22:     chr19-1248117-1249915  TCF3 GABPA      GABPA   NA          NA NA
# 23:     chr19-1248117-1249915  TCF3  TCF3       TCF3   NA          NA NA
# 24:   chr19-16330414-16330745  KLF2 GABPA      GABPA   NA          NA NA
# 25:     chr19-1811367-1813285  TCF3  REST       REST   NA          NA NA
# 26:   chr21-25606969-25608605 GABPA  CTCF       CTCF   NA          NA NA
# 27:  chr3-128456923-128457438 GATA2 GATA2      GATA2   NA          NA NA
# 28:  chr3-128467072-128467484 GATA2  ELF1       ELF1   NA          NA NA
# 29:  chr3-128492759-128493410 GATA2  EGR1       EGR1   NA          NA NA
# 30:  chr4-102500393-102502452 NFKB1  EGR1       EGR1   NA          NA NA
# 31:  chr4-102500393-102502452 NFKB1  ELF1       ELF1   NA          NA NA
# 32:  chr4-102500393-102502452 NFKB1 GABPA      GABPA   NA          NA NA
# 33:  chr4-102500393-102502452 NFKB1  ZEB1       ZEB1   NA          NA NA
# 34:  chr4-102627735-102628166 NFKB1  CTCF       CTCF   NA          NA NA
# 35:  chr4-102826505-102828592 NFKB1  REST       REST   NA          NA NA
# 36:    chr4-56653359-56653854  REST GABPA      GABPA   NA          NA NA
# 37:    chr4-56976457-56978161  REST  REST       REST   NA          NA NA
# 38:    chr4-56976457-56978161  REST  ZEB1       ZEB1   NA          NA NA
# 39:  chr5-138148859-138149080  EGR1  ELF1       ELF1   NA          NA NA
# 40:  chr5-138148859-138149080  EGR1  KLF2       KLF2   NA          NA NA
# 41:    chr5-88826833-88827723 MEF2C  ELF1       ELF1   NA          NA NA
# 42:    chr5-88826833-88827723 MEF2C GABPA      GABPA   NA          NA NA
# 43:    chr5-88883031-88884804 MEF2C MEF2C      MEF2C   NA          NA NA
# 44:      chr6-6723717-6725088 RREB1  CTCF       CTCF   NA          NA NA
# 45:      chr6-6723717-6725088 RREB1 GABPA      GABPA   NA          NA NA
# 46:      chr6-6942588-6943650 RREB1 RREB1      RREB1   NA          NA NA
# 47:  chr8-102655285-102656514 KLF10  EGR1       EGR1   NA          NA NA
# 48:  chr8-102655285-102656514 KLF10  KLF2       KLF2   NA          NA NA

#network altered tfs asso
library(ggnet)
library(network)
library(sna)

net<-as.network(unique(regn_alt[(linked.tfs.alt)][,.(tf,gene)]),loops = T,directed = T)
ggnet2(net,
       size = "degree",
       label = T,
       arrow.size = 5,arrow.gap =0.02)


#add nalt tf norm
net %v% "nalt.tf.norm" = unique(regn_alt[tf%in%network.vertex.names(net)][,.(tf,gene,nalt.tf.norm)],by="tf")[network.vertex.names(net),on="tf"]$nalt.tf.norm
nalt.tf.norm<-unique(regn_alt[tf%in%network.vertex.names(net)][,.(tf,gene,nalt.tf.norm)],by="tf")[network.vertex.names(net),on="tf"]$nalt.tf.norm
ggnet2(net,
       size = "nalt.tf.norm",
       color = "nalt.tf.norm",color.legend ="nalt.tf.norm", 
       label = T,
       arrow.size = 5,arrow.gap =0.02)



 ggnet2(net, 
              size = 0,
              #color = RImp,
              label = F,
        arrow.size = 5,arrow.gap =0.02) +
  theme(legend.text = element_text(size = 10)) +
  geom_point(aes(fill = nalt.tf.norm), size = nalt.tf.norm, col = "grey", shape = 21) +
   geom_text(aes(label = network.vertex.names(net)),nudge_y=0.04 ) +
  scale_fill_continuous(name = "norm. alterations",
                        limits=c(0, 7), breaks=seq(0, 7, by= 7),
                        low = "floralwhite" ,high = "firebrick1")
 
 
 net2<-as.network(unique(regn_alt[(tf.alt)][,.(tf,gene)]),loops = T,directed = T)
ggnet2(net2,
       size = "degree",
       label = T,
       arrow.size = 5,arrow.gap =0.02)


#add methylation infos on node, on edge 

#same for chromatin change 

#same for DEGs

#bp assos

#beginer (upstream alteration)
#methylation

#targets with at 

#unstim + stim to increase power ?####

DefaultAssay(cbls)<-"peaks"
da_peaks <- FindMarkers(
  object = cbls,
  group.by = "group",
  ident.1 = "LGA", 
  ident.2 ="CTRL",
  subset.ident = "HSC",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = c('atac_peak_region_fragments')
)

da_peaks<-data.table(da_peaks,keep.rownames = "peak")
da_peaks[p_val_adj<0.05]
da_peaks[p_val<0.001] #not better than stim alone

#which genes ?
da_peaks_genes<-merge(da_peaks,res_peak_genes,all.x=T)
da_peaks_genes[p_val<0.05&!is.na(start)]


#effect of LGA on DEGs

