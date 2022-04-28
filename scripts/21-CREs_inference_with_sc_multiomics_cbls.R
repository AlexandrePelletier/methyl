out<-"outputs/21-CREs_inference_with_sc_multiomics_cbls"
dir.create(out)

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

# renv::install("bioc::biovizBase")

#renv::install('Bioconductor/BiocGenerics')
#renv::install('Bioconductor/GenomeInfoDb')
renv::install("bioc::EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)


# load the RNA and ATAC data
counts <- Read10X_h5("~/RUN/Run_699_single_cell/Output/cellranger-arc/multi2-cbl12-atac/outs/filtered_feature_bc_matrix.h5")
fragpath <- "~/RUN/Run_699_single_cell/Output/cellranger-arc/multi2-cbl12-atac/outs/atac_fragments.tsv.gz"
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
#genome(annotation) <- "hg38"

head(annotation)

# create a Seurat object containing the RNA adata
cbl12 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
cbl12[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
cbl12
# An object of class Seurat 
# 154050 features across 7841 samples within 2 assays 
# Active assay: ATAC (117449 features, 0 variable features)
#  1 other assay present: RNA

#QC
DefaultAssay(cbl12) <- "ATAC"

cbl12 <- NucleosomeSignal(cbl12)
cbl12 <- TSSEnrichment(cbl12)

VlnPlot(
  object = cbl12,
  group.by = "orig.ident",
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
cbl12 <- subset(
  x = cbl12,
  subset = nCount_ATAC < 50000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1.25 &
    TSS.enrichment > 1
)

cbl12 # 6502 cells

#peak calling 
# call peaks using MACS2
#reticulate::py_install(packages ="MACS2")

peaks <- CallPeaks(cbl12, macs2.path = "renv/python/virtualenvs/renv-python-3.7/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(cbl12),
  features = peaks,
  cells = colnames(cbl12)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
cbl12[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
saveRDS(cbl12,fp(out,"cbl12.rds"))

#Gene expression data processing
DefaultAssay(cbl12) <- "RNA"
cbl12 <- SCTransform(cbl12)
cbl12 <- RunPCA(cbl12)
cbl12 <- CellCycleScoring(cbl12,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)
cbl12$CC.Difference <- cbl12$S.Score - cbl12$G2M.Score
cbl12 <- RunUMAP(cbl12,dims = 1:30)
DimPlot(cbl12,group.by="Phase")
cbl12[["percent.mt"]] <- PercentageFeatureSet(object = cbl12, assay = "RNA",pattern = "^MT-")
FeaturePlot(cbl12,"percent.mt")
cbl12 <- SCTransform(cbl12,vars.to.regress = c("percent.mt","CC.Difference"))

#DNA accessibility data processing
DefaultAssay(cbl12) <- "peaks"
cbl12 <- FindTopFeatures(cbl12, min.cutoff = 5)
cbl12 <- RunTFIDF(cbl12)
cbl12 <- RunSVD(cbl12)

# load CBPs reference (hematomap)
hmap <- readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
DefaultAssay(hmap) <- "integrated"

DefaultAssay(cbl12) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = hmap,
  query = cbl12,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = hmap$cell_type,
  weight.reduction = cbl12[['pca']],
  dims = 1:50
)

cbl12 <- AddMetaData(
  object = cbl12,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(cbl12) <- "predicted.id"

# # set a reasonable order for cell types to be displayed when plotting
# levels(cbl12) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
#                   "CD8 Naive", "dnT",
#                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
#                  "NK Proliferating", "gdT",
#                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
#                  "CD14 Mono", "CD16 Mono",
#                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

#Joint UMAP visualization
#Using the weighted nearest neighbor methods in Seurat v4, we can compute a joint neighbor graph that represent both the gene expression and DNA accessibility measurements.

# build a joint neighbor graph using both assays
cbl12 <- FindMultiModalNeighbors(
  object = cbl12,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
cbl12 <- RunUMAP(
  object = cbl12,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(cbl12, label = TRUE, reduction = "umap") + NoLegend()

cbl12$cell_type_hmap<-Idents(cbl12)
cbl12[["lineage_hmap"]]<-sapply(as.character(cbl12@meta.data$cell_type_hmap), function(ct){
  if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
  else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
  else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
  else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
  else if(ct%in%c("CLP","proB"))return("Lymphoid")
  else return(ct)
  
})
DimPlot(cbl12,group.by = "lineage_hmap", label = TRUE, reduction = "umap") + NoLegend()

FeaturePlot(cbl12,c("AVP","MPO","CD99","VPREB1","GATA1"),reduction = "umap",max.cutoff = "q95")

saveRDS(cbl12,fp(out,"cbl12.rds"))
cbl12<-readRDS(fp(out,"cbl12.rds"))

# #map on hematomap
# cbl12<-SCTransform(cbl12,vars.to.regress=c("percent.mt","CC.Difference"),
#                   return.only.var.genes=F, 
#                   method = "glmGamPoi")
# hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")
# 
# anchors <- FindTransferAnchors(
#     reference = hmap,
#     query = cbl12,
#     k.filter = NA,
#     reference.reduction = "pca", 
#     reference.neighbors = "pca.annoy.neighbors", 
#     dims = 1:50
#   )
# 
# cbl12 <- MapQuery(
#     anchorset = anchors, 
#     query = cbl12,
#     reference = hmap, 
#     refdata = list(
#       cell_type = "cell_type", 
#       lineage = "lineage"),
#     reference.reduction = "pca",
#     reduction.model = "ref.umap"
#   )
# 
# DimPlot(cbl12, reduction = "ref.umap",
#         group.by =  "predicted.cell_type",
#         label = TRUE, repel = TRUE, label.size = 3) +
#   NoLegend()
# #doesnt work well because all cells are MPP and cluster in the center of the hmap
# 
# 
# saveRDS(cbl12,fp(out,"cbl12.rds"))


#Linking peaks to genes
#correlation between gene expression and accessibility at nearby peaks, and correcting for bias due to GC content, overall accessibility, and peak size
# first compute the GC content for each peak
DefaultAssay(cbl12)<-"peaks"
Idents(cbl12)<-"lineage_hmap"

cbl12 <- RegionStats(cbl12, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
cbl12 <- LinkPeaks(
  object = cbl12,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("EGR1", "KLF2","SOCS3","KLF13","MPO","GATA1","CD99","VPREB1")
)
#We can visualize these links using the CoveragePlot() function, or alternatively we could use the CoverageBrowser() function in an interactive analysis:

idents.plot <- c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

p1 <- CoveragePlot(
  object = cbl12,
  region = "SOCS3",
  features = "SOCS3",
  expression.assay = "SCT",
 idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = cbl12,
  region = "GATA1",
  features = "GATA1",
  expression.assay = "SCT",
 idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)


p1 <- CoveragePlot(
  object = cbl12,
  region = "EGR1",
  features = "EGR1",
  expression.assay = "SCT",
 idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = cbl12,
  region = "KLF2",
  features = "KLF2",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)


p1 <- CoveragePlot(
  object = cbl12,
  region = "VPREB1",
  features = "VPREB1",
  expression.assay = "SCT",
 idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = cbl12,
  region = "CD99",
  features = "CD99",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)


p1 <- CoveragePlot(
  object = cbl12,
  region = "KLF13",
  features = "KLF13",
  expression.assay = "SCT",
 idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = cbl12,
  region = "MPO",
  features = "MPO",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)

saveRDS(cbl12,fp(out,"cbl12.rds"))
