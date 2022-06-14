out<-"outputs/21-CREs_inference_with_sc_multiomics_cbls"
dir.create(out)

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

# renv::install("bioc::biovizBase")

#renv::install('Bioconductor/BiocGenerics')
#renv::install('Bioconductor/GenomeInfoDb')
#renv::install("bioc::EnsDb.Hsapiens.v86")
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


peak_links_dt<-data.table(as.data.frame(Links(cbl12)))
saveRDS(cbl12,fp(out,"cbl12.rds"))
fwrite(peak_links_dt,fp(out,"res_peaks_genes_links_important_genes.csv"))

#all peak genes links
#run 21A
#n of peak gene link, what are the limits, how improve with MULTI GRN ++

peak_links_dt<-fread(fp(out,"res_peaks_genes_links_all_genes.csv.gz"))

peak_links_dt #4590 peak gene link
unique(peak_links_dt,by="gene") #2664 genes

unique(peak_links_dt,by="peak") #4170 peaks

#n links in EGR1N
egr1n<-fread("outputs/16-GRN_final/egr1_KLF2_KLF4_network_tf_target_interactions.csv")
egr1n_genes<-unique(c(egr1n$target,egr1n$tf))
egr1n_genes# 123 genes

unique(peak_links_dt[gene%in%egr1n_genes],by="gene") #41/123 egr1n genes have linked peak

unique(peak_links_dt[gene%in%egr1n_genes]$gene)
#  [1] "ID3"     "TINAGL1" "LMNA"    "ATF3"    "ID2"     "ZFP36L2" "PELI1"   "CXCR4"   "ARL4C"  
# [10] "CTNNB1"  "NFKBIZ"  "TIPARP"  "SKIL"    "HES1"    "EGR1"    "IER3"    "BRD2"    "CDKN1A" 
# [19] "EZR"     "NDRG1"   "OTUD1"   "DDIT4"   "NFKB2"   "CD151"   "GRASP"   "BTG1"    "PNP"    
# [28] "FOS"     "CALM1"   "EIF5"    "AKAP13"  "SOCS3"   "NFATC1"  "MIDN"    "GADD45B" "IER2"   
# [37] "ZFP36"   "FOSB"    "ZNF667"  "PRNP"    "MAFF" 

#n of peaks by genes
table(peak_links_dt[gene%in%egr1n_genes]$gene)
#  AKAP13   ARL4C    ATF3    BRD2    BTG1   CALM1   CD151  CDKN1A  CTNNB1   CXCR4   DDIT4    EGR1 
#       1       2       2       1       1       1       1       2       1       2       1       5 
#    EIF5     EZR     FOS    FOSB GADD45B   GRASP    HES1     ID2     ID3    IER2    IER3    LMNA 
#       1       1       4       3       5       1       4       5       1       2       5       1 
#    MAFF    MIDN   NDRG1  NFATC1   NFKB2  NFKBIZ   OTUD1   PELI1     PNP    PRNP    SKIL   SOCS3 
#       2       2       3       1       1       2       1       2       1       1       1       2 
# TINAGL1  TIPARP   ZFP36 ZFP36L2  ZNF667 
#       1       6       4       1       1 

#merge with DMCs
peak_links_dt[,chr:=seqnames]
peak_links_dt[,start.peak:=start(peak)]
peak_links_dt[,end.peak:=end(peak),by="peak"]

res_cpgs<-fread("outputs/14-DMCs_atac_integr/res_cpgs_hg38.cs.gz")
cpgs_in_cres<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          peak_links_dt[,.(chr,start.peak,end.peak,peak)][order(chr,start.peak,end.peak)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peak","start","end"))

length(unique(cpgs_in_cres$cpg_id))#7854/750k cpgs
length(unique(cpgs_in_cres$peak)) #1.7k/4.2k peaks

fwrite(cpgs_in_cres,fp(out,"cpgs_in_CREs.csv.gz"))

#% DMCs overlapping CREs
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

cres_meth<-merge(cpgs_in_cres,res_meth,by="cpg_id")
nrow(cres_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]) #2.5% DMCs are in CREs 

#vs % CpGs overlapping 
nrow(cres_meth)/nrow(res_meth) #1% CpGs are in CREs

(nrow(cres_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]))/(nrow(cres_meth)/nrow(res_meth)) #2.3 foldenrichment


OR(set1 = res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id,
   set2= cres_meth$cpg_id,size_universe = nrow(res_meth)) #p = 3.627709e-15

#in EGRn
#DMCs overlapping CREs EGRN
cres_cpgs_egrn<-merge(cpgs_in_cres,peak_links_dt[gene%in%egr1n_genes][,-c("start","end")])
cres_meth_egrn<-merge(cres_cpgs_egrn,res_meth,by="cpg_id")
cres_meth_egrn[P.Value<0.001&abs(logFC)>25] #1DMCs only, in CRE of NFATC1

#vs CpGs overlapping 
nrow(cres_meth_egrn) #472 CpGs are in CREs EGRN

#increase Meth in CREs EGRN vs CREs ?
cres_meth[,in_EGRN:=cpg_id%in%cres_meth_egrn$cpg_id]

ggplot(cres_meth)+geom_boxplot(aes(x=in_EGRN,y=-log10(P.Value)*logFC)) #nop


 #CCL : only 5k identified CREs, (2.5% DMCs and only 1% CpG are in it),
#41/123 Egr1 genes have CREs, only 472 CpGs and 1DMCs fall in it.
##==> dont capture stress responsive methylation impacted CREs.

#NEED 1) do scMulti on +/- Stim 
# 2) improve CREs indentif, need know which celltype this CREs are active, and build a GRN based on it.

#1) 
#can we use HTO to multiplex conditons ???
#=> cbl12 HTO
#see 21B


#Annexe####
#add TF motif on cbl12
cbl12<-readRDS("outputs/21-CREs_inference_with_sc_multiomics_cbls/cbl12.rds")
library(JASPAR2020)
library(TFBSTools)

#renv::install("bioc::BSgenome.Hsapiens.UCSC.hg38")
#renv::install("bioc::motifmatchr")

library(BSgenome.Hsapiens.UCSC.hg38)

# Get a list of motif position frequency matrices from the JASPAR database
#?getMatrixSet
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
#?AddMotifs
DefaultAssay(cbl12)<-"peaks"
cbl12 <- AddMotifs(
  object = cbl12,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

saveRDS(cbl12,fp(out,"cbl12.rds"))

