
out<-"outputs/21C-Improved_CREs_inference"
dir.create(out)

library(Seurat)
library(Signac)
source("scripts/utils/new_utils.R")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

#get toy seurat object 
##pbmc mulltiomics
system("mkdir ../singlecell/datasets/pbmc_multiomics")
system("wget -O ../singlecell/datasets/pbmc_multiomics/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5 https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
system("wget -O ../singlecell/datasets/pbmc_multiomics/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
system("wget -O ../singlecell/datasets/pbmc_multiomics/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi")


#CREs inference with Signac
#preprocessing####
# load the RNA and ATAC data
counts <- Read10X_h5("../singlecell/datasets/pbmc_multiomics/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "../singlecell/datasets/pbmc_multiomics/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
pbmc

#QC
DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc

# call peaks using MACS2
peaks <- CallPeaks(pbmc, macs2.path = "renv/python/virtualenvs/renv-python-3.7/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

#Gene expression data processing
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

#DNA accessibility data processing
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

#Annotating cell types
library(SeuratDisk)

# load PBMC reference
reference <- LoadH5Seurat("../singlecell/datasets/pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                 "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)
DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
saveRDS(pbmc,fp(out,"pbmc_multiomics.rds"))

#Linking peaks to genes####
DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)
DefaultAssay(pbmc)<-"SCT"
#How that work ? ####
# For each gene,
# we computed the Pearson correlation coefficient r between
# the gene expression and the accessibility of each peak within
# 500 kb of the gene TSS. For each peak, we then computed
# a background set of expected correlation coefficients given
# properties of the peak by randomly sampling 200 peaks located on a different chromosome to the gene, matched for
# GC content, accessibility, and sequence length (MatchRegionStats function in Signac

# We then computed the Pearson
# correlation between the expression of the gene and the set of
# background peaks. A z-score was computed for each peak
# as z = (r − µ)/σ, where µ was the background mean correlation coefficient and σ was the standard deviation of the
# background correlation coefficients for the peak. We computed a p-value for each peak using a one-sided z-test, and
# retained peak-gene links with a p-value < 0.05 and a Pearson
# correlation coefficient > 0.05 or < -0.05

#Issue :
#1) don't take into account drop out effect. 
#=> better if smooth expression / DNA access with k Neighbours cells ?

#2) pearson correlation, so linear relationship,
#=>  better with Spearman/kendall ? other correlation methods ? (e.g random forest / GENIE method)

#3) 200 sampling, so pvalue low can't be <0.005 
#=> better with more sampling, others pvalue strategy ?

#4) Hard to have significant correlation for little cluster CREs : 
#random sampling of matching criteria OCRs (aspecially accessibility) will tend to choose cluster specific OCRs ?


#But, what's it better ??
#better is greater TRUE Positive/FALSE Negative ratio (sensitivity) while Lower FALSE Positive / TRUE Negative ratio (specificity)
#=> Greater sensitivity / specificity ratio
#not possible to estimate False Negative or True Negative, so will focus on TRUE positive rate / FALSE positive rate.
#what's a TRUE Positive ? => assume peak gene link TRUE if peak in close proximity to TSS (+/- 1 kb) 
#=> TRUE Positive rate = number of peaks in TSS region asso to gene expression / number of peaks in TSS region
#what's a FALSE positive ? => peak in an other Chromosome while associated => FALSE positive rate = pvalue threshold



#CREs inference improvement####
#0) what are the TPr, FPr, and AUC of signac method ?
#cbl12
cbl12_links<-fread("outputs/21-CREs_inference_with_sc_multiomics_cbls/res_peaks_genes_links_all_genes.csv.gz")

pvalue_threshold<-0.05
cbl12_links[,sig.link:=pvalue<pvalue_threshold&score>0.05]
cbl12_links[,n_sig.link:=sum(sig.peak)]


cbl12_links[,gene.start:=start]
cbl12_links[,peak.start:=start(peak)]
cbl12_links[,tss_dist:=peak.start-gene.start]
plot(density(cbl12_links$tss_dist))

cbl12_links[,peaks_TSS:=abs(tss_dist)<1000] 

length(unique(cbl12_links[(peaks_TSS)]$peak))#2176/4170
length(unique(cbl12_links$peak))#4170 => 52% significant peaks are < 1kb
cbl12_links[,n_peaks_TSS:=sum(peaks_TSS)]

cbl12_links[,n_TP:=pvalue<pvalue_threshold&peaks_TSS]

cbl12_links[,TP_rate:=n_TP/n_peaks_TSS]
cbl12_links[,FP_rate:=pvalue_threshold]

#pbmcs
#run 21Ca-pbmc_multi_cres_inference.R


#1) take into account drop out effect ?####

#=> cellranger pipeline do that with matrix smoothing: 
#matrix smoothing : The value of a feature in a given cell is enhanced by “borrowing” values of the same feature from “neighboring” cells. 
#=> determined by applying K-nearest neighbor (KNN) algorithm on the gene expression PCA reduced dimensions.
#The information “borrowing” is achieved by a weighted sum of signals from 30 closest neighboring cells.
#The smoothing weight is defined as the Gaussian kernel transformation of the euclidean distance between two cells
#in PCA space. The Gaussian kernel transformation guarantees that the smoothing weight is high only when two cells 
#have a highly similar gene expression profile and quickly decays to zero when the similarity between cells decreases.
#This is necessary to prevent over-smoothing in KNN-based methods (Andrews and Hemberg, F1000 Research 2019).
#a) compared with 10X results.
#=> TPr, FPr, and AUC of cellranger method 

