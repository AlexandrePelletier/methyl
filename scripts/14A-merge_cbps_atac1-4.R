
#merge datasets and peaks (ref = https://satijalab.org/signac/articles/merging.html)

out<-"outputs/14-DMCs_atac_integr/"


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("scripts/utils/new_utils.R")
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
# read in peak sets
peaks.atac1 <- read.table(
  file = "../atac/datasets/run_542_atac1/peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.atac2 <- read.table(
  file = "../atac/datasets/run_542_atac2/peaks.bed",
  col.names = c("chr", "start", "end")
  )
peaks.atac3 <- read.table(
  file = "../atac/datasets/run_542_atac3/peaks.bed",
  col.names = c("chr", "start", "end")
  )

peaks.atac4 <- read.table(
  file = "../atac/datasets/run_542_atac4/peaks.bed",
  col.names = c("chr", "start", "end")
  )

# convert to genomic ranges
gr.atac1 <- makeGRangesFromDataFrame(peaks.atac1)
gr.atac2 <- makeGRangesFromDataFrame(peaks.atac2)
gr.atac3 <- makeGRangesFromDataFrame(peaks.atac3)
gr.atac4 <- makeGRangesFromDataFrame(peaks.atac4)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.atac1, gr.atac2, gr.atac3, gr.atac4))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.atac1 <- read.table(
  file = "../atac/datasets/run_542_atac1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.atac2 <- read.table(
  file = "../atac/datasets/run_542_atac2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac3 <- read.table(
  file = "../atac/datasets/run_542_atac3/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac4 <- read.table(
  file = "../atac/datasets/run_542_atac4/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.atac1 <- md.atac1[md.atac1$passed_filters > 500, ]
md.atac2 <- md.atac2[md.atac2$passed_filters > 500, ]
md.atac3 <- md.atac3[md.atac3$passed_filters > 500, ]
md.atac4 <- md.atac4[md.atac4$passed_filters > 500, ] 

# create fragment objects
frags.atac1 <- CreateFragmentObject(
  path = "../atac/datasets/run_542_atac1/fragments.tsv.gz",
  cells = rownames(md.atac1)
)

frags.atac2 <- CreateFragmentObject(
  path = "../atac/datasets/run_542_atac2/fragments.tsv.gz",
  cells = rownames(md.atac2)
)

frags.atac3 <- CreateFragmentObject(
  path = "../atac/datasets/run_542_atac3/fragments.tsv.gz",
  cells = rownames(md.atac3)
)

frags.atac4 <- CreateFragmentObject(
  path = "../atac/datasets/run_542_atac4/fragments.tsv.gz",
  cells = rownames(md.atac4)
)

cbp.atac1.counts <- FeatureMatrix(
  fragments = frags.atac1,
  features = combined.peaks,
  cells = rownames(md.atac1)
)

cbp.atac2.counts <- FeatureMatrix(
  fragments = frags.atac2,
  features = combined.peaks,
  cells = rownames(md.atac2)
)

cbp.atac3.counts <- FeatureMatrix(
  fragments = frags.atac3,
  features = combined.peaks,
  cells = rownames(md.atac3)
)

cbp.atac4.counts <- FeatureMatrix(
  fragments = frags.atac4,
  features = combined.peaks,
  cells = rownames(md.atac4)
)

cbp.atac1_assay <- CreateChromatinAssay(cbp.atac1.counts, fragments = frags.atac1)
cbp.atac1 <- CreateSeuratObject(cbp.atac1_assay, assay = "ATAC")

cbp.atac2_assay <- CreateChromatinAssay(cbp.atac2.counts, fragments = frags.atac2)
cbp.atac2 <- CreateSeuratObject(cbp.atac2_assay, assay = "ATAC")

cbp.atac3_assay <- CreateChromatinAssay(cbp.atac3.counts, fragments = frags.atac3)
cbp.atac3 <- CreateSeuratObject(cbp.atac3_assay, assay = "ATAC")

cbp.atac4_assay <- CreateChromatinAssay(cbp.atac4.counts, fragments = frags.atac4)
cbp.atac4 <- CreateSeuratObject(cbp.atac4_assay, assay = "ATAC")


# add information to identify dataset of origin
cbp.atac1$dataset <- 'cbp.atac1'
cbp.atac2$dataset <- 'cbp.atac2'
cbp.atac3$dataset <- 'cbp.atac3'
cbp.atac4$dataset <- 'cbp.atac4'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = cbp.atac1,
  y = list(cbp.atac2, cbp.atac3, cbp.atac4),
  add.cell.ids = c("atac1", "atac2", "atac3", "atac4")
)
combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

saveRDS(combined,fp(out,"cbps_atac1-4_merged.rds"))

DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
ggsave(fp(out,"umap_combined.png"))

CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)

ggsave(fp(out,"region_check_coverage.png"))


#QC filtering
source("../methyl/scripts/utils/new_utils.R")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

cbps<-readRDS(fp(out,"cbps_atac1-4_merged.rds"))
cbps
# 126925 features across 232072 samples within 1 assay 
# Active assay: ATAC (126925 features, 126920 variable features)
#  2 dimensional reductions calculated: lsi, umap
head(cbps@meta.data)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)


# change to UCSC style because the data was mapped to hg38/UCSC
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

saveRDS(annotations,"../atac/ref/gene_annotations_hg38_GRanges.rds")

# add the gene information to the object
Annotation(cbps) <- annotations

#add metadata to the object

md.atac1 <- read.table(
  file = "../atac/datasets/run_542_atac1/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.atac2 <- read.table(
  file = "../atac/datasets/run_542_atac2/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac3 <- read.table(
  file = "../atac/datasets/run_542_atac3/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.atac4 <- read.table(
  file = "../atac/datasets/run_542_atac4/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

head(md.atac1)
head(cbps@meta.data)
rownames(md.atac1)<-paste0("atac1_",rownames(md.atac1))
rownames(md.atac2)<-paste0("atac2_",rownames(md.atac2))
rownames(md.atac3)<-paste0("atac3_",rownames(md.atac3))
rownames(md.atac4)<-paste0("atac4_",rownames(md.atac4))

head(md.atac1)
md.cbps<-Reduce(rbind,list(md.atac1,md.atac2,md.atac3,md.atac4))
head(md.cbps)

cbps <- AddMetaData(object = cbps, metadata = md.cbps)

#QC Metrics

# compute nucleosome signal score per cell
#approximate ratio of mononucleosomal to nucleosome-free fragments (stored as nucleosome_signal)
cbps <- NucleosomeSignal(object = cbps)

# compute TSS enrichment score per cell : ratio of fragments centered at the TSS to fragments in TSS-flanking regions
cbps <- TSSEnrichment(object = cbps, fast = FALSE)



#add fraction of reads in peaks (cellular sequencing depth / complexity)
cbps$pct_reads_in_peaks <- cbps$peak_region_fragments / cbps$passed_filters * 100
#and blacklist ratio (reads which are often associated with artefactual signal
cbps$blacklist_ratio <- cbps$blacklist_region_fragments / cbps$peak_region_fragments

#validate that TSS enrichment scores compute really represent enrichment for cells around TSS :
cbps$high.tss <- ifelse(cbps$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(cbps, group.by = 'high.tss') + NoLegend()
TSSPlot(cbps, group.by = 'dataset') + NoLegend()

#fragment length periodicity for all the cells :
#need mononucleosomal / nucleosome-free ratio < 4 to have a good Fragment length profile for ATACseq exp :
sum(cbps$nucleosome_signal > 4) #143
VlnPlot(cbps,"nucleosome_signal",group.by="dataset")
cbps$nucleosome_group <- ifelse(cbps$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = cbps, group.by = 'nucleosome_group')


VlnPlot(cbps,"peak_region_fragments",pt.size = 0.1,log=T,group.by="dataset")+geom_hline(yintercept = 5000)

cbps<-subset(cbps,peak_region_fragments>5000)
cbps #126925 features across 9699 samples
head(cbps@meta.data)

VlnPlot(
  object = cbps,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

VlnPlot(cbps,"pct_reads_in_peaks",group.by="dataset",pt.size = 0.1)+geom_hline(yintercept = 15)


VlnPlot(cbps,"TSS.enrichment",pt.size = 0.1,group.by="dataset")+geom_hline(yintercept = 2)

VlnPlot(cbps,"nucleosome_signal",pt.size = 0.1,group.by="dataset")+geom_hline(yintercept = 1)

VlnPlot(cbps,"blacklist_ratio",pt.size = 0.1,group.by="dataset")+geom_hline(yintercept = 0.0015)


table(subset(
  x = cbps,
  subset = peak_region_fragments > 5000 & 
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.0015 &
    nucleosome_signal < 1&
    TSS.enrichment > 2& TSS.enrichment <10 
)$dataset)

cbps<-subset(
  x = cbps,
  subset = peak_region_fragments > 5000 & 
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.0015 &
    nucleosome_signal < 1&
    TSS.enrichment > 2& TSS.enrichment <10 
)

cbps #126925 features across 8733 samples 

#Normalization and linear dimensional reduction
#1) term frequency-inverse document frequency (TF-IDF) normalization.
#     across cells to correct for differences in cellular sequencing depth, 
#     and across peaks to give higher values to more rare peaks
        
cbps <- RunTFIDF(cbps)
cbps <- FindTopFeatures(cbps, min.cutoff = 'q0') #use only the top n% of features (peaks) for dimensional reduction, or remove features present in less than n cells
cbps <- RunSVD(cbps) #singular value decomposition (SVD) on the TD-IDF matrix
#TF-IDF + SVD = latent semantic indexing (LSI),

# first LSI component often captures sequencing depth (technical variation) rather than biological variation
DepthCor(cbps) #++ correl comp.1 and depth => perform downstream steps without this component

#Non-linear dimension reduction and clustering
#same than for scRNA-seq
cbps <- RunUMAP(object = cbps, reduction = 'lsi', dims = 2:30)

DimPlot(object = cbps, group.by = 'dataset') #batch effect
p1<-DimPlot(object = cbps, group.by = 'dataset') +ggtitle("unintegrated")

#Integration with Harmony
library(harmony)
#error https://github.com/satijalab/seurat/issues/1849
#need adding assay=assay.use
# and @scale.data slot

head(cbps@assays$ATAC@scale.data)

cbps <- RunHarmony(
  object = cbps,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
cbps <- RunUMAP(cbps, dims = 2:30, reduction = 'harmony',reduction.name = "humap")
p2 <- DimPlot(cbps, group.by = 'dataset',reduction = "humap", pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")
p1 + p2

cbps <- FindNeighbors(object = cbps, reduction = 'harmony', dims = 2:30)
cbps <- FindClusters(object = cbps, verbose = FALSE, algorithm = 3)
DimPlot(object = cbps, reduction = "humap",label = TRUE) + NoLegend()


saveRDS(cbps,fp(out,"cbps_atac1-4_merged_qc.rds"))

