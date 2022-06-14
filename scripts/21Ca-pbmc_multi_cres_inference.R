out<-"outputs/21C-Improved_CREs_inference/"
dir.create(out)

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

pbmc<-readRDS(fp(out,"pbmc_multiomics.rds"))
DefaultAssay(pbmc)<-"peaks"
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  score_cutoff = 0,
  pvalue_cutoff = 0,
  expression.assay = "SCT")

saveRDS(Links(pbmc),fp(out,"pbmc_multiomics_peaks_genes_links.rds"))

peak_links_dt<-data.table(as.data.frame(Links(pbmc)))
fwrite(peak_links_dt,fp(out,"res_pbmc_multiomics_peaks_genes_links.csv.gz"))

message("Success !")
