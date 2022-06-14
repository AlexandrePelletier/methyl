
out<-"outputs/21C-Improved_CREs_inference/"

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

cbl12<-readRDS("outputs/21-CREs_inference_with_sc_multiomics_cbls/cbl12.rds")
DefaultAssay(cbl12)<-"peaks"

cbl12 <- LinkPeaks(
  object = cbl12,
  peak.assay = "peaks",
    score_cutoff = 0,
  pvalue_cutoff = 0,
  expression.assay = "SCT")

saveRDS(Links(cbl12),fp(out,"cbl12_multiomics_peaks_genes_links.rds"))

peak_links_dt<-data.table(as.data.frame(Links(cbl12)))
fwrite(peak_links_dt,fp(out,"res_cbl12_multiomics_peaks_genes_links.csv.gz"))

message("Success !")
