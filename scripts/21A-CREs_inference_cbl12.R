
out<-"outputs/21-CREs_inference_with_sc_multiomics_cbls"

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

cbl12<-readRDS(fp(out,"cbl12.rds"))
DefaultAssay(cbl12)<-"peaks"
Idents(cbl12)<-"lineage_hmap"

cbl12 <- LinkPeaks(
  object = cbl12,
  peak.assay = "peaks",
  expression.assay = "SCT")

saveRDS(cbl12,fp(out,"cbl12.rds"))

peak_links_dt<-data.table(as.data.frame(Links(cbl12)))
fwrite(peak_links_dt,fp(out,"res_peaks_genes_links_all_genes.csv.gz"))

