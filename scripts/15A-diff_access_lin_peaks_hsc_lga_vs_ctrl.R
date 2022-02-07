out<-'outputs/14-DMCs_atac_integr/'

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
DefaultAssay(atacs)<-"lin_peaks"
Idents(atacs)<-"predicted.id"
peaks_hsc_lga<-FindMarkers(atacs,
                           subset.ident = "HSC",
                          group.by = "group",
                          ident.1 = 'lga',
                          ident.2 = 'ctrl',
                          only.pos = F,
                          min.pct = 0.05,
                          logfc.threshold = 0,
                          test.use = "LR",
                          latent.vars = "peak_region_fragments")
peaks_hsc_lga<-data.table(peaks_hsc_lga,keep.rownames = "peak")

fwrite(peaks_hsc_lga,fp(out,"differential_peaks_accessibility_lga_vs_ctrl_hsc_logFC0.csv.gz"))

message("Success !")
