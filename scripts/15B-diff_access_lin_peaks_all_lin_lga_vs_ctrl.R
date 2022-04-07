out<-'outputs/15-chromatin_change_LGA_vs_Ctrl/'

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
DefaultAssay(atacs)<-"lin_peaks"
Idents(atacs)<-"predicted.id"
atacs$group<-ifelse(str_detect(atacs$dataset,"1|3"),"ctrl","lga")
lins<-c("HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")
peaks_lga<-Reduce(rbind,lapply(lins,function(lin)data.table(FindMarkers(atacs,
                           subset.ident = lin,
                          group.by = "group",
                          ident.1 = 'lga',
                          ident.2 = 'ctrl',
                          only.pos = F,
                          min.pct = 0.05,
                          logfc.threshold = 0,
                          test.use = "LR",
                          latent.vars = "peak_region_fragments"),
                          keep.rownames = "peak")[,lineage:=lin]))

fwrite(peaks_lga,fp(out,"differential_peaks_accessibility_lga_vs_ctrl_by_lin_logFC0.csv.gz"))

message("Success !")
