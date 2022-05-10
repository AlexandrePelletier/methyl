### Project Setup ==================================================================================
library(here)
out <- here("outputs", "10-SCENIC")
dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  #renv::install("bioc::AUCell")
  library(Seurat)
  source("scripts/utils/new_utils.R")
  library(AUCell)
  library(SCENIC)
  source("scripts/utils/scenic_utils.R")
})



### Analysis =======================================================================================

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
regulons_list<-readRDS("outputs/10-SCENIC/cbps_14k/regulons_list.rds")

#with AUCell
cells_rankings <- AUCell_buildRankings(as.matrix(cbps@assays$RNA@counts),nCores = 20)

cells_AUC <- AUCell_calcAUC(regulons_list, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
auc_mat<-getAUC(cells_AUC)
cbps@assays[["TF_AUC"]] <- CreateAssayObject(auc_mat)
DefaultAssay(cbps)<-"TF_AUC"

saveRDS(cbps,fp(out,"cbps_with_regulons_activity.rds"))


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
