out<-"outputs/E01-scIntegration_QC"
dir.create(out)
source("scripts/utils/new_utils.R")

#assess quality of integration by evaluating biological conservation over batch correction
#using  scib package

reticulate::py_install("scib")

#anndata
obj
