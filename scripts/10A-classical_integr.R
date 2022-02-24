
source("../methyl/scripts/utils/new_utils.R")
library(Seurat)
options(future.globals.maxSize = 50000 * 1024^2)
sample_name<-"cbps0-8_clean"
out<-"outputs/10A-classical_integr/"
dir.create(out)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps.rds")
DefaultAssay(cbps)<-"RNA"
cbps_list<-SplitObject(cbps,split.by = "orig.ident")

cbps_list<-lapply(cbps_list, SCTransform,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F,
                  method = "glmGamPoi")


features <- SelectIntegrationFeatures(object.list = cbps_list, nfeatures = 3000,assay = rep("SCT",length(cbps_list)))
length(features)
cbps_list <- PrepSCTIntegration(object.list = cbps_list, anchor.features = features)
cbps_list <- lapply(X = cbps_list, FUN = RunPCA, features = features)

cbps.anchors <- FindIntegrationAnchors(object.list = cbps_list, normalization.method = "SCT", 
    anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20,
    reference = which(names(cbps_list)%in%c("cbp0_ctrl","cbp0_lga","cbp2","cbp3","cbp4")))



cbps <- IntegrateData(anchorset = cbps.anchors, normalization.method = "SCT", dims = 1:50)

cbps <- RunPCA(cbps,reduction.name = 'integrated.pca')

DimHeatmap(cbps,dims = 1:6,cells=500,balanced = T,reduction = 'integrated.pca')

cbps <- RunUMAP(cbps, reduction = "integrated.pca",reduction.name = "integrated.umap", dims = c(1:50))

rm(cbps.anchors)

DimPlot(cbps,group.by = "lineage_hmap",reduction = "integrated.umap",label=T)

DimPlot(cbps,group.by = "cell_type_hmap",reduction = "integrated.umap",label=T)

saveRDS(cbps,fp(out,paste0(sample_name,".rds")))


message("Success !")
