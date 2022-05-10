out<-"outputs/05-make_hematomap/"
library(Seurat)
source("scripts/utils/new_utils.R")

batch<-"CBP7c"
matrix_path<-"~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-c/outs/filtered_feature_bc_matrix/"

#CBP7b contain CBP537 (CTRL M)
sample_male<-"ctrlM537"


cbp<-Read10X(matrix_path)

cbp<-CreateSeuratObject(cbp,project = batch)
cbp<-PercentageFeatureSet(cbp,pattern = "MT-",col.name = "percent.mt")
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt")) 
cbp<-subset(cbp,nCount_RNA<35000&nFeature_RNA<6000&nFeature_RNA>200&percent.mt<15)

cbp #33538 features across 6594  samples
cbp[["sample"]]<-sample_male

saveRDS(subset(cbp,sample==sample_male),fp(out,ps(sample_male,".rds")))

