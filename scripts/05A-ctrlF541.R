out<-"outputs/05-make_hematomap/"
library(Seurat)
source("scripts/utils/new_utils.R")

batch<-"CBP6-c"
matrix_path<-"~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-c/outs/filtered_feature_bc_matrix/"

#CBP6c contain CBP541 (CTRL F) et 526 (LGA M)
cbp<-Read10X(matrix_path)

cbp<-CreateSeuratObject(cbp,project = batch)
cbp<-PercentageFeatureSet(cbp,pattern = "MT-",col.name = "percent.mt")
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt")) 
cbp<-subset(cbp,nCount_RNA<20000&nFeature_RNA<4000&nFeature_RNA>200&percent.mt<25)


#DEMultiplex
VlnPlot(cbp,c("XIST","RPS4Y1"),group.by = "orig.ident")
FeaturePlot(cbp,c("XIST","RPS4Y1"))

sum(cbp@assays$RNA@data["XIST",]>0) 
sum(cbp@assays$RNA@data["RPS4Y1",]>0) 
sum(cbp@assays$RNA@data["XIST",]>0&cbp@assays$RNA@data["RPS4Y1",]>0) 
sum(cbp@assays$RNA@data["XIST",]==0&cbp@assays$RNA@data["RPS4Y1",]==0) 


cbp@meta.data[cbp@assays$RNA@data["XIST",]>0&
                  cbp@assays$RNA@data["RPS4Y1",]>0,"sample"]<-"Doublet"

cbp@meta.data[cbp@assays$RNA@data["XIST",]>0&
                  cbp@assays$RNA@data["RPS4Y1",]==0,"sample"]<-"ctrlF541"

cbp@meta.data[cbp@assays$RNA@data["XIST",]==0&
                  cbp@assays$RNA@data["RPS4Y1",]>0,"sample"]<-"lgaM526"

cbp@meta.data[cbp@assays$RNA@data["XIST",]==0&
                  cbp@assays$RNA@data["RPS4Y1",]==0,"sample"]<-"Negative"

table(cbp@meta.data$sample)
# ctrlF541  Doublet  lgaM526 Negative 
# 2438      493     1381      208 



saveRDS(subset(cbp,sample=="ctrlF541"),fp(out,"ctrlF541.rds"))

