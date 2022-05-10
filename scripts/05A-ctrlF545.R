out<-"outputs/05-make_hematomap/"
library(Seurat)
source("scripts/utils/new_utils.R")


batch<-"CBP6-b"
matrix_path<-"~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-b/outs/filtered_feature_bc_matrix/"

#CBP6b contain CBP545 (CTRL F) et 533 (LGA M)
cbp<-Read10X(matrix_path)

cbp<-CreateSeuratObject(cbp,project = batch)


cbp<-PercentageFeatureSet(cbp,pattern = "MT-",col.name = "percent.mt")
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt")) 
cbp<-subset(cbp,nCount_RNA<20000&nFeature_RNA<4000&nFeature_RNA>200&percent.mt<25)
cbp 
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
                  cbp@assays$RNA@data["RPS4Y1",]==0,"sample"]<-"ctrlF545"

cbp@meta.data[cbp@assays$RNA@data["XIST",]==0&
                  cbp@assays$RNA@data["RPS4Y1",]>0,"sample"]<-"lgaM533"

cbp@meta.data[cbp@assays$RNA@data["XIST",]==0&
                  cbp@assays$RNA@data["RPS4Y1",]==0,"sample"]<-"Negative"

table(cbp@meta.data$sample)
# ctrlF545  Doublet  lgaM533 Negative 
# 798      151     2924       94


saveRDS(subset(cbp,sample=="ctrlF545"),fp(out,"ctrlF545.rds"))

