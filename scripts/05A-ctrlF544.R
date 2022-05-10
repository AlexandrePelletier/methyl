out<-"outputs/05-make_hematomap/"
library(Seurat)
source("scripts/utils/new_utils.R")

batch<-"CBP6-a"
matrix_path<-"~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-a/outs/filtered_feature_bc_matrix/"

#CBP6a contain CBP544 (CTRL F) et CBP556 (LGA M)
cbp<-Read10X(matrix_path)

cbp<-CreateSeuratObject(cbp,project = batch)
cbp #33538 features across 9590

cbp<-PercentageFeatureSet(cbp,pattern = "MT-",col.name = "percent.mt")
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt")) #weird there is 2 groups of cells ~ number of genes : at ~2k and at ~500 genes

cbp<-subset(cbp,nCount_RNA<20000&nFeature_RNA<4000&percent.mt<25)
cbp #33538 features across 8862 

#DEMultiplex
VlnPlot(cbp,c("XIST","RPS4Y1"),group.by = "orig.ident")
FeaturePlot(cbp,c("XIST","RPS4Y1"))

sum(cbp@assays$RNA@data["XIST",]>0) #2555 female
sum(cbp@assays$RNA@data["RPS4Y1",]>0) #3734 male
sum(cbp@assays$RNA@data["XIST",]>0&cbp@assays$RNA@data["RPS4Y1",]>0) #345 doublet
sum(cbp@assays$RNA@data["XIST",]==0&cbp@assays$RNA@data["RPS4Y1",]==0) #2918 negative


cbp@meta.data[cbp@assays$RNA@data["XIST",]>0&
                  cbp@assays$RNA@data["RPS4Y1",]>0,"sample"]<-"Doublet"

cbp@meta.data[cbp@assays$RNA@data["XIST",]>0&
                  cbp@assays$RNA@data["RPS4Y1",]==0,"sample"]<-"ctrlF544"

cbp@meta.data[cbp@assays$RNA@data["XIST",]==0&
                  cbp@assays$RNA@data["RPS4Y1",]>0,"sample"]<-"lgaM556"

cbp@meta.data[cbp@assays$RNA@data["XIST",]==0&
                  cbp@assays$RNA@data["RPS4Y1",]==0,"sample"]<-"Negative"

table(cbp@meta.data$sample)
# ctrlF544  Doublet  lgaM556 Negative 
# 2210      345     3389     2918 

VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt"),group.by = "sample") 
#there is still 2 groups of cells ~ number of genes : at ~2k and at ~500 genes, in the 2 sample
#+ no diff for percen.mt between group

saveRDS(subset(cbp,sample=="ctrlF544"),fp(out,"ctrlF544.rds"))

