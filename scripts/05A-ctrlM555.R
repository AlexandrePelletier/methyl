out<-"outputs/05-make_hematomap/"
library(Seurat)

batch<-"CBP7a"
matrix_path<-"~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-a/outs/filtered_feature_bc_matrix/"

#CBP7a contain CBP554 (LGA F) et 555 (CTRL M)
cbp<-Read10X(matrix_path)

cbp<-CreateSeuratObject(cbp,project = batch)
cbp<-PercentageFeatureSet(cbp,pattern = "MT-",col.name = "percent.mt")
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt")) 
cbp<-subset(cbp,nCount_RNA<30000&nFeature_RNA<5000&percent.mt<25)

#normalisation
cbp <-SCTransform(cbp)

#DEMultiplex

VlnPlot(cbp,c("XIST","RPS4Y1"),group.by = "orig.ident")
FeaturePlot(cbp,c("XIST","RPS4Y1"))

cbp@meta.data[cbp@assays$SCT@data["XIST",]>1&
                  cbp@assays$SCT@data["RPS4Y1",]>1,"sex.ID"]<-"Doublet"

cbp@meta.data[cbp@assays$SCT@data["XIST",]>1&
                  cbp@assays$SCT@data["RPS4Y1",]<=1,"sex.ID"]<-"lgaF554"

cbp@meta.data[cbp@assays$SCT@data["XIST",]<=1&
                  cbp@assays$SCT@data["RPS4Y1",]>1,"sex.ID"]<-"ctrlM555"

cbp@meta.data[cbp@assays$SCT@data["XIST",]<=1&
                  cbp@assays$SCT@data["RPS4Y1",]<=1,"sex.ID"]<-"Negative"

table(cbp@meta.data$sex.ID)
# ctrlM555  Doublet  lgaF554 Negative 
#     1978      332     3416      555 

cbp[["sample"]]<-cbp[["sex.ID"]]

saveRDS(subset(cbp,sample=="ctrlM555"),fp(out,"ctrlM555.rds"))

