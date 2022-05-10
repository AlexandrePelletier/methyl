out<-"outputs/05-make_hematomap/"
library(Seurat)
source("scripts/utils/new_utils.R")

batch<-"CBP7b"
matrix_path<-"~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-b/outs/filtered_feature_bc_matrix/"

#CBP7b contain CBP521 (LGA F) et 518 (CTRL M)
sample_female<-"lgaF521"
sample_male<-"ctrlM518"


cbp<-Read10X(matrix_path)

cbp<-CreateSeuratObject(cbp,project = batch)
cbp<-PercentageFeatureSet(cbp,pattern = "MT-",col.name = "percent.mt")
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt")) 
cbp<-subset(cbp,nCount_RNA<30000&nFeature_RNA<5000&percent.mt<15)

cbp #33538 features across 6594  samples
#normalisation

cbp <-SCTransform(cbp)

#DEMultiplex

FeaturePlot(cbp,c("XIST","RPS4Y1"))

sum(cbp@assays$SCT@data["XIST",]>0&cbp@assays$SCT@data["RPS4Y1",]>0) #1975 doublet
sum(cbp@assays$SCT@data["XIST",]==0&cbp@assays$SCT@data["RPS4Y1",]==0) #67 negative

#a lot of doublet, check if bimodale (vrai doublet et background)

plot(density(cbp@assays$SCT@data["XIST",cbp@assays$SCT@data["XIST",]>0]))
abline(v=1.25)


sum(cbp@assays$SCT@data["XIST",]<1.25&cbp@assays$SCT@data["XIST",]>0&cbp@assays$SCT@data["RPS4Y1",]>0)
sum(cbp@assays$SCT@data["XIST",]<1.25&cbp@assays$SCT@data["XIST",]>0) #1083/1430, so its background

plot(density(cbp7a@assays$SCT@data["RPS4Y1",cbp7a@assays$SCT@data["RPS4Y1",]>0]))
abline(v=1)
cbp@assays$SCT@misc$xist_thres<-1.25
cbp@assays$SCT@misc$rps4y_thres<-1

cbp@meta.data[cbp@assays$SCT@data["XIST",]>cbp@assays$SCT@misc$xist_thres&
                  cbp@assays$SCT@data["RPS4Y1",]>cbp@assays$SCT@misc$rps4y_thres,"sex.ID"]<-"Doublet"
cbp@meta.data[cbp@assays$SCT@data["XIST",]>cbp@assays$SCT@misc$xist_thres&
                  cbp@assays$SCT@data["RPS4Y1",]<=cbp@assays$SCT@misc$rps4y_thres,"sex.ID"]<-sample_female

cbp@meta.data[cbp@assays$SCT@data["XIST",]<=cbp@assays$SCT@misc$xist_thres&
                  cbp@assays$SCT@data["RPS4Y1",]>cbp@assays$SCT@misc$rps4y_thres,"sex.ID"]<-sample_male

cbp@meta.data[cbp@assays$SCT@data["XIST",]<=cbp@assays$SCT@misc$xist_thres&
                  cbp@assays$SCT@data["RPS4Y1",]<=cbp@assays$SCT@misc$rps4y_thres,"sex.ID"]<-"Negative"

table(cbp@meta.data$sex.ID)
# ctrlM518  Doublet  lgaF521 Negative 
# 2042      161     3506      885 

VlnPlot(cbp,c("nCount_RNA","nFeature_RNA","percent.mt"),group.by = "sex.ID") 
# no diff for percen.mt between group
VlnPlot(cbp,c("nCount_RNA","nFeature_RNA"),group.by = "sex.ID",log = T) 


DimPlot(cbp,label=T,group.by = "sex.ID")+p2
cbp[["sample"]]<-cbp[["sex.ID"]]

cbp_s<-subset(cbp,sample%in%c(sample_female,sample_male))


saveRDS(subset(cbp,sample==sample_male),fp(out,ps(sample_male,".rds")))

