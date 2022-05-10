out<-"outputs/05-make_hematomap/"
library(Seurat)
source("scripts/utils/new_utils.R")


#CTRL CBP547
#loading
data <- read.csv("../singlecell/datasets/CTRL_CD34/CBP547.csv",row.names = 1)
data<-data[rowSums(data)>30,]
dim(data) #12008  6965
head(data[,1:10])

#translation in hgnc_symbolID
library(data.table)
ref<-fread("ref/ENSEMBL_ID_to_SYMBOL.csv")
data$hgnc_symbol<-ref$hgnc_symbol[match(rownames(data),ref$ensembl_gene_id)]
data<-data[!duplicated(data$hgnc_symbol),]
data<-data[!is.na(data$hgnc_symbol),]
rownames(data)<-data$hgnc_symbol
data<-data[rownames(data)!="",]
data<-data[,-c("hgnc_symbol")]
dim(data) #7517 6965
head(data[,(ncol(data)-5):ncol(data)])

#Seurat object
sample <- CreateSeuratObject(counts = data, project = "CD34_CTRL_CBP547", min.cells = 3, min.features = 200)
sample #7517 features across 6925 samples


#QC / degraded or doublet cells removal
library(ggplot2)
library(patchwork)

sample[["percent.mt"]] <- PercentageFeatureSet(object = sample, pattern = "^MT-")
head(x = sample@meta.data, 5)

vps<-VlnPlot(object = sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vp1<-vps[[1]]+geom_hline(yintercept = 2800,colour="red")
vp2<-vps[[2]]+geom_hline(yintercept = 12000,colour="red")
vp3<-vps[[3]]+geom_hline(yintercept = 9,colour="red")

vp1+vp2+vp3

sample <- subset(sample,  nFeature_RNA < 2800 & nCount_RNA<12000 & percent.mt < 9) 
sample #7517 features across 6428 samples

saveRDS(sample,fp(out,"ctrlF547.rds"))

