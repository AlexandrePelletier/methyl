library(Seurat)
#renv::install("satijalab/seurat-data")
library(SeuratData)
#renv::install("mojaveazure/seurat-disk")
library(SeuratDisk)

system("wget -O outputs/Tabula_Sapiens.h5ad.zip https://figshare.com/ndownloader/files/28846899")
system("unzip outputs/Tabula_Sapiens.h5ad.zip -d outputs")

Convert("outputs/TS_Skin.h5ad", dest = "h5seurat", overwrite = TRUE)
tabula <- LoadH5Seurat("outputs/TS_Skin.h5seurat",assays='raw_counts')
tabula
head(tabula[[]])
DimPlot(tabula,group.by="compartment",label=T)+NoLegend()
