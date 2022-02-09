library(Seurat)
#renv::install("satijalab/seurat-data")
library(SeuratData)
#renv::install("mojaveazure/seurat-disk")
library(SeuratDisk)

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

#convert to anndata
SaveH5Seurat(cbps, filename = "outputs/06-integr_singlecell_cbps/cbps.h5Seurat")
Convert("outputs/06-integr_singlecell_cbps/cbps.h5Seurat", dest = "h5ad")

#light version (with only RNA assay)
cbps_light<-cbps
DefaultAssay(cbps_light)<-"RNA"
cbps_light[["SCT"]]<-NULL
cbps_light[["SNP"]]<-NULL
cbps_light[["prediction.score.cell_type"]]<-NULL
cbps_light[["prediction.score.lineage"]]<-NULL
saveRDS(cbps_light,"outputs/06-integr_singlecell_cbps/cbps_light.rds")
SaveH5Seurat(cbps_light, filename = "outputs/06-integr_singlecell_cbps/cbps_light.h5Seurat")
Convert("outputs/06-integr_singlecell_cbps/cbps_light.h5Seurat", dest = "h5ad")



#cbps_4k version (with only RNA assay)
cbps_4k<-readRDS("outputs/06-integr_singlecell_cbps/cbps_4k.rds")
DefaultAssay(cbps_4k)<-"RNA"
cbps_4k[["SCT"]]<-NULL
cbps_4k[["SNP"]]<-NULL
cbps_4k[["prediction.score.cell_type"]]<-NULL
cbps_4k[["prediction.score.lineage"]]<-NULL
SaveH5Seurat(cbps_4k, filename = "outputs/06-integr_singlecell_cbps/cbps_4k.h5Seurat")
Convert("outputs/06-integr_singlecell_cbps/cbps_4k.h5Seurat", dest = "h5ad")


#light version (with only SCT assay)
cbps_light<-cbps
DefaultAssay(cbps_light)<-"SCT"
cbps_light[["RNA"]]<-NULL
cbps_light[["SNP"]]<-NULL
cbps_light[["prediction.score.cell_type"]]<-NULL
cbps_light[["prediction.score.lineage"]]<-NULL
saveRDS(cbps_light,"outputs/06-integr_singlecell_cbps/cbps_sct_light.rds")
SaveH5Seurat(cbps_light, filename = "outputs/06-integr_singlecell_cbps/cbps_sct_light.h5Seurat")
Convert("outputs/06-integr_singlecell_cbps/cbps_sct_light.h5Seurat", dest = "h5ad")


cbps_atac<-readRDS("../atac/")
cbps_light<-cbps
DefaultAssay(cbps_light)<-"SCT"
cbps_light[["RNA"]]<-NULL
cbps_light[["SNP"]]<-NULL
cbps_light[["prediction.score.cell_type"]]<-NULL
cbps_light[["prediction.score.lineage"]]<-NULL

cbps_light[["SCT"]]@scale.data<-NULL

cbps_tr<-CreateSeuratObject( cbps_light[["SCT"]]@data,
                            meta.data =cbps_light@meta.data,
                            assay = "SCT")
min(cbps_tr@assays$SCT@data)

saveRDS(cbps_tr,"outputs/06-integr_singlecell_cbps/cbps_sct_light.rds")
SaveH5Seurat(cbps_tr,filename = "outputs/06-integr_singlecell_cbps/cbps_sct_light.h5Seurat")
Convert("outputs/06-integr_singlecell_cbps/cbps_sct_light.h5Seurat", dest = "h5ad")


#atac data
mat_a<-readRDS("../atac/outputs/cbps_merged/gene_activities.rds")
dim(mat_a)
mtd_a<-fread("../atac/outputs/cbps_merged/metadata.csv")
unique(mtd_a$dataset)

mat_g<-cbps[["RNA"]]@counts
mtd_g<-data.table(cbps_light@meta.data,keep.rownames = "bc")

inter_genes<-intersect(rownames(mat_a),rownames(mat_g))
length(inter_genes) #18897
#colnames(mat_g)<-paste0(colnames(mat_g),"_rna")
#colnames(mat_a)<-paste0(colnames(mat_a),"_atac")
mat_merge<-cbind(mat_a[inter_genes,mtd_a$bc],mat_g[inter_genes,])

mtd_a[,omics:="atac"]
mtd_g[,omics:="rna"]

mtd_a[,batch:=dataset]

mtd_merge<-rbind(mtd_a,mtd_g,fill=T)

cbps_m<-CreateSeuratObject( mat_merge,
                            meta.data =data.frame(mtd_merge,row.names = "bc"),
                            assay = "RNA")
saveRDS(cbps_m,"outputs/scGANs_batch_correction/cbps_merge.rds")
cbps_m<-readRDS("outputs/scGANs_batch_correction/cbps_merge.rds")

cbps_m@assay

SaveH5Seurat(cbps_m,filename = "outputs/scGANs_batch_correction/cbps_merge.h5Seurat")
Convert("outputs/scGANs_batch_correction/cbps_merge.h5Seurat", dest = "h5ad")

#Hematomap conversion
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

DefaultAssay(hmap)<-"RNA"
hmap[["SCT"]]<-NULL
hmap[["integrated"]]<-NULL

saveRDS(hmap,"outputs/05-make_hematomap/hematomap_light.rds")
hmap<-readRDS("outputs/05-make_hematomap/hematomap_light.rds")
head(hmap[[]])
hmap$cell_type_hmap<-as.character(hmap$cell_type_hmap)
SaveH5Seurat(hmap, filename = "outputs/05-make_hematomap/hematomap_light.h5Seurat")
Convert("outputs/05-make_hematomap/hematomap_light.h5Seurat", dest = "h5ad")

