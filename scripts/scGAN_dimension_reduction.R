library("Seurat")
source("scripts/utils/new_utils.R")

out<-"outputs/scGANs_batch_correction"
dir.create(out)
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

#without batch effect correction
#RNA
DefaultAssay(cbps)<-"RNA"
#VariableFeatures(cbps)<-Reduce(union,lapply(SplitObject(cbps,split.by = "orig.ident"),function(x)VariableFeatures(FindVariableFeatures(x))))
cbps<-FindVariableFeatures(cbps)
cbps<-NormalizeData(cbps)
cbps<-ScaleData(cbps)
cbps<-RunPCA(cbps,assay = "RNA",reduction.name = "pca_null")

# cbps<-FindNeighbors(cbps,dims = 1:15,reduction = "pca_null",graph.name = c("nn","z_snn"))
# cbps<-FindClusters(cbps,resolution = 0.6,graph.name = "z_snn")
# cbps$z_rna_clusters<-Idents(cbps)
cbps<-RunUMAP(cbps,reduction = "pca_null",dims = 1:30,
              reduction.name = "umap_null",reduction.key = "nUMAP_" )

DimPlot(cbps,reduction = "umap_null",group.by = c("batch"))
DimPlot(cbps,reduction = "umap_null",group.by = c("lineage_hmap"))
DimPlot(cbps,reduction = "umap_null",group.by = "hto")
DimPlot(cbps,reduction = "umap_null",group.by = "group")
FeaturePlot(cbps,c("MPO","LTB"),max.cutoff = "q95",reduction = "zumap")

#SCT
DefaultAssay(cbps)<-"SCT"
#VariableFeatures(cbps)<-Reduce(union,lapply(SplitObject(cbps,split.by = "orig.ident"),function(x)VariableFeatures(FindVariableFeatures(x))))
cbps<-RunPCA(cbps,assay = "SCT",reduction.name = "pca_sct",features = rownames(cbps))

# cbps<-FindNeighbors(cbps,dims = 1:15,reduction = "pca_null",graph.name = c("nn","z_snn"))
# cbps<-FindClusters(cbps,resolution = 0.6,graph.name = "z_snn")
# cbps$z_rna_clusters<-Idents(cbps)
cbps<-RunUMAP(cbps,reduction = "pca_null",dims = 1:30,
              reduction.name = "umap_null",reduction.key = "nUMAP_" )

DimPlot(cbps,reduction = "umap_null",group.by = c("batch"))
DimPlot(cbps,reduction = "umap_null",group.by = c("lineage_hmap"),label = T)
DimPlot(cbps,reduction = "umap_null",group.by = "hto")
DimPlot(cbps,reduction = "umap_null",group.by = "group")
FeaturePlot(cbps,c("MPO","LTB"),max.cutoff = "q95",reduction = "zumap")



#with seurat based batxh effect corr
DimPlot(cbps,reduction = "ref.umap",group.by = c("batch"))
DimPlot(cbps,reduction = "ref.umap",group.by = c("lineage_hmap"),label = T)

#scGAN RNA 15 embeddings
z<-fread("outputs/scGANs_batch_correction/z_15_embeddings.csv",header = T,select = 2:17,
         col.names = c(paste0("z_",1:15),"cell_label"))

z_mat<-as.matrix(data.frame(z,row.names = "cell_label"))
colnames(z_mat)<-as.numeric(str_extract_all(colnames(z_mat),"[0-9]+"))

cbps[["z_rna"]]<-CreateDimReducObject(z_mat,key = "z_")
DefaultAssay(cbps)<-"RNA"
cbps<-FindNeighbors(cbps,dims = 1:15,reduction = "z_rna",graph.name = c("z_nn","z_snn"))
cbps<-FindClusters(cbps,resolution = 0.6,graph.name = "z_snn")
cbps$z_rna_clusters<-Idents(cbps)
cbps<-RunUMAP(cbps,reduction = "z_rna",dims = 1:15,
              reduction.name = "zumap",reduction.key = "zUMAP_" )

DimPlot(cbps,reduction = "zumap",group.by = "batch")
DimPlot(cbps,reduction = "zumap",group.by = c("lineage_hmap",
                                             "z_rna_clusters"),label = T)
DimPlot(cbps,reduction = "zumap",group.by = "hto")
DimPlot(cbps,reduction = "zumap",group.by = "group")
FeaturePlot(cbps,c("MPO","LTB"),max.cutoff = "q95",reduction = "zumap")


#test all cells with 15 embeddings sct norm
z<-fread("outputs/scGANs_batch_correction/z_embeddings_sct.csv",header = T,select = 2:17,
         col.names = c(paste0("z_",1:15),"cell_label"))
dim(z)
z_mat<-as.matrix(data.frame(z,row.names = "cell_label"))
colnames(z_mat)<-as.numeric(str_extract_all(colnames(z_mat),"[0-9]+"))
#cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps[["z_sct"]]<-CreateDimReducObject(z_mat,key = "z_")
DefaultAssay(cbps)<-"SCT"
cbps<-FindNeighbors(cbps,dims = 1:15,reduction = "z_sct",graph.name = c("z_nn","z_snn"))
cbps<-FindClusters(cbps,resolution = 0.6,graph.name = "z_snn")
cbps$z_snn_clusters<-Idents(cbps)
cbps<-RunUMAP(cbps,reduction = "z_sct",dims = 1:15,
              reduction.name = "zumap",reduction.key = "zUMAP_" )

DimPlot(cbps,reduction = "zumap",group.by = "batch")
DimPlot(cbps,reduction = "zumap",group.by = c("lineage_hmap",
                                             "z_snn_clusters"),label = T)
DimPlot(cbps,reduction = "zumap",group.by = "hto")
DimPlot(cbps,reduction = "zumap",group.by = "group")
FeaturePlot(cbps,c("MPO","LTB"),reduction = "zumap")

#MPO and LTB / Myeloid and Lymphoid are mixed 

#hmap scGAN
z<-fread("outputs/scGANs_batch_correction/z_embeddings_hmap.csv",header = T,select = 2:17,
         col.names = c(paste0("z_",1:15),"cell_label"))
dim(z)
z_mat<-as.matrix(data.frame(z,row.names = "cell_label"))
colnames(z_mat)<-as.numeric(str_extract_all(colnames(z_mat),"[0-9]+"))
hmap<-readRDS("outputs/05-make_hematomap/hematomap_light.rds")

hmap[["z_rna"]]<-CreateDimReducObject(z_mat,key = "z_")
DefaultAssay(hmap)<-"RNA"
hmap<-FindNeighbors(hmap,dims = 1:15,reduction = "z_rna",graph.name = c("z_nn","z_snn"))
hmap<-FindClusters(hmap,resolution = 0.6,graph.name = "z_snn")
hmap$z_snn_clusters<-Idents(hmap)
hmap<-RunUMAP(hmap,reduction = "z_rna",dims = 1:15,
              reduction.name = "zumap",reduction.key = "zUMAP_" )

DimPlot(hmap,reduction = "zumap",group.by = "orig.ident")
DimPlot(hmap,reduction = "zumap",group.by = c("lineage",
                                             "z_snn_clusters"),label = T)

FeaturePlot(hmap,c("MPO","LTB"),reduction = "zumap",max.cutoff = "q95")

#calculate Batch entropy for all
#Batch entropy(BE) and clustering quality based on the Adjusted Rand Index (ARI)
#BE 50 nearest neighbors of the sampled cell based on the cell embeddings (which
#vary depending on the method).We measured the batch label frequency of
#the 50 nearest neighbor cells,


BatchEntropy<-function(seurat_obj,reduction,batch="batch",dims=1:30,n.neighb=50,n.sampling=100,n.repet=50){
  shannon_entropy<-function(batch_fq,n_batches){
    return(-sum(batch_fq * log(batch_fq))/log(n_batches))

  }
  
  seurat_obj<-FindNeighbors(seurat_obj,
                    dims = dims,
                    reduction = reduction,
                    k.param = n.neighb,
                    graph.name='nn',
                    return.neighbor = T,
                    compute.SNN = F)
  
  n_batches<-length(unique(seurat_obj@meta.data[,batch]))
  
  return(mean(sapply(1:n.repet,function(x){
    set.seed(x)
    cells_idx<-sample(1:ncol(seurat_obj),size = n.sampling)
    
    mean(sapply(cells_idx,function(c){
      ngb_idx<-seurat_obj@neighbors$nn@nn.idx[c,]
      batch_fq<-table(seurat_obj@meta.data[ngb_idx,batch])/length(ngb_idx)
    return(shannon_entropy(batch_fq,n_batches))
      }))
    
    })))
  
  
  }

res_be<-data.table(reduction=c('pca_null','pca_sct','ref.pca','z_rna','z_sct'))

res_be[,batch_entropy:=BatchEntropy(cbps,reduction,dims=1:15),"reduction"]

res_be[,cell_type_entropy:=BatchEntropy(cbps,reduction,batch = "lineage_hmap",dims=1:15),"reduction"]

p1<-ggplot(res_be)+geom_col(aes(x=reduction,y=batch_entropy,fill=reduction))
p2<-ggplot(res_be)+geom_col(aes(x=reduction,y=cell_type_entropy,fill=reduction))
p1+p2+plot_layout(guides = "collect")

fwrite(res_be,fp(out,"res_batch_entropy180122.csv"))
