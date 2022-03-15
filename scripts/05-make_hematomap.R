source("../methyl/scripts/utils/new_utils.R")
library(Seurat)
options(future.globals.maxSize = 50000 * 1024^2)
sample_name<-"hematomap_ctrls_sans_stress"
out<-"outputs/05-make_hematomap/"
dir.create(out)

hmap_list<-list(ctrl0=readRDS("../singlecell/outputs/01-Analyses_Individuelles/CD34_CTRL_CBP547/2020-05-28_seurat_obj2_final.rds"),
              ctrlF544=subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP6-a/cbp6a.rds"),sample=="ctrlF544"),
              ctrlF545=subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP6-b/cbp6b.rds"),sample=="ctrlF545"),
              ctrlF541=subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP6-c/cbp6c.rds"),sample=="ctrlF541"),
              ctrlhmap55=subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/cbp7a/cbp7a_singlet.rds"),sample=="ctrlhmap55"),
              ctrlhmap18=subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/cbp7b/cbp7b_singlet.rds"),sample=="ctrlhmap18"),
              ctrlhmap37=readRDS("../singlecell/outputs/01-Analyses_Individuelles/cbp7c/cbp7c_singlet.rds")
              )

lapply(hmap_list, function(x)head(x@meta.data))
hmap_list<-lapply(hmap_list,function(x){
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
  })

# renv::install("bioc::glmGamPoi")
cbps_list<-lapply(cbps_list, SCTransform,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi")

features <- SelectIntegrationFeatures(object.list = cbps_list, nfeatures = 3000)
length(features)
cbps_list <- PrepSCTIntegration(object.list = cbps_list, anchor.features = features)
cbps_list <- lapply(X = cbps_list, FUN = RunPCA, features = features)
hmap.anchors <- FindIntegrationAnchors(object.list = cbps_list, normalization.method = "SCT", 
    anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
hmap <- IntegrateData(anchorset = hmap.anchors, normalization.method = "SCT", dims = 1:50)
hmap <- RunPCA(hmap, verbose = FALSE)
hmap <- RunUMAP(hmap, reduction = "pca", dims = 1:50)
# Visualization
p1 <- DimPlot(hmap, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(hmap, reduction = "umap", group.by = "cell_type", label = TRUE, 
    repel = TRUE)
p1 + p2



rm(hmap.anchors)
hmap@meta.data[is.na(hmap@meta.data$sample),"sample"]<-hmap@meta.data[is.na(hmap@meta.data$sample),"orig.ident"]
hmap@meta.data[hmap@meta.data$sample=="CD34_CTRL_CBP547","sample"]<-"ctrlF547"
hmap@meta.data[hmap@meta.data$sample=="CD34_LGA_CBP552","sample"]<-"lgaF552"

table(hmap@meta.data$sample)
hmap[["group"]]<-str_extract(hmap@meta.data$sample,"ctrl|lga|iugr")
hmap[["sex"]]<-str_extract(hmap@meta.data$sample,"M|F")
hmap[["group_sex"]]<-paste0(hmap@meta.data$group,hmap@meta.data$sex)


saveRDS(hmap,fp(out,paste0(sample_name,".rds")))
hmap<-readRDS(fp(out,paste0(sample_name,".rds")))

#clustering
hmap <- FindNeighbors(object = hmap, dims = 1:50)

hmap<- FindClusters(hmap,resolution = 0.6,
                               algorithm = 4) 
p1<-DimPlot(hmap,label = T)
p2<-DimPlot(hmap,label = T,group.by = "sample")
p3<-DimPlot(hmap,label = T,group.by = "orig.ident")
p4<-DimPlot(hmap,label = T,group.by = "Phase")
p_all<-(p1+p2)/(p3+p4)
p_all
ggsave(fp(out,"umap_SCT_percent.mt_CC.Difference_rpca_integrated_leiden_res0.6.png"),plot = p_all,width = 10,height = 10)

DimPlot(hmap,label = T,group.by = "seurat_clusters")
#markers identif 
Idents(hmap)<-"seurat_clusters"
DefaultAssay(hmap)<-"integrated"
markers<-FindAllMarkers(hmap, min.pct = 0.3, only.pos = TRUE, logfc.threshold = 0.4)
source("../singlecell/scripts/utils/seurat_utils.R")
source("../singlecell/scripts/utils/scoreCluster.R")
markers<-scoreMarquageCluster(markers,hmap,seuil = "intraClusterFixe",filtreMin = 2)
markers<-annotMarkers(markers)
fwrite()
View(markers)
feat1<-c("ID1","EGR1","IRF1","GATA2","GATA1","MPO","KLF2","LTB","VPREB1")
FeaturePlot(hmap,feat1)

head(markers[cell_type=="GMP"],100)

feat2<-c("ID1","ID2","DUSP2", #LT-HSC
           "EGR1","AVP", #HSC-1
         "CXCL8","NFKBIA","ZFP36", #HSC-2
           "MLLT3","CDK6", #MPP
           "SELL","CD99", #LMPP
           "LTB", #CLP
         "VPREB1","IGLL1", # proB
           "IGHM","CD37", #B cell
           "KLF2","TSC22D3", #
           "TNFAIP3","CD7", #T cell
           "IRF1","STAT1", #MkP
           "GATA2", #MPP-Ery
           "GATA1", #EMP
           "BIRC5","MKI67","TOP2A", #ErP-cycle
           "HDC","TFRC","BLVRB","KLF1", #Mast
           "PLEK","HBD", #Mk/Er
           "MPO","CEBPA","CTSG","AZU1", #GMP
         "CST3","CD83") #DC
DotPlot(hmap,features = feat2)

Idents(hmap)<-"seurat_clusters"
hmap<-RenameIdents(hmap,
             "15"="LT-HSC",
             "16"="proB",
             "8"="CLP",
             "1"="LMPP",
             "9"="HSC-4",
             "5"="GMP",
             "13"="GMP-cycle",
             "11"="ErP-cycle",
             "7"="EMP-cycle",
             "6"="EMP",
             "17"="T cell",
             "19"="DC",
             "14"="B cell",
             "20"="Mk/Er",
             "18"="18",
             "12"="HSC-3",
             "3"="MPP-Ery",
             "2"="HSC-1",
             "4"="MPP",
             "10"="HSC-2")
DimPlot(hmap,label = T)
DimPlot(hmap,label = T,cells.highlight = WhichCells(hmap,idents = "HSC-4"))

hmap[["cell_type"]]<-Idents(hmap)
Project(hmap)<-"hmap"
hmap[[paste("cell_type",hmap@project.name,sep="_")]]<-hmap$cell_type
head(hmap[[]])

mtd<-data.table(hmap@meta.data,keep.rownames = "bc")
ct<-unique(mtd[,.(seurat_clusters,cell_type)])
ct[,cluster:=as.factor(seurat_clusters)]
markers<-merge(markers,ct[,.(cluster,cell_type)])





#hSC-2 really HSC ?
DefaultAssay(hmap)<-"SCT"
head(markers[cell_type=="HSC-CXCL8"],20)
FeaturePlot(hmap,c("EGR1"),max.cutoff = 2) #yes
FeaturePlot(hmap,c("DUSP2"),max.cutoff = 2) #yes

#HSC-4 are protT because
head(markers[cell_type=="proT"],20)
#based on https://www.nature.com/articles/s41586-019-1652-y
#yes because KLF2 and TSC22D3 express in innate T cell cluster

#HSC-4 are MkP because :
head(markers[cell_type=="MkP"],20)
#yes with STAT3 and its effector IRF1 trigger megakaryopoiesis (https://www.jci.org/articles/view/33010)

#annot lineage based on lineage markers

hmap[["lineage"]]<-sapply(as.character(hmap@meta.data$cell_type), function(ct){
  if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
  else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
  else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
  else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
  else if(ct%in%c("CLP","proB"))return("Lymphoid")
  else return(ct)
  
})
DimPlot(hmap,label = T,group.by="cell_type")
DimPlot(hmap,label = T,group.by="lineage")

lin<-unique(mtd[,.(seurat_clusters,lineage)])
lin[,cluster:=as.numeric(seurat_clusters)]
markers<-merge(markers,lin[,.(cluster,lineage)])
fwrite(markers,fp(out,paste0(sample_name,"SCT_Leiden_res0.6_markers.csv.gz")),sep=";")

Idents(hmap)<-"lineage"
m_lin<-FindAllMarkers(hmap, min.pct = 0.3, only.pos = TRUE, logfc.threshold = 0.4)
m_lin2<-scoreMarquageCluster(m_lin,hmap,seuil = "intraClusterFixe",filtreMin = 0)
m_lin3<-annotMarkers(m_lin2)
m_lin3[cluster=='MPP/LMPP'&MarqueurPops!=""]
fwrite(m_lin3,fp(out,"markers_lineage_annotated.csv.gz"))
hmap<-RunUMAP(hmap,dims=1:50,
                   reduction.name="ref.umap",
                   reduction.key = "refUMAP_",
                   return.model=TRUE,
                   )

saveRDS(hmap,fp(out,paste0(sample_name,".rds")))

#compute the first 50 neighbors in the PCA space of the reference.
# store this information in the spca.annoy.neighbors object within the reference Seurat object
#and also cache the annoy index data structure (via cache.index = TRUE)

hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
hmap
DefaultAssay(hmap)<-"integrated"

hmap <- FindNeighbors(
  object = hmap,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

SaveAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = fp(out,"reftmp.idx"))
