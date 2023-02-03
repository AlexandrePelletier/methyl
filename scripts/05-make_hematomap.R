source("../methyl/scripts/utils/new_utils.R")
library(Seurat)
options(future.globals.maxSize = 50000 * 1024^2)
sample_name<-"hematomap_ctrls_sans_stress"
out<-"outputs/05-make_hematomap/"
dir.create(out)

#first, need preprocess independently the 7 CBPs CTRL datasets
#run 04A files

hmap_list<-list(ctrlF547=readRDS(fp(out,"ctrlF547.rds"),
              ctrlF544=readRDS(fp(out,"ctrlF544.rds")),
              ctrlF545=readRDS(fp(out,"ctrlF545.rds")),
              ctrlF541=readRDS(fp(out,"ctrlF541.rds")),
              ctrlM555=readRDS(fp(out,"ctrlM555.rds")),
              ctrlM518=readRDS(fp(out,"ctrlM518.rds")),
              ctrlM537=readRDS(fp(out,"ctrlM537.rds")),
              ))


lapply(hmap_list, function(x)head(x@meta.data))

?SCTransform

#Normalisation 
hmap_list<-lapply(hmap_list, SCTransform,
                  return.only.var.genes=F, 
                  method = "glmGamPoi")

#Corriger pour le cycle cellulaire
hmap_list<-lapply(hmap_list,function(x){
  x <- CellCycleScoring(x,s.features = cc.genes$s.genes,
                                 g2m.features = cc.genes$g2m.genes,
                                 set.ident = TRUE,
                                 search=TRUE)
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
  })

# renv::install("bioc::glmGamPoi")
hmap_list<-lapply(hmap_list, SCTransform,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi")

#Alternative : NormalizeData
#?NormalizeData

#Integration


#find common expressed genes
features <- SelectIntegrationFeatures(object.list = hmap_list, nfeatures = 3000)
length(features)
#check normalization in all datasets
hmap_list <- PrepSCTIntegration(object.list = hmap_list, anchor.features = features)

#
hmap_list <- lapply(X = hmap_list, FUN = RunPCA, features = features)

#find link between cell of different library/dataset
hmap.anchors <- FindIntegrationAnchors(object.list = hmap_list, normalization.method = "SCT", 
    anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)

#integrate data
hmap <- IntegrateData(anchorset = hmap.anchors, normalization.method = "SCT", dims = 1:50)

hmap <- RunPCA(hmap,assay="integrated",verbose = FALSE)
hmap <- RunUMAP(hmap, reduction = "pca", dims = 1:50)

hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

# Visualization
DimPlot(hmap, reduction = "umap", group.by = "orig.ident")

p1 <- DimPlot(hmap, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(hmap, reduction = "umap", group.by = "batch", label = TRUE, 
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
DefaultAssay(hmap)<-"integrated"
hmap <- FindNeighbors(object = hmap, dims = 1:50)

hmap<- FindClusters(hmap,resolution = 0.6,
                               algorithm = 4) 
DimPlot(hmap,group.by = "seurat_clusters")
?FindClusters
p1<-DimPlot(hmap,label = T)
p2<-DimPlot(hmap,label = T,group.by = "sample")
p3<-DimPlot(hmap,label = T,group.by = "orig.ident")
p4<-DimPlot(hmap,label = T,group.by = "Phase")
p_all<-(p1+p2)/(p3+p4)
p_all
ggsave(fp(out,"umap_SCT_percent.mt_CC.Difference_rpca_integrated_leiden_res0.6.png"),plot = p_all,width = 10,height = 10)

DimPlot(hmap,label = T,group.by = "seurat_clusters")


#markers identif ####
hmap
Idents(hmap)<-"seurat_clusters"
DefaultAssay(hmap)<-"integrated"

#
markers<-FindMarkers(hmap,ident.1 = "1",ident.2 ="3" ) #pour identifier marqueurs du cluster 1
?FindMarkers

markers<-FindAllMarkers(hmap, min.pct = 0.3, only.pos = TRUE, logfc.threshold = 0.4)
source("../singlecell/scripts/utils/seurat_utils.R")
source("../singlecell/scripts/utils/scoreCluster.R")

#home made function
markers<-scoreMarquageCluster(markers,hmap,seuil = "intraClusterFixe",filtreMin = 2)
markers<-annotMarkers(markers)
#home made function

#fwrite()
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
head(hmap@meta.data)
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

m_lin<-fread("outputs/05-make_hematomap/markers_lineage_annotated.csv.gz")
lapply(split(m_lin,by="cluster"),function(c)c[order(-score,-avg_log2FC)][1:10])
lapply(split(m_lin,by="cluster"),function(dt)dt[order(-score,-avg_log2FC)][MarqueurPops!=""][1:10])



m_lin<-fread("outputs/05-make_hematomap/markers_lineage_annotated.csv.gz")
lapply(split(m_lin,by="cluster"),function(c)c[order(-score,-avg_log2FC)][1:10])
lapply(split(m_lin,by="cluster"),function(dt)dt[order(-score,-avg_log2FC)][MarqueurPops!=""][1:10])



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
