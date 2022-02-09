library("Seurat")
source("scripts/utils/new_utils.R")
out<-"outputs/IL6_HSC_stimulation_validation_with_expanded_data"
dir.create(out)

#Is there HSC in expanded ?
expanded<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/expanded_CTRL_CBP547/2020-05-28_seurat_obj1_final.rds")
DimPlot(expanded,label=T)

FeaturePlot(expanded,c("SOCS3","ID1","HES1","EGR1"))
FeaturePlot(expanded,c("CDK6"))
FeaturePlot(expanded,c("DUSP2","DUSP1"))



#map expanded on hmap
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
DefaultAssay(hmap)<-"integrated"
hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

exp_list<-list(expanded,
        readRDS("../singlecell/outputs/01-Analyses_Individuelles/expanded2_CTRL_CBP547/2020-05-29_seurat_obj1_final.rds")
)

exp_list

exp_list<-lapply(exp_list,function(x){
  #message("calculate CC.Difference for",x@project.name)
  if(!"S.Score"%in%colnames(x@meta.data)){
    x<-SCTransform(x,method = "glmGamPoi")
    x <- CellCycleScoring(x,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)


  }
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
  })

exp_list<-lapply(exp_list, SCTransform,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi")


anchors <- list()
for (i in 1:length(exp_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = hmap,
    query = exp_list[[i]],
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )
}


for (i in 1:length(exp_list)) {
  exp_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = exp_list[[i]],
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )
}


# Merge the queries
exp <- merge(exp_list[[1]], exp_list[2:length(exp_list)],merge.dr = c("ref.pca","ref.umap"))
p<-DimPlot(exp, reduction = "ref.umap", group.by =  "predicted.cell_type", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
ggsave(fp(out,"predicted_cell_type.png"),plot=p)

DimPlot(exp, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
