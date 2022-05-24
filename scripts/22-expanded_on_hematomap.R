source("scripts/utils/new_utils.R")
library(Seurat)

out<-"outputs/22-expanded_on_hematomap"
dir.create(out)

hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
DefaultAssay(hmap)<-"integrated"
hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

#expanded 

expanded<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/expanded_CTRL_CBP547/2020-05-28_seurat_obj1_final.rds")

expanded<-SCTransform(expanded,method = "glmGamPoi")
expanded <- CellCycleScoring(expanded,s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

expanded$CC.Difference <- expanded$S.Score - expanded$G2M.Score
expanded<-SCTransform(expanded,method = "glmGamPoi",
                      vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F)

anchors <- FindTransferAnchors(
    reference = hmap,
    query = expanded,
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )

expanded <- MapQuery(
    anchorset = anchors, 
    query = expanded,
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )

DimPlot(expanded, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
saveRDS(expanded,fp(out,"expanded1.rds"))
#exp 3438
exp3438<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/expandCD34CD38_CTRL_CBP547/2020-05-29_seurat_obj1_final.rds")

exp3438<-SCTransform(exp3438,method = "glmGamPoi")
exp3438 <- CellCycleScoring(exp3438,s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

exp3438$CC.Difference <- exp3438$S.Score - exp3438$G2M.Score
exp3438<-SCTransform(exp3438,method = "glmGamPoi",
                      vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F)

anchors <- FindTransferAnchors(
    reference = hmap,
    query = exp3438,
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )

exp3438 <- MapQuery(
    anchorset = anchors, 
    query = exp3438,
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )

DimPlot(exp3438, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
p1<-DimPlot(exp3438, reduction = "ref.umap", group.by =  "predicted.cell_type", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
p2<-FeaturePlot(exp3438, c("GATA1","MPO","VPREB1","CD99"),reduction = "ref.umap") + NoLegend()
p1+p2
p2<-FeaturePlot(exp3438, c("GATA1","HBD","LTB","AZU1"),reduction = "ref.umap") + NoLegend()
p1+p2

saveRDS(exp3438,fp(out,"expanded3438.rds"))
Idents(exp3438)<-"predicted.cell_type"
erPCycle_markers<-FindMarkers(exp3438,ident.1 = "ErP-cycle",only.pos = T)
erPCycle_markers #markers de prolif + HMGB1, HSProt, LDHA et B, BIRC5, KIF23...
#cellules ++ differencié ? 

#cluster late branched, because are dividing cells ?
FeaturePlot(exp3438, c("AVP","SOCS3","ID1","HES1"),reduction = "ref.umap") 
FeaturePlot(exp3438, c("JUNB","ID2","EGR1","HES1"),reduction = "ref.umap") 
