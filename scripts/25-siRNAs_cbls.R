#sirnas experiments

out<-"outputs/25-siRNAs_KLF2_cbl"
dir.create(out)
library(Seurat)
source("scripts/utils/new_utils.R")

#siCTRL
mat<-Read10X("~/RUN/Run_727_single-cell/output/count/humain/sict/outs/filtered_feature_bc_matrix/")

sict<-CreateSeuratObject(mat,project = "siCTRL_cbls")
sict #1919 cells
sict[["percent.mt"]]<-PercentageFeatureSet(sict,pattern = "MT-")
VlnPlot(object = sict, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

sict<-subset(sict,nFeature_RNA>2500&nCount_RNA>10000&percent.mt<8)
sict#1533 cells

#siKLF2
mat<-Read10X("~/RUN/Run_727_single-cell/output/count/humain/siklf2/outs/filtered_feature_bc_matrix/")

sikl<-CreateSeuratObject(mat,project = "siKLF2_cbls")
sikl #2836 cells
sikl[["percent.mt"]]<-PercentageFeatureSet(sikl,pattern = "MT-")
VlnPlot(object = sikl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

sikl<-subset(sikl,nFeature_RNA>2500&nCount_RNA>10000&percent.mt<8)
sikl#2116  cells

#merge
sict$sirna<-"CTRL"
sikl$sirna<-"KLF2"

sicbl<-merge(sict,sikl)

#norm
sicbl_list<-SplitObject(sicbl,split.by = "orig.ident")

sicbl_list<-lapply(sicbl_list,"SCTransform")

#KLF2 downregulation?
sicbl<-merge(sicbl_list[[1]],sicbl_list[[2]])
sicbl#3649 feature

VlnPlot(sicbl,features = "KLF2")

#clustering
VariableFeatures(sicbl)<-intersect(VariableFeatures(sicbl_list[[1]]),VariableFeatures(sicbl_list[[2]]))
sicbl<-RunPCA(sicbl)
sicbl<-FindNeighbors(sicbl,dims = 1:30)
sicbl<-RunUMAP(sicbl,dims = 1:30)
sicbl<-FindClusters(sicbl,resolution = 0.6)

sicbl <- CellCycleScoring(sicbl,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)
DimPlot(sicbl,label=T)#cluster ~ cell cycle

head(sicbl[[]])
sicbl_list<-SplitObject(sicbl,split.by = "orig.ident")

sicbl_list<-lapply(sicbl_list, SCTransform,vars.to.regress=c("S.Score","G2M.Score","percent.mt"))
sicbl<-merge(sicbl_list[[1]],sicbl_list[[2]])

#clustering
VariableFeatures(sicbl)<-intersect(VariableFeatures(sicbl_list[[1]]),VariableFeatures(sicbl_list[[2]]))
sicbl<-RunPCA(sicbl)
sicbl<-FindNeighbors(sicbl,dims = 1:30)
sicbl<-RunUMAP(sicbl,dims = 1:30)
DimPlot(sicbl,label=T)#ok
FeaturePlot(sicbl,features = "percent.mt")#ok

sicbl<-FindClusters(sicbl,resolution = 0.6)
DimPlot(sicbl,label=T)#
FeaturePlot(sicbl,features =c("GATA1","VPREB1","LTB","MPO","SOCS3","AVP"))#ok

saveRDS(sicbl,fp(out,"sicbls.rds"))

sicbl<-readRDS(fp(out,"sicbls.rds"))
sicbl #3649

sicbl_list<-SplitObject(sicbl,split.by = "orig.ident")


#label cells based on hmap
library(parallel)
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

DefaultAssay(hmap)<-"integrated"
#hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

sicbl_list<-lapply(sicbl_list,function(x){
    x<-SCTransform(x)
    x <- CellCycleScoring(x,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)


  
  x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
  },mc.cores = 2)

for (i in 1:length(sicbl_list)) {
  sicbl_list[[i]] <- SCTransform(sicbl_list[[i]],vars.to.regress = c("percent.mt","CC.Difference"))
  sicbl_list[[i]] <- RunPCA(sicbl_list[[i]])
}

for (i in 1:length(sicbl_list)) {
  # transfer cell type labels from reference to query
  transfer_anchors <- FindTransferAnchors(
    reference = hmap,
    query = sicbl_list[[i]],
    normalization.method = "SCT",
    reference.reduction = "pca",
    recompute.residuals = FALSE,
    dims = 1:30
  )
  
  predictions <- TransferData(
    anchorset = transfer_anchors, 
    refdata = hmap$cell_type,
    weight.reduction = sicbl_list[[i]][['pca']],
    dims = 1:30
  )
  
  sicbl_list[[i]] <- AddMetaData(
    object = sicbl_list[[i]],
    metadata = predictions
  )
  
  # set the cell identities to the cell type predictions
  Idents(sicbl_list[[i]]) <- "predicted.id"
  
    sicbl_list[[i]]$cell_type_hmap<-Idents(sicbl_list[[i]])
    sicbl_list[[i]][["lineage_hmap"]]<-sapply(as.character(sicbl_list[[i]]@meta.data$cell_type_hmap), function(ct){
      if(ct%in%c("HSC-1","HSC-2","HSC-3","HSC-4"))return("HSC")
      else if(ct%in%c("MPP","MPP-Ery","LMPP"))return("MPP/LMPP")
      else if(ct%in%c("EMP","EMP-cycle","ErP-cycle"))return("Erythro-Mas")
      else if(ct%in%c("GMP","GMP-cycle"))return("Myeloid")
      else if(ct%in%c("CLP","proB"))return("Lymphoid")
      else return(ct)
  })

}

# Merge the queries
sicbls <- merge(sicbl_list[[1]], sicbl_list[[2]])
saveRDS(sicbls,fp(out,"sicbls2.rds"))

# unsupervised clustering
features<-SelectIntegrationFeatures(sicbl_list, nfeatures = 3000)

sicbls <- RunPCA(sicbls, assay = "SCT",features=features)
sicbls <- RunUMAP(sicbls, reduction = 'pca', dims = 1:30,reduction.name = "umap")
DimPlot(sicbls, group.by = c('sirna',"lineage_hmap"),reduction = "umap", label = TRUE)
DimPlot(sicbls, group.by = c('cell_type_hmap',"lineage_hmap"),reduction = "umap", label = TRUE)

FeaturePlot(sicbls,features =c("GATA1","VPREB1","LTB","MPO","SOCS3","AVP"))#ok
FeaturePlot(sicbls,features =c("GATA1","GATA2"))#ok

FeaturePlot(sicbls,features =c("SOCS3","AVP","EGR1","ID1"))#ok
FeaturePlot(sicbls,features =c("AVP","CDK6"))#ok

#harmony ? no
# library(harmony)
# cbls_stim<-RunHarmony(cbls_stim,reduction = "pca",group.by.vars = "condition",assay.use = "SCT")
# cbls_stim <- RunUMAP(cbls_stim, reduction = 'harmony', dims = 1:30,reduction.name = "humap")
# DimPlot(cbls_stim, group.by = c('condition',"lineage_hmap"),reduction = "humap", label = TRUE)

#DimPlot(cbls_stim, group.by = 'cell_type_hmap',reduction = "humap", label = TRUE)
#DimPlot(cbls_stim, group.by = 'lineage_hmap',reduction = "humap", label = TRUE)
table(sicbls$sirna,sicbls$lineage_hmap)
  #     B cell   DC Erythro-Mas  HSC Lymphoid Mk/Er MPP/LMPP Myeloid
  # CTRL      0    2         373   31       40     8      994      85
  # KLF2      1    5         495   34       57    31     1376     117
VlnPlot(sicbls, features=c("EGR1","MPO","LTB","GATA2"), group.by = 'lineage_hmap')

#markers KLF2 with this few HSC ?
Idents(sicbls)<-'lineage_hmap'
mklf2<-FindMarkers(sicbls,group.by = "sirna",ident.1 = "KLF2",subset.ident = "HSC")
fwrite(data.table(mklf2,keep.rownames = "gene"),fp(out,'degs_klf2_hsc_restreint.csv.gz'))


DimPlot(sicbls, reduction = "umap", group.by =  c("Phase","orig.ident")) 
#Phase +
FeaturePlot(sicbls,features =c("AVP","CDK6"),min.cutoff = "q10")#ok
#HSC > cluster du haut
sicbls<-FindNeighbors(sicbls,dims = 1:30,reduction = "pca")

sicbls<-FindClusters(sicbls,resolution = 1.2)

DimPlot(sicbls, reduction = "umap", group.by = c( 'seurat_clusters',"lineage_hmap"), label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
DimPlot(sicbls, reduction = "umap", group.by = c( 'seurat_clusters',"cell_type_hmap"), label = TRUE, repel = TRUE, label.size = 3) 
FeaturePlot(sicbls,features =c("MPO","CDK6"),min.cutoff = "q10")#ok
FeaturePlot(sicbls,features =c("HBD","GATA1","PLEK"),min.cutoff = "q10")#ok

newids<-c(
  "10"="HSC",
  "3"="MPP",
  "2"="MPP-Ery",
  "4"="MPP-Ery",
  "0"="EMP",
  "8"="ErP/GMP-cycle",
  "11"="GMP",
  "1"="LMPP",
  "5"="LMPP",
  "7"="LMPP",
  "6"="LMPP",
  "9"="CLP"
)
sicbls<-RenameIdents(sicbls,newids)
DimPlot(sicbls, reduction = "umap",  label = TRUE, label.size = 4) 
FeaturePlot(sicbls,features =c("MPO","LTB",'GATA2'),min.cutoff = "q10")#ok

sicbls$cell_type2<-Idents(sicbls)


#KLF2 reudction in MPP ?
VlnPlot(sicbls,"KLF2",group.by = "cell_type2") #nop


#diff of lin ?
ggplot(sicbls@meta.data)+geom_bar(aes(x=sirna,fill=cell_type2),position = "fill")

#DEGs between sirnas ?
Idents(sicbls)<-"cell_type2"
res<-FindMarkers(sicbls,group.by = "sirna",ident.1 = "KLF2",ident.2 = "CTRL",subset.ident = "HSC",logfc.threshold = 0)
res<-data.table(res,keep.rownames = "gene")

#volcano
ggplot(res,aes(x=avg_log2FC,y=-log10(p_val),col=p_val_adj<0.05))+
  geom_point(size=1) +
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-1,1))+
  theme_minimal() +
  theme(legend.position = "bottom")

res[p_val_adj<0.05] #32
degs_list<-split(res[p_val_adj<0.05]$gene,res[p_val_adj<0.05]$avg_log2FC>0)
names(degs_list)<-c("down","up")
degs_list
# $down
#  [1] "LRRC75A"  "CTNNB1"   "CD63"     "HNRNPH1"  "CDC42"    "TMSB10"  
#  [7] "PLCG2"    "YWHAE"    "SH3BGRL3" "B2M"      "PRDX6"    "RBX1"    
# [13] "JTB"      "RPL36AL"  "MYL12B"   "ATP1B3"   "SSR4"     "MIF"     
# [19] "PHF5A"    "MYL6"     "VMP1"     "NDUFB8"   "SKP1"     "FAU"     
# [25] "COX6C"    "TXN"      "DDX6"     "PPP1CB"   "PPIA"     "SRSF6"   
# [31] "NORAD"   
# 
# $up
# [1] "CD81"


library(gprofiler2)
res_gos<-gost(query = degs_list$down,organism = "hsapiens")
gprofiler2::gostplot(res_gos) #klf15 motig
resgos<-data.table(res_gos$result)
lapply(split(resgos,by = "source"),head)

#save
saveRDS(sicbls,fp(out,"sicbls2.rds"))
fwrite(degs,fp(out,"degs_siklf2_vs_ctrl_hsc.csv.gz"))
fwrite(resgos,fp(out,"res_gost_degs_siklf2_vs_ctrl_hsc.csv.gz"))

#enrichment for regulons KLF2
regulons<-fread("outputs/16-GRN_final/tf_target_interactions.csv")
intersect(regulons[tf=="KLF2"]$target,res[p_val_adj<0.05]$gene) # "PLCG2" 1/32 degs..
#donc nop

#gsea lik


res_gos<-gost(query = res[order(-log10(p_val)*sign(avg_log2FC))]$gene,organism = "hsapiens",ordered_query = T)
gprofiler2::gostplot(res_gos) 

resgos<-data.table(res_gos$result)
resgos[p_value<1e-15][source=="TF"][1:20]
#gsea with regulons
library(fgsea)
res[,score:=sign(avg_log2FC)*rank(-p_val)]

genes_rank<-res$score
names(genes_rank)<-res$gene
genes_rank<-sort(genes_rank,decreasing = T)
head(genes_rank)
tail(genes_rank)

res_gsea<-fgsea(pathways=split(regulons$target,regulons$tf),
      stats=genes_rank,scoreType = "std")
res_gsea[order(pval)][1:30] #CTCF, ELF1, but alsoKLF2!

res_gsea[,size.regulon:=length(split(regulons$target,regulons$tf)[[pathway]]),by="pathway"]
head(res_gsea[order(padj)],30)

fwrite(res_gsea,fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_sign_p_val.csv.gz"))


#fig
lead_genes<-res_gsea[pathway=="KLF2"]$leadingEdge[[1]]
length(lead_genes)#63
regulons[tf=="KLF2"]#188
res[gene%in%lead_genes]


ggplot(res,aes(x=avg_log2FC,y=-log10(p_val),col=gene%in%lead_genes,alpha=gene%in%lead_genes))+
  geom_point(size=1) +
  scale_color_manual(values = c("grey","red")) +
  scale_x_continuous(limits=c(-1,1))+
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(res_gsea,aes(x=NES,y=-log10(pval),col=padj<0.05))+
  geom_point(aes(size=size))+
  scale_color_manual(values = c("grey","red")) +
  geom_label_repel(aes(label=ifelse(padj<0.05,pathway,"")))

#most statistically different
res[,score:=rank(-log10(p_val))]

genes_rank<-res$score
names(genes_rank)<-res$gene
genes_rank<-sort(genes_rank,decreasing = T)
res_gsea<-fgsea(pathways=split(regulons$target,regulons$tf),
      stats=genes_rank,scoreType = "pos")
res_gsea[order(pval)][1:10] #KLF2 in top5!
lead_genes<-res_gsea[pathway=="KLF2"]$leadingEdge[[1]]
length(lead_genes)#63
regulons[tf=="KLF2"]#188
res[gene%in%lead_genes]
#        gene        p_val  avg_log2FC pct.1 pct.2    p_val_adj  score
#  1:    PLCG2 3.568329e-08 -0.55100628 1.000 1.000 0.0006284185 9886.0
#  2:     ARF6 1.380192e-04 -0.57357839 0.688 0.877 1.0000000000 9767.0
#  3:    RBM39 4.257180e-04 -0.36311841 1.000 1.000 1.0000000000 9697.0
#  4:      UBC 4.371398e-03 -0.54468093 1.000 1.000 1.0000000000 9409.0
#  5:     PIM3 9.075945e-03 -0.27463334 0.584 0.754 1.0000000000 9230.0
#  6:     EIF1 9.245998e-03 -0.29452590 1.000 1.000 1.0000000000 9224.0
#  7:    HLA-C 9.907604e-03 -0.29355551 1.000 1.000 1.0000000000 9195.0
#  8:     RHOC 1.008053e-02 -0.47747807 0.545 0.692 1.0000000000 9192.0
#  9:   TIPARP 1.144470e-02 -0.37914383 0.831 0.938 1.0000000000 9151.0
# 10:    HLA-E 1.182079e-02 -0.31530221 1.000 1.000 1.0000000000 9139.0
# 11:    ITM2B 2.297621e-02 -0.19986792 1.000 1.000 1.0000000000 8854.0
# 12:     KLF7 2.804374e-02 -0.29071238 0.117 0.262 1.0000000000 8744.0
# 13:    SRSF5 3.498782e-02 -0.20521728 1.000 1.000 1.0000000000 8607.0
# 14:    RAB21 4.437118e-02 -0.37111073 0.701 0.785 1.0000000000 8430.0
# 15:     BTG1 4.707513e-02 -0.25927669 0.935 0.954 1.0000000000 8374.0
# 16:    KLF10 5.571895e-02 -0.36198832 0.234 0.354 1.0000000000 8258.0
# 17:     LMNA 5.767270e-02 -0.71975674 0.325 0.462 1.0000000000 8215.0
# 18:      JUN 5.969008e-02 -0.30538674 0.208 0.338 1.0000000000 8184.0
# 19:     WSB1 6.454420e-02 -0.29351328 0.818 0.831 1.0000000000 8102.0
# 20:     SKIL 7.129324e-02 -0.29081163 0.805 0.846 1.0000000000 7999.0
# 21:    ARL4C 8.344200e-02 -0.31153292 0.130 0.246 1.0000000000 7848.0
# 22:    SRSF7 8.518643e-02 -0.19642878 0.987 1.000 1.0000000000 7825.0
# 23:     CNN3 8.867850e-02 -0.25849391 0.662 0.692 1.0000000000 7773.0
# 24:     TOB1 9.175405e-02  0.09888302 0.104 0.031 1.0000000000 7732.5
# 25:     UBL3 9.703169e-02 -0.24441873 0.429 0.523 1.0000000000 7662.0
res[p_val<max(res[gene%in%lead_genes]$p_val)] #5k plus petit

res_gsea[,size.regulon:=length(split(regulons$target,regulons$tf)[[pathway]]),by="pathway"]
fwrite(res_gsea,fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_padj_rank.csv.gz"))
fread(fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_padj_rank.csv.gz"))[order(pval)]


#impact on subpop
mtd<-data.table(sicbls@meta.data,keep.rownames = "bc")
mtd[,n.sample:=.N,by="sirna"]

mtd[,pct.ct:=.N/n.sample,.(sirna,cell_type2)]
ggplot(unique(mtd,by=c("sirna","cell_type2")))+
  geom_col(aes(x=cell_type2,y=pct.ct,fill=sirna),position = "dodge")

#pseudotime
  library(SeuratWrappers)
  library(monocle3)
  library(Matrix)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cbp.cds, time_bin="HSC"){
  cell_ids <- which(colData(cbp.cds)[, "cell_type2"] == time_bin )
  
  closest_vertex <-
  cbp.cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cbp.cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cbp.cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


sicbls.cds <- as.cell_data_set(sicbls)

sicbls.cds <- cluster_cells(cds = sicbls.cds, reduction_method = "UMAP")
sicbls.cds <- learn_graph(sicbls.cds, use_partition = TRUE)

sicbls.cds <- order_cells(sicbls.cds, root_pr_nodes=get_earliest_principal_node(sicbls.cds))


plot_cells(
  cds = sicbls.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points= TRUE,
  label_leaves=FALSE
)

sicbls <- AddMetaData(
  object = sicbls,
  metadata = sicbls.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime"
)

FeaturePlot(sicbls,"pseudotime",label=T,cols = c("cyan","darkorange"))
VlnPlot(sicbls,"pseudotime",group.by = "sirna")

ggplot(sicbls@meta.data)+
  geom_boxplot(aes(x=cell_type2,y=pseudotime,fill=cell_type2))

ggplot(sicbls@meta.data)+
  geom_density(aes(x=pseudotime,fill=sirna),alpha=0.5)
wilcox.test(sicbls$pseudotime[sicbls$sirna=="KLF2"],sicbls$pseudotime[sicbls$sirna=="CTRL"]) #p-value = 0.0007288

mean(sicbls$pseudotime[sicbls$sirna=="KLF2"])-mean(sicbls$pseudotime[sicbls$sirna=="CTRL"]) #+0.07

mean(sicbls$pseudotime[sicbls$sirna=="KLF2"])/mean(sicbls$pseudotime[sicbls$sirna=="CTRL"]) #1.01

quantile(sicbls$pseudotime[sicbls$sirna=="KLF2"],0.1)/quantile(sicbls$pseudotime[sicbls$sirna=="CTRL"],0.1) # q10dt= 0.23

mtd<-data.table(sicbls@meta.data)
mtd[,q10pseudo:=quantile(pseudotime,0.1),by="sirna"]
ggplot(unique(mtd,by="sirna"))+geom_col(aes(x=sirna,y=q10pseudo,fill=sirna))

saveRDS(sicbls,fp(out,"sicbls2.rds"))
