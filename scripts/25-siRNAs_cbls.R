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
sicbl#3649 samples

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
sicbls<-readRDS(fp(out,"sicbls2.rds"))

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
DimPlot(sicbls, reduction = "umap",  label = TRUE,label.size = 5) 
FeaturePlot(sicbls,features =c("MPO","LTB",'GATA2'),min.cutoff = "q10")#ok

sicbls$cell_type2<-Idents(sicbls)

DimPlot(sicbls, reduction = "umap",  label = TRUE, label.size = 4) 


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
fwrite(res,fp(out,"degs_siklf2_vs_ctrl_hsc.csv.gz"))
fwrite(resgos,fp(out,"res_gost_degs_siklf2_vs_ctrl_hsc.csv.gz"))

#enrichment for regulons KLF2
regulons<-fread("outputs/16-GRN_final/tf_target_interactions.csv")
intersect(regulons[tf=="KLF2"]$target,res[p_val_adj<0.05]$gene) # "PLCG2" 1/32 degs..
#donc nop

#gsea lik
res<-fread(fp(out,"degs_siklf2_vs_ctrl_hsc.csv.gz"))


res_gos<-gost(query = res[order(-log10(p_val)*sign(avg_log2FC))]$gene,custom_bg=res[order(-log10(p_val)*sign(avg_log2FC))]$gene,organism = "hsapiens",ordered_query = T)
resgos<-data.table(res_gos$result)

resgos[p_value<0.001&str_detect(source,"GO:BP")]
gprofiler2::gostplot(res_gos) 

resgos[p_value<1e-15][source=="TF"][1:20]
#gsea with regulons
library(fgsea)
res[,score:=sign(avg_log2FC)*rank(-p_val)]

genes_rank<-res$score
names(genes_rank)<-res$gene
genes_rank<-sort(genes_rank,decreasing = T)
head(genes_rank)
tail(genes_rank)

res_gsea<-fgsea(pathways=split(regulons[!(extended)]$target,regulons[!(extended)]$tf),
      stats=genes_rank,scoreType = "std")
res_gsea[order(pval)][1:10] #CTCF, ELF1, but also EGR1 and KLF2!
#     pathway        pval       padj   log2err         ES       NES size                                  leadingEdge
#  1:   STAT3 0.001052680 0.04762403 0.4550599 -0.4758871 -1.728147   45      DDX17,ETV6,TNRC18,STAT3,RUNX1,HOXB4,...
#  2:    KLF2 0.001322890 0.04762403 0.4550599 -0.4011175 -1.574541   82          PLCG2,ARF6,RBM39,UBC,EIF1,HLA-C,...
#  3:   CREB1 0.004052747 0.09726593 0.4070179 -0.3037233 -1.310681  315 ATP1B3,SRSF6,TMSB4X,TRA2A,AKIRIN2,YTHDC2,...
#  4:    JUND 0.005560023 0.10008041 0.4070179 -0.3885119 -1.525335   81     SRSF6,H3F3B,AKIRIN2,ARF6,G3BP2,TRA2B,...
#  5:   NR3C1 0.006972387 0.10040237 0.4070179 -0.7209525 -1.665708    7                  TRA2B,FBXW7,TCF4,CCDC6,UBL3
#  6:    JUNB 0.012517858 0.15021430 0.3807304 -0.3525416 -1.397748   94      EMP3,ANXA5,H3F3B,AKIRIN2,CD44,RAB8B,...
#  7:  HOXA10 0.020784540 0.21378384 0.3524879 -0.3490524 -1.392606   97       HOXB3,SSBP2,MEIS1,SOCS2,RUNX1,SVIL,...
#  8:   NFKB2 0.025743153 0.21995927 0.3524879 -0.3583167 -1.399739   76         PLCG2,B2M,EMP3,S100A10,BRI3,CD59,...
#  9:     FOS 0.027494908 0.21995927 0.2765006 -0.3520112 -1.382973   83       H3F3B,AKIRIN2,EEF1A1,SBDS,CD44,UBC,...
# 10:     MYC 0.031762295 0.22868852 0.2572065 -0.3746260 -1.409962   60      MYL12A,S100A6,RAN,RAP2B,ZNF706,YBX1,...
res_gsea[,size.regulon:=length(split(regulons[!(extended)]$target,regulons[!(extended)]$tf)[[pathway]]),by="pathway"]
head(res_gsea[order(padj)],30)

fwrite(res_gsea,fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_sign_p_val.csv.gz"))

res_gsea<-fread(fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_sign_p_val.csv.gz"))

#with extended
res_gsea2<-fgsea(pathways=split(regulons$target,regulons$tf),
      stats=genes_rank,scoreType = "std")
res_gsea2[order(pval)][1:10] #CTCF, ELF1, but also EGR1 and KLF2!
#     pathway        pval       padj   log2err         ES       NES size                                  leadingEdge
#  1:   STAT3 0.001052680 0.04762403 0.4550599 -0.4758871 -1.728147   45      DDX17,ETV6,TNRC18,STAT3,RUNX1,HOXB4,...
#  2:    KLF2 0.001322890 0.04762403 0.4550599 -0.4011175 -1.574541   82          PLCG2,ARF6,RBM39,UBC,EIF1,HLA-C,...
#  3:   CREB1 0.004052747 0.09726593 0.4070179 -0.3037233 -1.310681  315 ATP1B3,SRSF6,TMSB4X,TRA2A,AKIRIN2,YTHDC2,...
#  4:    JUND 0.005560023 0.10008041 0.4070179 -0.3885119 -1.525335   81     SRSF6,H3F3B,AKIRIN2,ARF6,G3BP2,TRA2B,...
#  5:   NR3C1 0.006972387 0.10040237 0.4070179 -0.7209525 -1.665708    7                  TRA2B,FBXW7,TCF4,CCDC6,UBL3
#  6:    JUNB 0.012517858 0.15021430 0.3807304 -0.3525416 -1.397748   94      EMP3,ANXA5,H3F3B,AKIRIN2,CD44,RAB8B,...
#  7:  HOXA10 0.020784540 0.21378384 0.3524879 -0.3490524 -1.392606   97       HOXB3,SSBP2,MEIS1,SOCS2,RUNX1,SVIL,...
#  8:   NFKB2 0.025743153 0.21995927 0.3524879 -0.3583167 -1.399739   76         PLCG2,B2M,EMP3,S100A10,BRI3,CD59,...
#  9:     FOS 0.027494908 0.21995927 0.2765006 -0.3520112 -1.382973   83       H3F3B,AKIRIN2,EEF1A1,SBDS,CD44,UBC,...
# 10:     MYC 0.031762295 0.22868852 0.2572065 -0.3746260 -1.409962   60      MYL12A,S100A6,RAN,RAP2B,ZNF706,YBX1,...
res_gsea[,size.regulon:=length(split(regulons[!(extended)]$target,regulons[!(extended)]$tf)[[pathway]]),by="pathway"]
head(res_gsea[order(padj)],30)

fwrite(res_gsea,fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_sign_p_val.csv.gz"))


#fig
res_gsea<-fread(fp(out,"res_gsea_degs_klf2_hsc_regulons_filtered_sign_p_val.csv.gz"))

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

ggplot(res_gsea,aes(x=NES,y=-log10(pval),col=padj<0.1))+
  geom_point(aes(size=size.regulon))+
  scale_color_manual(values = c("grey","red")) +
  geom_label_repel(aes(label=ifelse(padj<0.1,pathway,"")))+theme_minimal()


#impact on subpop
table(sicbls$cell_type2,sicbls$sirna)
  #               CTRL KLF2
  # HSC             65   77
  # MPP            162  178
  # MPP-Ery        312  449
  # EMP            247  319
  # ErP/GMP-cycle   99  142
  # GMP             55   82
  # LMPP           534  783
  # CLP             59   86
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


#figures

table(sicbls$sirna)
# CTRL KLF2 
# 1533 2116 

round(table(sicbls$cell_type2)/ncol(sicbls) * 100)
          # HSC           MPP       MPP-Ery           EMP ErP/GMP-cycle 
          #   4             9            21            16             7 
          # GMP          LMPP           CLP 
          #   4            36             4 

#vs unincubated cells 
mtdo<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
round(table(mtdo$lineage_hmap)/nrow(mtd) * 100)


#subpop shift

mtdsi<-data.table(sicbls@meta.data)
mtdsi[,pct.hscmpp:=sum(cell_type2%in%c("HSC","MPP"))/.N,by="sirna"]

ggplot(unique(mtdsi,by="sirna"))+geom_col(aes(x=sirna,y=pct.hscmpp,fill=sirna))
chisq.test(c(sum(mtdsi$sirna=="KLF2"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="KLF2"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   p = c(sum(mtdsi$sirna=="CTRL"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="CTRL"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   rescale.p = T)$p.value #0.0035

chisq.test(c(sum(mtdsi$sirna=="KLF2"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="KLF2"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   y = c(sum(mtdsi$sirna=="CTRL"&mtdsi$cell_type2%in%c("HSC","MPP")),sum(mtdsi$sirna=="CTRL"&!(mtdsi$cell_type2%in%c("HSC","MPP")))),
                                   )$p.value #0.0035

#degs regulons KLF2
regulons[tf=="KLF2"&!(extended)]
DefaultAssay(sicbls)<-"SCT"

res_klf2<-data.table(FindMarkers(sicbls,subset.ident = "HSC",
            ident.1 = "KLF2",ident.2 = "CTRL",group.by = "sirna",
            features =regulons[tf=="KLF2"&!(extended)]$target,
            logfc.threshold = 0,test.use = "t"),
            keep.rownames = "gene")

res_klf2[,padj:=p.adjust(p_val,method = 'BH')]

res_klf2[padj<0.05]
#     gene        p_val avg_log2FC pct.1 pct.2   p_val_adj         padj
# 1: PLCG2 4.280812e-07 -0.5510063 1.000 1.000 0.007538939 3.510266e-05
# 2:  ARF6 9.416765e-05 -0.5735784 0.688 0.877 1.000000000 3.860874e-03
# 3: RBM39 8.927806e-04 -0.3631184 1.000 1.000 1.000000000 2.112964e-02
# 4:   UBC 1.030714e-03 -0.5446809 1.000 1.000 1.000000000 2.112964e-02
res_klf2[p_val<0.05]
#      gene        p_val avg_log2FC pct.1 pct.2   p_val_adj         padj
#  1:  PLCG2 4.280812e-07 -0.5510063 1.000 1.000 0.007538939 3.510266e-05
#  2:   ARF6 9.416765e-05 -0.5735784 0.688 0.877 1.000000000 3.860874e-03
#  3:  RBM39 8.927806e-04 -0.3631184 1.000 1.000 1.000000000 2.112964e-02
#  4:    UBC 1.030714e-03 -0.5446809 1.000 1.000 1.000000000 2.112964e-02
#  5:  HLA-E 3.503152e-03 -0.3153022 1.000 1.000 1.000000000 5.745170e-02
#  6:   EIF1 6.374270e-03 -0.2945259 1.000 1.000 1.000000000 8.011478e-02
#  7:   RHOC 7.533594e-03 -0.4774781 0.545 0.692 1.000000000 8.011478e-02
#  8: TIPARP 7.816077e-03 -0.3791438 0.831 0.938 1.000000000 8.011478e-02
#  9:  HLA-C 1.382889e-02 -0.2935555 1.000 1.000 1.000000000 1.259966e-01
# 10:  ITM2B 3.388506e-02 -0.1998679 1.000 1.000 1.000000000 2.494340e-01
# 11:  KLF10 3.677988e-02 -0.3619883 0.234 0.354 1.000000000 2.494340e-01
# 12:   LMNA 4.014808e-02 -0.7197567 0.325 0.462 1.000000000 2.494340e-01
# 13:  SRSF7 4.150240e-02 -0.1964288 0.987 1.000 1.000000000 2.494340e-01
# 14:  SRSF5 4.468353e-02 -0.2052173 1.000 1.000 1.000000000 2.494340e-01
# 15:  RAB21 4.562816e-02 -0.3711107 0.701 0.785 1.000000000 2.494340e-01
#not work

#res unsupervised analysis
res[p_val_adj<0.05]

res_gsea<-fread("outputs/25-siRNAs_KLF2_cbl/res_gsea_degs_klf2_hsc_regulons_filtered_sign_p_val.csv.gz")
res_gsea[order(pval)]
res_gsea[,top10:=rank(pval)<11]
ggplot(res_gsea,aes(x=NES,y=-log10(pval),col=padj<0.1))+
  geom_point(aes(size=size.regulon))+
  scale_color_manual(values = c("grey","red")) +
  geom_label_repel(aes(label=ifelse(top10,pathway,"")))


#aucell regulons
  library(AUCell)
  library(SCENIC)
  source("scripts/utils/scenic_utils.R")

cells_rankings <- AUCell_buildRankings(as.matrix(sicbls@assays$RNA@counts),nCores = 20)
regulons_list<-split(regulons[!(extended)]$target,regulons[!(extended)]$tf)
cells_AUC <- AUCell_calcAUC(regulons_list, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
auc_mat<-getAUC(cells_AUC)
sicbls@assays[["TF_AUC"]] <- CreateAssayObject(auc_mat)
DefaultAssay(sicbls)<-"TF_AUC"

saveRDS(sicbls,fp(out,"sicbls2_with_regulons_activity.rds"))

FeaturePlot(sicbls,"KLF2",min.cutoff = "q10")

VlnPlot(sicbls,"KLF2",group.by = "cell_type2")+legend(x="cell_type",y="Regulon activity")

VlnPlot(sicbls,"KLF2",group.by = "cell_type2",split.by = "sirna")


res<-Reduce(rbind,lapply(levels(sicbls), function(ct){
  return(data.table(FindMarkers(sicbls,subset.ident = ct,logfc.threshold = 0,
                                group.by = "sirna",ident.1 = "KLF2"),keep.rownames = "regulon")[,cell_type2:=ct])
}))

res[p_val_adj<0.05]
res[p_val<0.05&cell_type2=="HSC"]

