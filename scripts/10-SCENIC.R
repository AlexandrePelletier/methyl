#scenic on cbps0-8 
out<-"outputs/10-SCENIC"
dir.create(out)
library(Seurat)
#need first batch correct the counts with seuratv3
#run 10A-



#run do_scenic_on_cbps0-8_clean.r with working directory = out #[to update]

#get regulon activity by cells 
#run 10B-get_regul_activity


#added integrated and SCENIC assa cbps_map_on_hmap
cbps<-readRDS(fp(out,"cbps_with_regulons_activity.rds"))
saveRDS(cbps@assays$TF_AUC,fp(out,"TF_AUC_assay.rds"))
VlnPlot(cbps, c("STAT3","GATA1","SPI1"),group.by="lineage",pt.size = 0) 

DefaultAssay(cbps)<-"TF_AUC"

#valid clustering with scenic
#cluster with tf mat [to update]
cbps<-FindVariableFeatures(cbps)
cbps<-ScaleData(cbps)
cbps<-RunPCA(cbps,reduction.name = "scenic.pca")
cbps<-RunUMAP(cbps,reduction.name = "scenic.pca",dims=1:30)
cbps <- FindNeighbors(object = cbps,reduction = "scenic.pca", dims = 1:30,graph.name = "tf_graph")

cbps<- FindClusters(cbps,resolution = 0.6,
                               algorithm = 1,graph.name = "tf_graph") 
p1<-DimPlot(cbps,label = T)
p2<-DimPlot(cbps,group.by="cell_type_hmap",label = T)
p1+p2

#if doesnt work, binarize version 

#find Markers by cell_type
Idents(cbps)<-"cell_type_hmap"
DefaultAssay(cbps)<-"TF_AUC"
tf_markers<-FindAllMarkers(cbps,logfc.threshold = 0)
tf_markers<-data.table(tf_markers)[,cell_type:=cluster][,regulon:=gene][,-c('cluster','gene')]
fwrite(tf_markers,fp(out,"tf_markers_cell_type.csv.gz"),sep=";")

saveRDS(cbps,"../singlecell/outputs/cbps0_8.rds")

#[end to do update]


#tf activity stim vs not  by lineage
Idents(cbps)<-"lineage_hmap"
i<-0
for(lin in levels(cbps)){
print(lin)

tf_diff<-data.table(FindMarkers(cbps,assay="TF_AUC",logfc.threshold=0,
                                subset.ident = lin,group.by = "hto",
                     ident.1 = TRUE,ident.2 = FALSE),keep.rownames = "regulon")
i<-i+1

tf_diff[,lineage:=lin]

if(i==1){
  tf_diff_merge<-copy(tf_diff)
}else{
  tf_diff_merge<-rbind(tf_diff_merge,tf_diff,fill=T)
}
}

fwrite(tf_diff_merge,fp(out,"regulon_activity_HTO_vs_not_by_lineage.csv.gz"),sep=";")



#tf activity LGA vs ctrl by lineage and hto
Idents(cbps)<-"lineage_hmap"


i<-0
for(hto_ in c(TRUE,FALSE)){
  for(lin in levels(cbps)){
  print(lin)
  
  tf_diff<-data.table(FindMarkers(cbps,assay="TF_AUC",logfc.threshold=0,
                                  subset.ident = lin,group.by = "group_hto",
                       ident.1 = paste0("lga",hto_),ident.2 = paste0("ctrl",hto_)),keep.rownames = "regulon")
  i<-i+1
  
  tf_diff[,lineage:=lin]
  tf_diff[,hto:=hto_]

  if(i==1){
    tf_diff_merge<-copy(tf_diff)
  }else{
    tf_diff_merge<-rbind(tf_diff_merge,tf_diff,fill=T)
    }
  

  
  }
}


fwrite(tf_diff_merge,fp(out,"regulon_activity_lga_vs_ctrl_HTO_by_lineage.csv.gz"),sep=";")

#compare auc by group [to update]
auc_score<-merge(data.table(cbps@meta.data[,colnames(cbps@meta.data)!="bc"],keep.rownames = 'bc'),melt.data.table(data.table(as.matrix(cbps@assays$SCENIC@counts),keep.rownames = "regulon"),variable.name ="bc" ,value.name ="score" ))
auc_score[,score_scaled:=scale(score),by=.(regulon)]
auc_score[,avg_score:=mean(score_scaled),by=.(regulon,group,hto,lineage_hmap)]
auc_score[,score_diff_lga.ctrl:=unique(avg_score[group=="lga"])-unique(avg_score[group=="ctrl"]),by=.(regulon,lineage_hmap,hto)]
auc_score_diff<-unique(auc_score[,.(regulon, lineage_hmap, hto, score_diff_lga.ctrl)],by=c("regulon","lineage_hmap","hto"))[!is.na(regulon)]
auc_score_diff[lineage_hmap=="HSC"&order(score_diff_lga.ctrl)]$regulon


tfs_of_interest<-c("STAT3e","STAT3","EGR1","JUN","FOSB","ARID5A")
tf_int_dt<-data.table(t(as.matrix(cbps@assays$SCENIC@data[tfs_of_interest,])),keep.rownames = "cell")
tf_int_dt<-melt(tf_int_dt,id.vars = "cell",variable.name ="regulon",value.name = "activity" )
for(tf in tfs_of_interest){
  print(tf)
  print(ggplot(tf_int_dt[regulon==tf])+geom_density(aes(x=activity)))
  tf_int_dt[regulon==tf,activ.thr:=as.numeric(readline("threshold: "))]
  
}

tf_int_dt[,activated:=activity>activ.thr]
mtd<-data.table(cbps@meta.data,keep.rownames = "cell")
tf_int_dt_mtd<-merge(tf_int_dt,mtd[,.(cell,sample,sample_hto,group,batch,hto,lineage_hmap,cell_type_hmap)],by="cell")
tf_int_dt_mtd[,pct.activ:=sum(activated)/.N,,by=.(sample_hto,lineage_hmap,regulon)]

fwrite(tf_int_dt_mtd,fp(out,"activated_tf_of_interest_by_cells.csv.gz"))





