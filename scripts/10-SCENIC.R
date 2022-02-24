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

tf_diff<-data.table(FindMarkers(cbps,assay="TF_AUC",
                 logfc.threshold = 0,
                 min.pct = 0,
                 mean.fxn = rowMeans,
                 fc.name = "avg_diff",
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
                                  min.pct = 0,
                                  mean.fxn = rowMeans,
                                  fc.name = "avg_diff",
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

head(tf_diff_merge[order(p_val_adj)][!str_detect(regulon,"e$")],20)
tf_diff_merge[order(p_val_adj)][!str_detect(regulon,"e$")][lineage=="HSC"&p_val_adj<0.001&abs(avg_diff)>0.03]
#    regulon        p_val    avg_diff pct.1 pct.2    p_val_adj lineage   hto
# 1:    CHD1 7.297385e-76  0.04498438 1.000 0.999 1.824346e-73     HSC FALSE
# 2:    ELF1 3.532168e-59  0.03256795 0.999 0.999 8.830419e-57     HSC FALSE
# 3:  ARID5A 4.279509e-57 -0.05573158 1.000 1.000 1.069877e-54     HSC  TRUE
# 4:   KLF10 4.323276e-56 -0.03207762 1.000 1.000 1.080819e-53     HSC  TRUE
# 5:    EGR1 6.308418e-54 -0.05621167 1.000 1.000 1.577104e-51     HSC  TRUE
# 6:    FOSB 3.068542e-50 -0.05625283 1.000 1.000 7.671355e-48     HSC  TRUE
# 7:    KLF4 2.574408e-48 -0.04466093 1.000 1.000 6.436020e-46     HSC  TRUE
# 8:    KLF2 2.044410e-46 -0.03816241 1.000 1.000 5.111024e-44     HSC  TRUE
# 9:     JUN 3.100966e-39 -0.03612825 1.000 1.000 7.752415e-37     HSC  TRUE

tf_diff_merge[order(p_val_adj)][!str_detect(regulon,"e$")][lineage=="MPP/LMPP"&p_val_adj<0.001&abs(avg_diff)>0.03]
#    regulon         p_val   avg_diff pct.1 pct.2     p_val_adj  lineage   hto
# 1:    JUND 5.139971e-198 0.04005085 1.000 1.000 1.284993e-195 MPP/LMPP  TRUE
# 2:    CHD1 3.941817e-112 0.03857185 0.999 0.999 9.854543e-110 MPP/LMPP FALSE
# 3:  BCLAF1  3.263674e-63 0.03607896 0.999 0.997  8.159185e-61 MPP/LMPP  TRUE
# 4:   FOSL2  1.576938e-56 0.04148842 0.999 0.998  3.942344e-54 MPP/LMPP  TRUE
# 5:    CHD1  1.538326e-41 0.03761871 0.999 1.000  3.845816e-39 MPP/LMPP  TRUE
tf_diff_merge[order(p_val_adj)][!str_detect(regulon,"e$")][lineage=="Myeloid"&p_val_adj<0.001&abs(avg_diff)>0.03]
#    regulon        p_val    avg_diff pct.1 pct.2    p_val_adj lineage  hto
# 1:    KLF4 4.352730e-12 -0.05246811 1.000 1.000 1.088183e-09 Myeloid TRUE
# 2:   RUNX2 7.882590e-11  0.03216305 1.000 1.000 1.970647e-08 Myeloid TRUE
# 3:   KLF10 9.408117e-11 -0.03245992 1.000 1.000 2.352029e-08 Myeloid TRUE
# 4: ZSCAN18 6.706685e-09 -0.03396890 0.775 0.825 1.676671e-06 Myeloid TRUE
# 5:    KLF2 8.625406e-09 -0.03770843 1.000 1.000 2.156352e-06 Myeloid TRUE
# 6:  BCL11A 3.318704e-07  0.03600228 0.997 0.982 8.296760e-05 Myeloid TRUE

tf_diff_merge[order(p_val_adj)][!str_detect(regulon,"e$")][lineage=="Lymphoid"&p_val_adj<0.001&abs(avg_diff)>0.03]
#     regulon        p_val    avg_diff pct.1 pct.2    p_val_adj  lineage  hto
#  1:   HOXA6 1.470871e-18 -0.04672259 1.000 1.000 3.677179e-16 Lymphoid TRUE
#  2:    KLF4 2.826992e-17 -0.03786592 1.000 1.000 7.067480e-15 Lymphoid TRUE
#  3:    JUNB 1.321600e-16 -0.03412035 1.000 1.000 3.304001e-14 Lymphoid TRUE
#  4:     FOS 3.102292e-16 -0.03720474 1.000 1.000 7.755730e-14 Lymphoid TRUE
#  5:    EGR1 9.248277e-16 -0.06765369 1.000 1.000 2.312069e-13 Lymphoid TRUE
#  6:    FOSB 3.391711e-15 -0.06710951 1.000 1.000 8.479278e-13 Lymphoid TRUE
#  7:     JUN 7.510019e-15 -0.04062413 1.000 1.000 1.877505e-12 Lymphoid TRUE
#  8:  ARID5A 1.126868e-14 -0.07054619 1.000 1.000 2.817171e-12 Lymphoid TRUE
#  9:    KLF2 1.435747e-13 -0.03344990 1.000 1.000 3.589368e-11 Lymphoid TRUE
# 10:   STAT3 2.531644e-13 -0.03299205 1.000 1.000 6.329110e-11 Lymphoid TRUE
# 11: SMARCA4 4.082356e-11  0.04631834 1.000 1.000 1.020589e-08 Lymphoid TRUE
# 12:   FOSL2 5.794451e-11 -0.04357780 0.991 0.997 1.448613e-08 Lymphoid TRUE
# 13:   CEBPB 6.476611e-11 -0.03232078 1.000 1.000 1.619153e-08 Lymphoid TRUE
# 14:   DDIT3 1.818829e-09 -0.03094514 0.994 0.994 4.547073e-07 Lymphoid TRUE
# 15:  BCL11A 2.128758e-09  0.03303592 1.000 1.000 5.321894e-07 Lymphoid TRUE
# 16:    IRF1 3.982996e-09 -0.04357225 0.947 0.948 9.957489e-07 Lymphoid TRUE
# 17:    TCF3 7.919680e-09  0.07655157 1.000 1.000 1.979920e-06 Lymphoid TRUE
# 18:     MAF 9.841950e-09 -0.03700230 0.840 0.876 2.460487e-06 Lymphoid TRUE
tf_diff_merge[order(p_val_adj)][!str_detect(regulon,"e$")][lineage=="Erythro-Mas"&p_val_adj<0.001&abs(avg_diff)>0.03]
#    regulon        p_val    avg_diff pct.1 pct.2    p_val_adj     lineage   hto
# 1:    ELF1 9.358429e-51  0.03698465 0.999 0.998 2.339607e-48 Erythro-Mas FALSE
# 2:  STAT5A 3.997241e-28 -0.03473640 0.923 0.943 9.993103e-26 Erythro-Mas FALSE
# 3:  BCLAF1 3.443776e-07  0.03188338 1.000 0.994 8.609440e-05 Erythro-Mas  TRUE
# 4:    CHD1 9.396316e-07  0.03158779 1.000 0.991 2.349079e-04 Erythro-Mas  TRUE

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





