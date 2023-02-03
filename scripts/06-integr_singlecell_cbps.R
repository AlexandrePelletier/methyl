#integr CBPs datasets thanks to hematomap
out<-"outputs/06-integr_singlecell_cbps"
dir.create(out)
source("scripts/utils/new_utils.R")
source("../singlecell/scripts/utils/HTO_utils.R")
library(Seurat)

####QC filtering and Demultiplexing  data ####
# remove low quality/outlyers cells (doublet with nGene/RNA hi/lo or percent.mt Hi)
# threshold : if not 2 subpop :  median +/- 4* median absolute deviation(mad), else cutoff based on distribution 
#cbp4####
sample<-"cbp4"
umis<- Read10X("~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/filtered_feature_bc_matrix/")$`Gene Expression`

cbp4_all <- CreateSeuratObject(counts = umis,project = sample)
cbp4_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp4_all, pattern = "^MT-")

VlnPlot(object = cbp4_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
p1<-VlnPlot(object = cbp4_all, features ="percent.mt")+geom_hline(yintercept = median(cbp4_all$percent.mt)+ 4*mad(cbp4_all$percent.mt) )

p2<-VlnPlot(object = cbp4_all, features ="nCount_RNA")+geom_hline(yintercept = 60000)

p3<-VlnPlot(object = cbp4_all, features ="nFeature_RNA")+geom_hline(yintercept = 7000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp4_QC_cells_metrics.png"))

cbp4_qc<-subset(cbp4_all,percent.mt<median(cbp4_all$percent.mt)+ 4*mad(cbp4_all$percent.mt)&nCount_RNA<60000&nFeature_RNA<7000)

#reassign samples
cbp4.htos<-as.matrix(Read10X("~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/filtered_feature_bc_matrix/")$`Antibody Capture`)
rownames(cbp4.htos)<-c("ctrlM555",
                          "ctrlM518",
                          "ctrlM537",
                          "lgaF551",
                          "lgaF543")

cbp4_qc[["HTO"]] <- CreateAssayObject(counts = cbp4.htos[,colnames(cbp4_qc)])

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cbp4_qc <- NormalizeData(cbp4_qc, assay = "HTO", normalization.method = "CLR")
cbp4_qc <- HTODemux(cbp4_qc, assay = "HTO",positive.quantile = 0.95)
table(cbp4_qc$HTO_classification.global)
 # Doublet Negative  Singlet 
 #    1017     2373     2487 


#sex based recovery
cbp4_qc<-checkHTOSex(cbp4_qc,gene_male="RPS4Y1",gene_female="XIST")
# calculating pct of real singlet male/female cells expressing sex marker ( save in 'misc' of 'HTO' assay):
# for male :   82 %  express the male gene
# for female : 97% express the female gene

#as ++ doublet/singlet, and only 82% male cells express RPS4y1, split in 2 based hi / lo RNA count
VlnPlot(object = cbp4_all, features ="nCount_RNA")+geom_hline(yintercept = 8000)
cbp4_qc_lo<-subset(cbp4_qc,nCount_RNA<8000)
cbp4_qc_lo <- NormalizeData(cbp4_qc_lo, assay = "HTO", normalization.method = "CLR")
cbp4_qc_lo <- HTODemux(cbp4_qc_lo, assay = "HTO",positive.quantile = 0.9999)
table(cbp4_qc_lo$HTO_classification.global)
 # Doublet Negative  Singlet 
 #     323      605      739 

cbp4_qc_lo<-checkHTOSex(cbp4_qc_lo,gene_male="RPS4Y1",gene_female="XIST")
# for male :  52 % express the male gene
# for female :  85 % express the female gene
cbp4_qc_lo<-sexBasedHTOAssign(cbp4_qc_lo)
cbp4_qc_lo_s<-subset(cbp4_qc_lo,new.HTO_classif.global=="Singlet")

cbp4_qc_hi<-subset(cbp4_qc,nCount_RNA>=8000)
cbp4_qc_hi <- NormalizeData(cbp4_qc_hi, assay = "HTO", normalization.method = "CLR")
cbp4_qc_hi <- HTODemux(cbp4_qc_hi, assay = "HTO",positive.quantile = 0.95)
table(cbp4_qc_hi$HTO_classification.global)
# Doublet Negative  Singlet 
#      880     1511     1819 

cbp4_qc_hi<-checkHTOSex(cbp4_qc_hi,gene_male="RPS4Y1",gene_female="XIST")
# for male :  100 % express the male gene
# for female :  100 % express the female gene
cbp4_qc_hi<-sexBasedHTOAssign(cbp4_qc_hi)
# Bad_HTO_assign        Doublet       Negative        Singlet 
#            112           1201            446           2451 
cbp4_qc_hi_s<-subset(cbp4_qc_hi,new.HTO_classif.global=="Singlet")

#merge the 2

cbp4_qc_s<-merge(cbp4_qc_hi_s,cbp4_qc_lo_s)

cbp4_qc_s$sample<-cbp4_qc_s$new.ID
table(cbp4_qc_s$sample)
# ctrlM518 ctrlM537 ctrlM555  lgaF543  lgaF551 
#      678      470      619     1013      745 
saveRDS(cbp4_qc_s,fp(out,"cbp4.rds"))


#cbp2####
sample<-"cbp2"
umis<- Read10X("~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp2b_tri/single_cell_barcode_539_HTO_cbp2b/outs/filtered_feature_bc_matrix/")$`Gene Expression`

cbp2_all <- CreateSeuratObject(counts = umis,project = sample)
cbp2_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp2_all, pattern = "^MT-")

VlnPlot(object = cbp2_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
p1<-VlnPlot(object = cbp2_all, features ="percent.mt")+geom_hline(yintercept = median(cbp2_all$percent.mt)+ 4*mad(cbp2_all$percent.mt) )

p2<-VlnPlot(object = cbp2_all, features ="nCount_RNA")+geom_hline(yintercept = 60000)

p3<-VlnPlot(object = cbp2_all, features ="nFeature_RNA")+geom_hline(yintercept = 6500 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp2_QC_cells_metrics.png"))

cbp2_qc<-subset(cbp2_all,percent.mt<median(cbp2_all$percent.mt)+ 4*mad(cbp2_all$percent.mt)&nCount_RNA<60000&nFeature_RNA<6500)

#reassign samples

cbp2.htos<-as.matrix(Read10X("~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp2b_tri/single_cell_barcode_539_HTO_cbp2b/outs/filtered_feature_bc_matrix/")$`Antibody Capture`)
rownames(cbp2.htos)<-c("ctrlF528",
                           "ctrlM539",
                           "iugrM553",
                           "iugrM558",
                           "lgaM549",
                           "lgaF532")

cbp2_qc[["HTO"]] <- CreateAssayObject(counts = cbp2.htos[,colnames(cbp2_qc)])

cbp2_qc <- NormalizeData(cbp2_qc, assay = "HTO", normalization.method = "CLR")
cbp2_qc <- HTODemux(cbp2_qc, assay = "HTO",positive.quantile = 0.9999)
table(cbp2_qc$HTO_classification.global)
 # Doublet Negative  Singlet 
 #    9536     1143       98 

#as ++ doublet split in 2 based hi / lo RNA count
VlnPlot(object = cbp2_qc, features ="nCount_RNA",group.by = "orig.ident")+geom_hline(yintercept = 6000)
cbp2_qc_lo<-subset(cbp2_qc,nCount_RNA<6000)
cbp2_qc_lo <- NormalizeData(cbp2_qc_lo, assay = "HTO", normalization.method = "CLR")
cbp2_qc_lo <- HTODemux(cbp2_qc_lo, assay = "HTO",positive.quantile = 0.95)
table(cbp2_qc_lo$HTO_classification.global)
 # Doublet Negative  Singlet 
 #    1749     1002       49 

cbp2_qc_lo<-checkHTOSex(cbp2_qc_lo,gene_male="RPS4Y1",gene_female="XIST")
# for male :  52 % express the male gene
# for female :  36 % express the female gene

cbp2_qc_hi<-subset(cbp2_qc,nCount_RNA>=6000)
cbp2_qc_hi <- NormalizeData(cbp2_qc_hi, assay = "HTO", normalization.method = "CLR")
cbp2_qc_hi <- HTODemux(cbp2_qc_hi, assay = "HTO",positive.quantile = 0.99)
table(cbp2_qc_hi$HTO_classification.global)
 # Doublet Negative  Singlet 
 #    1323     1844     4810 

cbp2_qc_hi<-checkHTOSex(cbp2_qc_hi,gene_male="RPS4Y1",gene_female="XIST")
# for male :  99 % express the male gene
# for female :  98 % express the female gene
cbp2_qc_hi<-sexBasedHTOAssign(cbp2_qc_hi)
# Bad_HTO_assign        Doublet       Negative        Singlet 
#             84           1596            503           5794 
#merge the 2

cbp2_qc<-merge(cbp2_qc_hi,cbp2_qc_lo)

#demultiplex on SNP info
mat_gt<-fread("../lineage_tracing/outputs/cbp2/25pct_det_and_variables_snps_barcodes_genotype_matrices_imputed.tsv")
mat_gt<-as.matrix(data.frame(mat_gt,row.names = "snp"))
dim(mat_gt)
colnames(mat_gt)<-str_replace(colnames(mat_gt),"\\.","-")

sum(colnames(cbp2_qc)%in%colnames(mat_gt)) #8292
cbp2_qc_snp<-cbp2_qc[,colnames(cbp2_qc)%in%colnames(mat_gt)]
cbp2_qc_snp[["SNP"]]<-CreateAssayObject(data=mat_gt[,colnames(cbp2_qc_snp)])
DefaultAssay(cbp2_qc_snp)<-"SNP"
cbp2_qc_snp<-FindVariableFeatures(cbp2_qc_snp)
cbp2_qc_snp<-ScaleData(cbp2_qc_snp)
cbp2_qc_snp<-RunPCA(cbp2_qc_snp)
cbp2_qc_snp<-RunUMAP(cbp2_qc_snp,dims = 1:10)
cbp2_qc_snp<-FindNeighbors(cbp2_qc_snp,dims = 1:10)
cbp2_qc_snp<-FindClusters(cbp2_qc_snp,resolution = 0.5)

DimPlot(cbp2_qc_snp,label=T,group.by = c("seurat_clusters","new.ID"))
new.idents<-c("lgaM549",
  "iugrM558",
  "ctrlM539",
  "lgaF532",
  "iugrM553",
  "lgaM549",
    "ctrlF528",
  "Doublet",
  "ctrlM539",
  "iugrM558",
  "Doublet",
  "Doublet",
  "Doublet",
  "Doublet"
  )
names(new.idents)<-levels(cbp2_qc_snp)
cbp2_qc_snp<-RenameIdents(cbp2_qc_snp,new.idents)
DimPlot(cbp2_qc_snp,label = T)
cbp2_qc_snp[["snp.ID"]]<-Idents(cbp2_qc_snp)

cbp2_qc_s<-subset(cbp2_qc_snp,snp.ID!="Doublet")

cbp2_qc_s$sample<-cbp2_qc_s$snp.ID
table(cbp2_qc_s$sample)
 # lgaM549 iugrM558 ctrlM539  lgaF532 iugrM553 ctrlF528 
 #    2459     1381     1358      929      768      692 
saveRDS(cbp2_qc_s,fp(out,"cbp2.rds"))


#cbp3####
sample<-"cbp3"
umis<- Read10X("../singlecell/datasets/cbp3/filtered_feature_bc_matrix/")

cbp3_all <- CreateSeuratObject(counts = umis,project = sample)
cbp3_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp3_all, pattern = "^MT-")

VlnPlot(object = cbp3_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
p1<-VlnPlot(object = cbp3_all, features ="percent.mt")+geom_hline(yintercept = median(cbp3_all$percent.mt)+ 4*mad(cbp3_all$percent.mt) )

p2<-VlnPlot(object = cbp3_all, features ="nCount_RNA")+geom_hline(yintercept = 60000)

p3<-VlnPlot(object = cbp3_all, features ="nFeature_RNA")+geom_hline(yintercept = 6500 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp3_QC_cells_metrics.png"))

cbp3_qc<-subset(cbp3_all,percent.mt<median(cbp3_all$percent.mt)+ 4*mad(cbp3_all$percent.mt)&nCount_RNA<60000&nFeature_RNA<6500)

#reassign samples
cbp3.htos<-as.matrix(Read10X("../singlecell/datasets/cbp3/HTO_CBP3/umi_count/",gene.column = 1))[1:3,]
rownames(cbp3.htos)<-c("ctrlF523",
                          "iugrF524",
                          "lgaF559")
colnames(cbp3.htos)<-paste0(colnames(cbp3.htos),"-1")
sum(colnames(cbp3_qc)%in%colnames(cbp3.htos)) #2570/2590
common_bc<-intersect(colnames(cbp3_qc),colnames(cbp3.htos))
cbp3_qc<-cbp3_qc[,common_bc]
cbp3_qc[["HTO"]] <- CreateAssayObject(counts = cbp3.htos[,common_bc])

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cbp3_qc <- NormalizeData(cbp3_qc, assay = "HTO", normalization.method = "CLR")
cbp3_qc <- HTODemux(cbp3_qc, assay = "HTO",positive.quantile = 0.98)
table(cbp3_qc$HTO_classification.global)
#++soublet need split in 2

VlnPlot(object = cbp3_all, features ="nFeature_RNA")+geom_hline(yintercept = 800 )

cbp3_qc_lo<-subset(cbp3_qc,nFeature_RNA<800)
cbp3_qc_lo <- NormalizeData(cbp3_qc_lo, assay = "HTO", normalization.method = "CLR")
cbp3_qc_lo <- HTODemux(cbp3_qc_lo, assay = "HTO",positive.quantile = 0.95)
table(cbp3_qc_lo$HTO_classification.global)
 # Doublet Negative  Singlet 
 #     321      411      108 
cbp3_qc_hi<-subset(cbp3_qc,nFeature_RNA>=800)
cbp3_qc_hi <- NormalizeData(cbp3_qc_hi, assay = "HTO", normalization.method = "CLR")
cbp3_qc_hi <- HTODemux(cbp3_qc_hi, assay = "HTO",positive.quantile = 0.98)
table(cbp3_qc_hi$HTO_classification.global)
 # Doublet Negative  Singlet 
 #     278      723      729 
cbp3_qc<-merge(cbp3_qc_hi,cbp3_qc_lo)
cbp3_qc_s<-subset(cbp3_qc,HTO_classification.global=="Singlet")

cbp3_qc_s$sample<-cbp3_qc_s$hash.ID
table(cbp3_qc_s$sample)
# ctrlF523 iugrF524  lgaF559 
#      206      368      263 
saveRDS(cbp3_qc_s,fp(out,"cbp3.rds"))


#cbp8####
sample<-"cbp8"
umis<- Read10X("~/RUN/Run_554_10x_standard/Output/cellranger_count_hto/single_cell_barcode_run_554_cbp8/outs/filtered_feature_bc_matrix/")$`Gene Expression`

cbp8_all <- CreateSeuratObject(counts = umis,project = sample)
cbp8_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp8_all, pattern = "^MT-")

VlnPlot(object = cbp8_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
p1<-VlnPlot(object = cbp8_all, features ="percent.mt")+geom_hline(yintercept = median(cbp8_all$percent.mt)+ 4*mad(cbp8_all$percent.mt) )

p2<-VlnPlot(object = cbp8_all, features ="nCount_RNA")+geom_hline(yintercept = 40000)

p3<-VlnPlot(object = cbp8_all, features ="nFeature_RNA")+geom_hline(yintercept = 6000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp8_QC_cells_metrics.png"))

cbp8_qc<-subset(cbp8_all,percent.mt<median(cbp8_all$percent.mt)+ 4*mad(cbp8_all$percent.mt)&nCount_RNA<40000&nFeature_RNA<6000)

#reassign samples
cbp8.htos<-as.matrix(Read10X("~/RUN/Run_554_10x_standard/Output/cellranger_count_hto/single_cell_barcode_run_554_cbp8/outs/filtered_feature_bc_matrix/")$`Antibody Capture`)
rownames(cbp8.htos)<-c("ctrlM503",
                       "ctrlM530",
                       "lgaM556",
                       "lgaM496")

cbp8_qc[["HTO"]] <- CreateAssayObject(counts = cbp8.htos[,colnames(cbp8_qc)])

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
cbp8_qc <- NormalizeData(cbp8_qc, assay = "HTO", normalization.method = "CLR")
cbp8_qc <- HTODemux(cbp8_qc, assay = "HTO",positive.quantile = 0.97)

table(cbp8_qc$HTO_classification.global)
 # Doublet Negative  Singlet 
 #     843     1367     3652 


cbp8_qc_s<-subset(cbp8_qc,HTO_classification.global=="Singlet")

cbp8_qc_s$sample<-cbp8_qc_s$hash.ID
saveRDS(cbp8_qc_s,fp(out,"cbp8.rds"))


#cbp0_ctrl####
sample<-"cbp0_ctrl"
umis<- as.matrix(data.frame(fread("../singlecell/datasets/cbp0/CBP547.csv"),row.names=1))
dim(umis) #33538  6965
head(rownames(umis))
gene_trad<-TransEnsembltoSymbol(rownames(umis))

umis_f<-umis[gene_trad[hgnc_symbol!=""]$ensembl_gene_id,]
rownames(umis_f)<-gene_trad[hgnc_symbol!=""]$hgnc_symbol
sum(duplicated(rownames(umis_f))) #8
cbp0_ctrl_all <- CreateSeuratObject(counts = umis_f,project = sample)
cbp0_ctrl_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp0_ctrl_all, pattern = "^MT-")

VlnPlot(object = cbp0_ctrl_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
p1<-VlnPlot(object = cbp0_ctrl_all, features ="percent.mt")+
  geom_hline(yintercept = median(cbp0_ctrl_all$percent.mt)+ 4*mad(cbp0_ctrl_all$percent.mt) )

p2<-VlnPlot(object = cbp0_ctrl_all, features ="nCount_RNA")+
  geom_hline(yintercept = median(cbp0_ctrl_all$nCount_RNA)-2*mad(cbp0_ctrl_all$nCount_RNA))+
  geom_hline(yintercept = 14000)

p3<-VlnPlot(object = cbp0_ctrl_all, features ="nFeature_RNA")+
  geom_hline(yintercept = median(cbp0_ctrl_all$nFeature_RNA)-2*mad(cbp0_ctrl_all$nFeature_RNA) )+
  geom_hline(yintercept = 2900)
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp0_ctrl_QC_cells_metrics.png"))

cbp0_ctrl_qc<-subset(cbp0_ctrl_all,percent.mt<median(cbp0_ctrl_all$percent.mt)+ 4*mad(cbp0_ctrl_all$percent.mt)&
                       nCount_RNA>median(cbp0_ctrl_all$nCount_RNA)-2*mad(cbp0_ctrl_all$nCount_RNA)&nCount_RNA<14000&
                       nFeature_RNA>median(cbp0_ctrl_all$nFeature_RNA)-2*mad(cbp0_ctrl_all$nFeature_RNA)&nFeature_RNA<2900)

cbp0_ctrl_qc #6478 cells

cbp0_ctrl_qc$sample<-"ctrlF547"

saveRDS(cbp0_ctrl_qc,fp(out,"cbp0_ctrl.rds"))



#cbp0_lga####
sample<-"cbp0_lga"
umis<- as.matrix(data.frame(fread("../singlecell/datasets/cbp0/CBP552.csv"),row.names=1))
dim(umis) #33538  3020
umis_f<-umis[gene_trad[hgnc_symbol!=""]$ensembl_gene_id,]
rownames(umis_f)<-gene_trad[hgnc_symbol!=""]$hgnc_symbol
sum(duplicated(rownames(umis_f))) #8
cbp0_lga_all <- CreateSeuratObject(counts = umis_f,project = sample)
cbp0_lga_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp0_lga_all, pattern = "^MT-")

VlnPlot(object = cbp0_lga_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
p1<-VlnPlot(object = cbp0_lga_all, features ="percent.mt")+geom_hline(yintercept = median(cbp0_lga_all$percent.mt)+ 4*mad(cbp0_lga_all$percent.mt) )

p2<-VlnPlot(object = cbp0_lga_all, features ="nCount_RNA")+geom_hline(yintercept = median(cbp0_lga_all$nCount_RNA)-2*mad(cbp0_lga_all$nCount_RNA))+geom_hline(yintercept = 20000)

p3<-VlnPlot(object = cbp0_lga_all, features ="nFeature_RNA")+geom_hline(yintercept = median(cbp0_lga_all$nFeature_RNA)-2*mad(cbp0_lga_all$nFeature_RNA) )+geom_hline(yintercept = 4000)
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp0_lga_QC_cells_metrics.png"))


cbp0_lga_qc<-subset(cbp0_lga_all,percent.mt<median(cbp0_lga_all$percent.mt)+ 4*mad(cbp0_lga_all$percent.mt)&nCount_RNA>median(cbp0_lga_all$nCount_RNA)-2*mad(cbp0_lga_all$nCount_RNA)&nCount_RNA<20000&nFeature_RNA>median(cbp0_lga_all$nFeature_RNA)-2*mad(cbp0_lga_all$nFeature_RNA)&nFeature_RNA<4000)

cbp0_lga_qc #2897 cells

cbp0_lga_qc$sample<-"lgaF552"

saveRDS(cbp0_lga_qc,fp(out,"cbp0_lga.rds"))



#cbp6a####
sample<-"cbp6a"
umis<- Read10X("~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-a/outs/filtered_feature_bc_matrix/")

cbp6a_all <- CreateSeuratObject(counts = umis,project = sample)
cbp6a_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp6a_all, pattern = "^MT-")

VlnPlot(object = cbp6a_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
thr_mt<-median(cbp6a_all$percent.mt)+ 4*mad(cbp6a_all$percent.mt)
thr_ge<-median(cbp6a_all$nFeature_RNA)- 2*mad(cbp6a_all$nFeature_RNA)
thr_rn<-median(cbp6a_all$nCount_RNA)-2*mad(cbp6a_all$nCount_RNA)
min(cbp6a_all$nCount_RNA)

p1<-VlnPlot(object = cbp6a_all, features ="percent.mt")+geom_hline(yintercept = thr_mt )

p2<-VlnPlot(object = cbp6a_all, features ="nCount_RNA")+geom_hline(yintercept = 30000)

p3<-VlnPlot(object = cbp6a_all, features ="nFeature_RNA")+geom_hline(yintercept = 220 )+geom_hline(yintercept = 5000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp6a_QC_cells_metrics.png"))

cbp6a_qc<-subset(cbp6a_all,percent.mt<median(cbp6a_all$percent.mt)+ 4*mad(cbp6a_all$percent.mt)&nCount_RNA<30000&nFeature_RNA>220&nFeature_RNA<5000)

#reassign samples
VlnPlot(cbp6a_qc,c("XIST","RPS4Y1"),group.by = "orig.ident")

sum(cbp6a_qc@assays$RNA@data["XIST",]>0) 
sum(cbp6a_qc@assays$RNA@data["RPS4Y1",]>0) 

cbp6a_qc@meta.data[cbp6a_qc@assays$RNA@data["XIST",]>0&
                  cbp6a_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"Doublet"

cbp6a_qc@meta.data[cbp6a_qc@assays$RNA@data["XIST",]>0&
                  cbp6a_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"ctrlF544"

cbp6a_qc@meta.data[cbp6a_qc@assays$RNA@data["XIST",]==0&
                  cbp6a_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"lgaM556"

cbp6a_qc@meta.data[cbp6a_qc@assays$RNA@data["XIST",]==0&
                  cbp6a_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"Negative"

table(cbp6a_qc@meta.data$sex.ID)
# ctrlF544  Doublet  lgaM556 Negative 
#     2098      456     3596     2435 

cbp6a_qc_s<-subset(cbp6a_qc,sex.ID!="Doublet"&sex.ID!="Negative")

cbp6a_qc_s$sample<-cbp6a_qc_s$sex.ID
table(cbp6a_qc_s$sample)
# ctrlF544  lgaM556 
#     2098     3596
saveRDS(cbp6a_qc_s,fp(out,"cbp6a.rds"))



#cbp6b####
sample<-"cbp6b"
umis<- Read10X("~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-b/outs/filtered_feature_bc_matrix/")

cbp6b_all <- CreateSeuratObject(counts = umis,project = sample)
cbp6b_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp6b_all, pattern = "^MT-")

VlnPlot(object = cbp6b_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
thr_mt<-median(cbp6b_all$percent.mt)+ 4*mad(cbp6b_all$percent.mt)
thr_ge<-median(cbp6b_all$nFeature_RNA)- 2*mad(cbp6b_all$nFeature_RNA)
thr_rn<-median(cbp6b_all$nCount_RNA)-2*mad(cbp6b_all$nCount_RNA)

p1<-VlnPlot(object = cbp6b_all, features ="percent.mt")+geom_hline(yintercept = thr_mt )

p2<-VlnPlot(object = cbp6b_all, features ="nCount_RNA")+geom_hline(yintercept = 30000)

p3<-VlnPlot(object = cbp6b_all, features ="nFeature_RNA")+geom_hline(yintercept = thr_ge )+geom_hline(yintercept = 5000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp6b_QC_cells_metrics.png"))

cbp6b_qc<-subset(cbp6b_all,percent.mt< thr_mt&nCount_RNA<30000&nFeature_RNA>thr_ge&nFeature_RNA<5000)

#reassign samples
VlnPlot(cbp6b_qc,c("XIST","RPS4Y1"),group.by = "orig.ident")

sum(cbp6b_qc@assays$RNA@data["XIST",]>0) 
sum(cbp6b_qc@assays$RNA@data["RPS4Y1",]>0) 

cbp6b_qc@meta.data[cbp6b_qc@assays$RNA@data["XIST",]>0&
                  cbp6b_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"Doublet"

cbp6b_qc@meta.data[cbp6b_qc@assays$RNA@data["XIST",]>0&
                  cbp6b_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"ctrlF545"

cbp6b_qc@meta.data[cbp6b_qc@assays$RNA@data["XIST",]==0&
                  cbp6b_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"lgaM533"

cbp6b_qc@meta.data[cbp6b_qc@assays$RNA@data["XIST",]==0&
                  cbp6b_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"Negative"

table(cbp6b_qc@meta.data$sex.ID)
# ctrlF545  Doublet  lgaM533 Negative 
#      875       66     2437      363 

cbp6b_qc_s<-subset(cbp6b_qc,sex.ID!="Doublet"&sex.ID!="Negative")

cbp6b_qc_s$sample<-cbp6b_qc_s$sex.ID
table(cbp6b_qc_s$sample)
# ctrlF545  lgaM533 
#      875     2437 
saveRDS(cbp6b_qc_s,fp(out,"cbp6b.rds"))


#cbp6c####
sample<-"cbp6c"
umis<- Read10X("~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-c/outs/filtered_feature_bc_matrix/")

cbp6c_all <- CreateSeuratObject(counts = umis,project = sample)
cbp6c_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp6c_all, pattern = "^MT-")

VlnPlot(object = cbp6c_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
thr_mt<-median(cbp6c_all$percent.mt)+ 4*mad(cbp6c_all$percent.mt)
thr_ge<-median(cbp6c_all$nFeature_RNA)- 2*mad(cbp6c_all$nFeature_RNA)
thr_rn<-median(cbp6c_all$nCount_RNA)-2*mad(cbp6c_all$nCount_RNA)

p1<-VlnPlot(object = cbp6c_all, features ="percent.mt")+geom_hline(yintercept = thr_mt )

p2<-VlnPlot(object = cbp6c_all, features ="nCount_RNA")+geom_hline(yintercept = 30000)

p3<-VlnPlot(object = cbp6c_all, features ="nFeature_RNA")+geom_hline(yintercept = 5000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp6c_QC_cells_metrics.png"))

cbp6c_qc<-subset(cbp6c_all,percent.mt< thr_mt&nCount_RNA<30000&nFeature_RNA<5000)

#reassign samples
VlnPlot(cbp6c_qc,c("XIST","RPS4Y1"),group.by = "orig.ident")

sum(cbp6c_qc@assays$RNA@data["XIST",]>0) 
sum(cbp6c_qc@assays$RNA@data["RPS4Y1",]>0) 

cbp6c_qc@meta.data[cbp6c_qc@assays$RNA@data["XIST",]>0&
                  cbp6c_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"Doublet"

cbp6c_qc@meta.data[cbp6c_qc@assays$RNA@data["XIST",]>0&
                  cbp6c_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"ctrlF541"

cbp6c_qc@meta.data[cbp6c_qc@assays$RNA@data["XIST",]==0&
                  cbp6c_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"lgaM526"

cbp6c_qc@meta.data[cbp6c_qc@assays$RNA@data["XIST",]==0&
                  cbp6c_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"Negative"

table(cbp6c_qc@meta.data$sex.ID)
# ctrlF541  Doublet  lgaM526 Negative 
#     2346      114     1550      426 

cbp6c_qc_s<-subset(cbp6c_qc,sex.ID!="Doublet"&sex.ID!="Negative")

cbp6c_qc_s$sample<-cbp6c_qc_s$sex.ID
table(cbp6c_qc_s$sample)
# ctrlF541  lgaM526 
#     2346     1550 

saveRDS(cbp6c_qc_s,fp(out,"cbp6c.rds"))



#cbp7a####
sample<-"cbp7a"
umis<- Read10X("~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-a/outs/filtered_feature_bc_matrix/")

cbp7a_all <- CreateSeuratObject(counts = umis,project = sample)
cbp7a_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp7a_all, pattern = "^MT-")

VlnPlot(object = cbp7a_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
thr_mt<-median(cbp7a_all$percent.mt)+ 4*mad(cbp7a_all$percent.mt)
thr_ge<-median(cbp7a_all$nFeature_RNA)- 2*mad(cbp7a_all$nFeature_RNA)
thr_rn<-median(cbp7a_all$nCount_RNA)-2*mad(cbp7a_all$nCount_RNA)

p1<-VlnPlot(object = cbp7a_all, features ="percent.mt")+geom_hline(yintercept = thr_mt )

p2<-VlnPlot(object = cbp7a_all, features ="nCount_RNA")+geom_hline(yintercept = 35000)

p3<-VlnPlot(object = cbp7a_all, features ="nFeature_RNA")+geom_hline(yintercept = 6000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp7a_QC_cells_metrics.png"))

cbp7a_qc<-subset(cbp7a_all,percent.mt<thr_mt&nCount_RNA<35000&nFeature_RNA<6000)

#reassign samples
VlnPlot(cbp7a_qc,c("XIST","RPS4Y1"),group.by = "orig.ident")
sum(cbp7a_qc@assays$RNA@data["XIST",]>0) 
sum(cbp7a_qc@assays$RNA@data["RPS4Y1",]>0) 
sum(cbp7a_qc@assays$RNA@data["XIST",]>0&cbp7a_qc@assays$RNA@data["RPS4Y1",]>0) #408 doublet


cbp7a_qc@meta.data[cbp7a_qc@assays$RNA@data["XIST",]>0&
                  cbp7a_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"Doublet"

cbp7a_qc@meta.data[cbp7a_qc@assays$RNA@data["XIST",]>0&
                  cbp7a_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"lgaF554"

cbp7a_qc@meta.data[cbp7a_qc@assays$RNA@data["XIST",]==0&
                  cbp7a_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"ctrlM555"

cbp7a_qc@meta.data[cbp7a_qc@assays$RNA@data["XIST",]==0&
                  cbp7a_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"Negative"

table(cbp7a_qc@meta.data$sex.ID)
# ctrlM555  Doublet  lgaF554 Negative 
#     1990      408     2800      573 

cbp7a_qc_s<-subset(cbp7a_qc,sex.ID!="Doublet"&sex.ID!="Negative")

cbp7a_qc_s$sample<-cbp7a_qc_s$sex.ID
table(cbp7a_qc_s$sample)
# ctrlM555  lgaF554 
#     1990     2800 

saveRDS(cbp7a_qc_s,fp(out,"cbp7a.rds"))


#cbp7b####
sample<-"cbp7b"
umis<- Read10X("~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-b/outs/filtered_feature_bc_matrix/")

cbp7b_all <- CreateSeuratObject(counts = umis,project = sample)
cbp7b_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp7b_all, pattern = "^MT-")

VlnPlot(object = cbp7b_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
thr_mt<-median(cbp7b_all$percent.mt)+ 4*mad(cbp7b_all$percent.mt)
thr_ge<-median(cbp7b_all$nFeature_RNA)- 2*mad(cbp7b_all$nFeature_RNA)
thr_rn<-median(cbp7b_all$nCount_RNA)-2*mad(cbp7b_all$nCount_RNA)

p1<-VlnPlot(object = cbp7b_all, features ="percent.mt")+geom_hline(yintercept = thr_mt )

p2<-VlnPlot(object = cbp7b_all, features ="nCount_RNA")+geom_hline(yintercept = 35000)

p3<-VlnPlot(object = cbp7b_all, features ="nFeature_RNA")+geom_hline(yintercept = 6000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp7b_QC_cells_metrics.png"))

cbp7b_qc<-subset(cbp7b_all,percent.mt<thr_mt&nCount_RNA<35000&nFeature_RNA<6000)

#reassign samples
VlnPlot(cbp7b_qc,c("XIST","RPS4Y1"),group.by = "orig.ident")
sum(cbp7b_qc@assays$RNA@data["XIST",]>0) 
sum(cbp7b_qc@assays$RNA@data["RPS4Y1",]>0) 
sum(cbp7b_qc@assays$RNA@data["XIST",]>0&cbp7b_qc@assays$RNA@data["RPS4Y1",]>0) #419 doublet


cbp7b_qc@meta.data[cbp7b_qc@assays$RNA@data["XIST",]>0&
                  cbp7b_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"Doublet"

cbp7b_qc@meta.data[cbp7b_qc@assays$RNA@data["XIST",]>0&
                  cbp7b_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"lgaF521"

cbp7b_qc@meta.data[cbp7b_qc@assays$RNA@data["XIST",]==0&
                  cbp7b_qc@assays$RNA@data["RPS4Y1",]>0,"sex.ID"]<-"ctrlM518"

cbp7b_qc@meta.data[cbp7b_qc@assays$RNA@data["XIST",]==0&
                  cbp7b_qc@assays$RNA@data["RPS4Y1",]==0,"sex.ID"]<-"Negative"

table(cbp7b_qc@meta.data$sex.ID)
# ctrlM518  Doublet  lgaF521 Negative 
#     2009      469     3576      532 

cbp7b_qc_s<-subset(cbp7b_qc,sex.ID!="Doublet"&sex.ID!="Negative")

cbp7b_qc_s$sample<-cbp7b_qc_s$sex.ID
table(cbp7b_qc_s$sample)
# ctrlM518  lgaF521 
#     2009     3576 

saveRDS(cbp7b_qc_s,fp(out,"cbp7b.rds"))

#cbp7c####
sample<-"cbp7c"
umis<- Read10X("~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-c/outs/filtered_feature_bc_matrix/")

cbp7c_all <- CreateSeuratObject(counts = umis,project = sample)
cbp7c_all[["percent.mt"]] <- PercentageFeatureSet(object = cbp7c_all, pattern = "^MT-")

VlnPlot(object = cbp7c_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#take the 4 four median absolute deviations above the median
thr_mt<-median(cbp7c_all$percent.mt)+ 4*mad(cbp7c_all$percent.mt)
thr_ge<-median(cbp7c_all$nFeature_RNA)- 2*mad(cbp7c_all$nFeature_RNA)
thr_rn<-median(cbp7c_all$nCount_RNA)-2*mad(cbp7c_all$nCount_RNA)


p1<-VlnPlot(object = cbp7c_all, features ="percent.mt")+geom_hline(yintercept = thr_mt )

p2<-VlnPlot(object = cbp7c_all, features ="nCount_RNA")+geom_hline(yintercept = 35000)

p3<-VlnPlot(object = cbp7c_all, features ="nFeature_RNA")+geom_hline(yintercept = 6000 )
p1|p2|p3+plot_layout(guides = "collect")
ggsave(fp(out,"cbp7c_QC_cells_metrics.png"))

cbp7c_qc<-subset(cbp7c_all,percent.mt<thr_mt&nCount_RNA<35000&nFeature_RNA<6000)
cbp7c_qc #2823 cells

cbp7c_qc$sample<-"ctrlM537"

saveRDS(cbp7c_qc,fp(out,"cbp7c.rds"))

#### INTEGRATION ####
#based on Hematomap  (see 01-Make_Hematomap )
options(future.globals.maxSize = 50000 * 1024^2)
out<-"outputs/06-integr_singlecell_cbps"
source("../methyl/scripts/utils/new_utils.R")
library(Seurat)
library(parallel)

hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")
hmap
DefaultAssay(hmap)<-"integrated"
hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

cbps_run<-c("cbp0_ctrl","cbp0_lga",
            paste0("cbp",2:4),
            ps("cbp6",c("a","b","c")),ps("cbp7",c("a","b","c")),
            "cbp8")
length(cbps_run)#12
  
cbps_list<-lapply(cbps_run, function(run)readRDS(fp(out,ps(run,".rds"))))
cbps_list

cbps_list<-mclapply(cbps_list,function(x){
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
  },mc.cores = 6)

cbps_list<-mclapply(cbps_list, SCTransform,vars.to.regress=c("percent.mt","CC.Difference"),
                  return.only.var.genes=F, 
                  method = "glmGamPoi",mc.cores = 6)


anchors <- list()
for (i in 1:length(cbps_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = hmap,
    query = cbps_list[[i]],
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )
}


for (i in 1:length(cbps_list)) {
  cbps_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = cbps_list[[i]],
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )
}


# Merge the queries
cbps <- merge(cbps_list[[1]], cbps_list[2:length(cbps_list)],merge.dr = c("ref.pca","ref.umap"))
p<-DimPlot(cbps, reduction = "ref.umap", group.by =  "predicted.cell_type", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
ggsave(fp(out,"predicted_cell_type.png"),plot=p)


DimPlot(cbps, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()



#add metadata

cbps[["group"]]<-str_extract(cbps@meta.data$sample,"ctrl|lga|iugr")
cbps[["sex"]]<-str_extract(cbps@meta.data$sample,"M|F")
cbps[["group_sex"]]<-paste0(cbps@meta.data$group,cbps@meta.data$sex)
cbps$hto<-cbps$orig.ident%in%c("cbp2","cbp3","cbp4","cbp8")
cbps$batch<-cbps$orig.ident

cbps[["group_hto"]]<-paste0(cbps@meta.data$group,cbps@meta.data$hto)
cbps[["sample_hto"]]<-paste0(cbps@meta.data$sample,cbps@meta.data$hto)
cbps[["ambigous"]]<-cbps@meta.data$sample%in%c("iugrM558","lgaF559")

cbps[["cell_type_hmap"]]<-cbps$predicted.cell_type
cbps[["lineage_hmap"]]<-cbps$predicted.lineage

cbps$differentiated<-cbps$lineage_hmap%in%c("Mk/Er","18","DC","T cell","B cell")


#denovo umap
cbps <- RunUMAP(cbps, reduction = 'ref.pca', dims = 1:30,reduction.name = "denovo.umap",n.components = 2)
DimPlot(cbps, group.by = 'cell_type_hmap',reduction = "denovo.umap", label = TRUE)
DimPlot(cbps, group.by = 'lineage_hmap',reduction = "denovo.umap", label = TRUE)

saveRDS(cbps,fp(out,"cbps.rds"))

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps.rds")

####CHECK INTEGR OK####
cbps_f<-subset(cbps,lineage_hmap!="18"&ambigous==F&group!="iugr")
rm(cbps)

#check qu'on a bien tout

#all
cbps_f#47995 cells
length(unique(cbps_f$sample_hto)) #27 samples_conditions
length(unique(cbps_f$sample)) #23 samples

#n samples by group
mtd<-data.table(cbps_f@meta.data,keep.rownames = "bc")
mts<-unique(mtd,by=c("sample_hto"))
table(mts$hto,mts$group)
  #       ctrl  lga
  # FALSE    7   6
  # TRUE     8   6

#n of cells
table(mtd$hto,mtd$group)
  #       ctrl   lga
  # FALSE 18520 16791
  # TRUE   5823  6861

#n of HSC 
table(mtd[lineage_hmap=="HSC"]$hto,mtd[lineage_hmap=="HSC"]$group)
  #      ctrl  lga
  # FALSE 3485 4200
  # TRUE  2093 1959

#check subpop distr
mtd$lineage_hmap<-factor(mtd$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Lymphoid","B cell","T cell","Erythro-Mas","Mk/Er","Myeloid","DC"))

ggplot(mtd)+geom_bar(aes(x=sample_hto,fill=lineage_hmap),position = "fill")
mtd[,n.sample:=.N,by="sample_hto"]
mtd[,pct.lin:=.N/n.sample,by=c("sample_hto","lineage_hmap")]
ggplot(unique(mtd,by=c("sample_hto","lineage_hmap")))+geom_boxplot(aes(x=hto,y=pct.lin,fill=group))+facet_wrap("lineage_hmap",scales = "free")

mtd[,pct.ct:=.N/n.sample,by=c("sample_hto","cell_type_hmap")]
ggplot(unique(mtd,by=c("sample_hto","cell_type_hmap")))+geom_boxplot(aes(x=hto,y=pct.ct,fill=group))+facet_wrap("cell_type_hmap",scales = "free")

unique(mtd[group=="ctrl"&pct.lin<0.1&lineage_hmap=="HSC"],by="sample_hto")

ggplot(mtd[sample%in%c("ctrlM555","ctrlM518","ctrlM537")])+
  geom_bar(aes(x=hto,fill=lineage_hmap),position = "fill")+
  facet_wrap("sample",scales = "free")
table(mtd[sample=="ctrlM537",.(lineage_hmap,hto)])

fwrite(mtd,fp(out,"metadata_cbps_filtered.csv.gz"))

mtd[,toy_data:=bc%in%sample(bc,ceiling(.N/10)),.(sample_hto,lineage_hmap)]
mtd[toy_data==T]
saveRDS(cbps_f[,mtd[toy_data==T]$bc],fp(out,"cbps_4k.rds"))

saveRDS(cbps_f,fp(out,"cbps_filtered.rds"))

#check good assignmenet

VlnPlot(cbps_f,"predicted.lineage.score",group.by = "predicted.lineage",pt.size = 0)
#attention a HSC2, 3 et 4, GMP cycle, et MPP Ery
