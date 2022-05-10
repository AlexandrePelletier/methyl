library(Seurat)
source("scripts/utils/new_utils.R")
out<-"outputs/05-make_hematomap/"

cbp<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/CD34_CTRL_CBP547/2020-05-28_seurat_obj2_final.rds")
saveRDS(cbp,fp(out,"ctrlF547.rds"))

cbp<-subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP6-a/cbp6a.rds"),sample=="ctrlF544")
saveRDS(cbp,fp(out,"ctrlF544.rds"))

cbp<-subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP6-b/cbp6b.rds"),sample=="ctrlF545")
saveRDS(cbp,fp(out,"ctrlF545.rds"))



cbp<-subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/CBP6-c/cbp6c.rds"),sample=="ctrlF541")
saveRDS(cbp,fp(out,"ctrlF541.rds"))

cbp<-subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/cbp7a/cbp7a_singlet.rds"),sample=="ctrlM555")
saveRDS(cbp,fp(out,"ctrlM555.rds"))

cbp<-subset(readRDS("../singlecell/outputs/01-Analyses_Individuelles/cbp7b/cbp7b_singlet.rds"),sample=="ctrlM518")
saveRDS(cbp,fp(out,"ctrlM518.rds"))


cbp<-readRDS("../singlecell/outputs/01-Analyses_Individuelles/cbp7c/cbp7c_singlet.rds")
saveRDS(cbp,fp(out,"ctrlM537.rds"))
