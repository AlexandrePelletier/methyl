out<-"outputs/17-LGA_vs_Ctrl_EGR1_interaction"
dir.create(out)
source("scripts/utils/new_utils.R")

#renv::install("arrow")

####Functions####
RunScenic<-function(seurat_object,
                    dir,
                    cisTarget_dir,
                    assay=NULL,
                    scenic.options=NULL,
                    run.coexpression=TRUE,
                    run.steps=1:3,
                    min.count=NULL,
                    min.cells=NULL,
                    data.title="SCENIC",
                    organism="hgnc",
                    verbose=TRUE,
                    nCores=20,
                    seed=123){
  require("withr")
  require("SCENIC")
  require("Seurat")
  
  #create scenic env
  dir.create(dir)
  withr::with_dir(dir,{
  
    dir.create("output")
    dir.create("int")
  
    if(is.null(scenic.options)){
      #data extraction
    if(is.null(assay))assay=DefaultAssay(seurat_object)
    
    exprMat<-as.matrix(seurat_object@assays[[assay]]@data)
    
    cellInfo<-seurat_object@meta.data
    saveRDS(cellInfo, file="int/cellInfo.Rds")
    data(defaultDbNames)
    dbs <- defaultDbNames[[organism]]
  
    #initialiser scenic
    scenicOptions <- initializeScenic(org=organism, dbDir=cisTarget_dir, dbs=dbs, datasetTitle=data.title, nCores=nCores) 
  
    #filtering : 
    #keeps only genes 1) with at least x UMI counts across all samples, and 2) detected in at least x% of the cells are kept
    if(is.null(min.count))min.count=0.03*ncol(exprMat)
    if(is.null(min.cells))min.cells=ncol(exprMat)*.01
  
    
    genesKept <- geneFiltering(as.matrix(seurat_object@assays$RNA@counts), scenicOptions=scenicOptions,
                               minCountsPerGene=min.count,
                               minSamples=min.cells) #(Adjust minimum values according to your dataset)
    
  
    exprMat_filtered <- exprMat[rownames(exprMat)%in%genesKept, ]
    message(nrow(exprMat_filtered)," genes kept after filtering")
    saveRDS(exprMat_filtered,"int/expr_mat_filtered.rds")
  
    }
    if(run.coexpression){
      
    
      #run correl
      runCorrelation(exprMat_filtered, scenicOptions)
  
      # Run GENIE3
      runGenie3(exprMat_filtered, scenicOptions) #run script run_genie3.R in a job because take long time 
  
    
    }
    
    scenicOptions@settings$verbose <- verbose
    scenicOptions@settings$nCores <- nCores
    scenicOptions@settings$seed <- seed
    saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
    
  
    ## run: 
    if(1%in%run.steps)scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
    if(2%in%run.steps)scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
    if(3%in%run.steps)scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)
  
    saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status
    return(scenicOptions)
  })
}



###ANALYSIS####
#LGA vs Ctrl GENIEbased tf-gene Interaction
#get tf-gene interaction df + regulons for 1 toy data
cbps_toy<-readRDS("outputs/06-integr_singlecell_cbps/cbps_4k.rds")
mtd<-data.table(cbps_toy@meta.data,keep.rownames = "bc")
mtd[,cells_choosed:=bc%in%sample(bc,20,replace = T),by=.(cell_type_hmap)] #20 cells max by sample_hto and cell_type


mtd[cells_choosed==T] #331 cells

cbps_toyf<-cbps_toy[,colnames(cbps_toy)%in%mtd[cells_choosed==T]$bc]

scenics<-RunScenic(cbps_toyf,dir = fp(out,"cbps_toy"),
                   cisTarget_dir = "/disks/DATATMP/PhD_AlexandrePelletier/marcos/ref/cisTarget_databases/human/")

tfg_toy<-data.table(readRDS("outputs/17-LGA_vs_Ctrl_EGR1_interaction/cbps_toy/int/1.4_GENIE3_linkList.Rds"))


#for 50*5k LGA cells
#run 17A-run_LGA_SCENIC_run
#for 50*5k CTRL cells
#run 17B

res_inter<-rbind(Reduce(rbind,lapply(1:23,function(i)data.table(readRDS(ps("outputs/17-LGA_vs_Ctrl_EGR1_interaction/lga_runs/run_",i,"/int/1.4_GENIE3_linkList.Rds")))[,run:=i]))[,group:="lga"],
                 Reduce(rbind,lapply(1:23,function(i)data.table(readRDS(ps("outputs/17-LGA_vs_Ctrl_EGR1_interaction/ctrl_runs/run_",i,"/int/1.4_GENIE3_linkList.Rds")))[,run:=i]))[,group:="ctrl"])

ggplot(res_inter[TF=="EGR1"&Target=="SOCS3"])+geom_boxplot(aes(x=group,y=weight))

ggplot(res_inter[TF=="KLF2"&Target=="SOCS3"])+geom_boxplot(aes(x=group,y=weight))

res_inter[,weight_norm:=scale(weight),by=.(TF,Target)]
ggplot(res_inter)+geom_boxplot(aes(x=group,y=weight_norm))

ggplot(res_inter)+geom_density(aes(x=log10(weight),col=group))

ggplot(res_inter)+geom_density(aes(x=weight_norm,col=group))


res_inter[,n.obs:=.N,.(group,TF,Target)]
plot(density(res_inter$n.obs))
res_inter[n.obs==23]
res_inter[n.obs>=20]

tftargets<-fread("outputs/16-GRN_final/tf_target_interactions.csv")
tftargets[,TF:=tf][,Target:=target]
tftargets[!(extended)]
res_interf<-merge(res_inter,tftargets[!(extended),.(TF,Target)])
res_interf[,inter_to_test:=TF%in%c("KLF2","KLF4","EGR1")]
res_interf[(inter_to_test),inter_to_test:=all(c(any(n.obs>=20&group=="lga"),any(n.obs>=20&group=="ctrl"))),by=.(TF,Target)]
unique(res_interf[(inter_to_test)],by = c("TF","Target")) #130
res_interf[(inter_to_test),pvalue:=wilcox.test(weight[group=="lga"],weight[group=="ctrl"])$p.value,by=.(TF,Target)]
res_interf[(inter_to_test),mean.lga:=mean(weight[group=="lga"]),by=.(TF,Target)]
res_interf[(inter_to_test),mean.ctrl:=mean(weight[group=="ctrl"]),by=.(TF,Target)]
res_interf[(inter_to_test),avg_diff:=mean.lga-mean.ctrl,by=.(TF,Target)]
res_interf[(inter_to_test),fold.change:=mean.lga/mean.ctrl,by=.(TF,Target)]
res_interf[(inter_to_test),avg_log2FC:=log2(fold.change),by=.(TF,Target)]

res_interf_stats<-unique(res_interf[(inter_to_test)],by=c("TF","Target"))
res_interf_stats[,padj:=p.adjust(pvalue,method = "BH")]
res_interf_stats[padj<0.05]#91/130

res_interf_stats[padj<0.05][order(-abs(avg_diff))]
res_interf_stats[padj<0.001&TF=="EGR1"][order(-abs(avg_diff))]

res_interf_stats[TF=='EGR1'&Target=="SOCS3"]
res_interf_stats[TF=='KLF2'&Target=="SOCS3"]

ggplot(res_interf_stats,aes(x=avg_log2FC,y=-log10(padj)))+
  geom_point(aes(col=padj<0.001))+
  geom_label_repel(aes(label=ifelse(padj<0.001,Target,"")),
                   max.overlaps = 500,
                   size = 3)+
  facet_wrap("TF",scales = "free_x")+
  scale_color_manual(values=c('grey','red'))
res_interf_stats[padj<0.001]

res_interf_stats[Target=="EGR1"]

#for samples independantly with RNA assay
#run 17C
tfgs<-fread(fp(out,"sample_runs/tf_targets_interactions_by_samples.csv.gz"))
tfgs

#EGR1, KLF2, KLF4, reduce interactions with target in LGA ?
tf_of_int<-c("EGR1","KLF2","KLF4")
tfgs_int<-tfgs[TF%in%tf_of_int]


#filter for targets validating with scATAC-seq
tftargets<-fread("outputs/16-GRN_final/tf_target_interactions.csv")
tftargets[,TF:=tf][,Target:=target]
tftargets[!(extended)]

tfgs_intf<-merge(tfgs_int,tftargets[!(extended),.(TF,Target)])
tfgs_intf[,group:=str_extract(sample,"ctrl|lga")]
tfgs_intf[,hto:=as.logical(str_extract(sample,"FALSE|TRUE"))]

#sample do not have same amount of cells, Is there a ncell influence on interac wieght ?
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtd[,avg.rna.count:=mean(nCount_RNA),by="sample_hto"]
mtd[,avg.gene.count:=mean(nFeature_RNA),by="sample_hto"]
tfgs_intf[,sample_hto:=sample]
tfgs_intfm<-merge(tfgs_intf,unique(mtd[,.(sample_hto,n.sample,avg.rna.count,avg.gene.count)]),by="sample_hto")

ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(n.sample),y=weight)) #seems to be a negative influence

#avg RNA count ?
ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=weight)) #positive influence

#==>need normalize weight or include in the model.

#normalize weight for boxplot
?scale
tfgs[,scaled.weight:=scale(weight,center = F),by="sample"]
tfgs_int<-tfgs[TF%in%tf_of_int]
tfgs_intf<-merge(tfgs_int,tftargets[!(extended),.(TF,Target)])
tfgs_intf[,group:=str_extract(sample,"ctrl|lga")]
tfgs_intf[,hto:=as.logical(str_extract(sample,"FALSE|TRUE"))]
tfgs_intf[,sample_hto:=sample]

tfgs_intfm<-merge(tfgs_intf,unique(mtd[,.(sample_hto,n.sample,avg.rna.count,avg.gene.count)]),by="sample_hto")

ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(n.sample),y=scaled.weight)) #seems to be a negative influence
ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=scaled.weight)) #positive influence


p1<-ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=scaled.weight)) 
p2<-ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=weight)) 
p1/p2 #doesnt change a lot

tfgs_intfm[,avg.scaled.weight:=mean(scaled.weight),by=c("sample_hto","TF")]
ggplot(tfgs_intfm)+geom_boxplot(aes(x=hto,y=avg.scaled.weight,fill=group))+facet_wrap("TF") 

#include in a model for stat
tfgs_intfm[,avg.weight:=mean(weight),by=c("sample_hto","TF")]

summary(lm(avg.weight~group*hto+avg.rna.count+n.sample,data = unique(tfgs_intfm[TF=="KLF4"],by="sample_hto")))
summary(lm(avg.weight~group*hto+avg.rna.count+n.sample,data = unique(tfgs_intfm[TF=="EGR1"],by="sample_hto")))
summary(lm(avg.weight~group*hto+avg.rna.count+n.sample,data = unique(tfgs_intfm[TF=="KLF2"],by="sample_hto")))


#scaled at regulons level
tfgsf<-merge(tfgs,tftargets[!(extended),.(TF,Target)])

tfgsf[,scaled.weight:=scale(weight,center = F),by="sample"]
tfgs_intf<-tfgsf[TF%in%tf_of_int]
tfgs_intf[,group:=str_extract(sample,"ctrl|lga")]
tfgs_intf[,hto:=as.logical(str_extract(sample,"FALSE|TRUE"))]
tfgs_intf[,sample_hto:=sample]

tfgs_intfm<-merge(tfgs_intf,unique(mtd[,.(sample_hto,n.sample,avg.rna.count,avg.gene.count)]),by="sample_hto")

ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(n.sample),y=scaled.weight)) #ok
ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=scaled.weight)) #


p1<-ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=scaled.weight)) 
p2<-ggplot(tfgs_intfm)+geom_boxplot(aes(x=as.factor(avg.rna.count),y=weight)) 
p1/p2 #better but not perfect

tfgs_intfm[,avg.scaled.weight:=mean(scaled.weight),by=c("sample_hto","TF")]
ggplot(tfgs_intfm)+geom_boxplot(aes(x=hto,y=avg.scaled.weight,fill=group))+facet_wrap("TF") 

#include in a model for stat
tfgs_intfm[,avg.weight:=mean(weight),by=c("sample_hto","TF")]

summary(lm(avg.weight~group*hto+avg.rna.count+n.sample,data = unique(tfgs_intfm[TF=="KLF4"],by="sample_hto")))
summary(lm(avg.weight~group*hto+avg.rna.count+n.sample,data = unique(tfgs_intfm[TF=="EGR1"],by="sample_hto")))
summary(lm(avg.weight~group*hto+avg.rna.count+n.sample,data = unique(tfgs_intfm[TF=="KLF2"],by="sample_hto")))


#CCL : overall doesnt see interaction differences between LGA and control, maybe too much noise/technical bias

ggplot(tfgs_intfm[TF=="EGR1"&Target=="SOCS3"])+geom_boxplot(aes(x=hto,y=weight,fill=group))

ggplot(tfgs_intfm[TF=="KLF2"&Target=="SOCS3"])+geom_boxplot(aes(x=hto,y=weight,fill=group))

