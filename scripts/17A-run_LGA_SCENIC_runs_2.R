#SCENIC for 50*5k LGA cells
out<-"outputs/17-LGA_vs_Ctrl_EGR1_interaction/lga_runs"
dir.create(out)

source("scripts/utils/new_utils.R")
library(SCENIC)
library(Seurat)

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


cbps<-readRDS("outputs/10A-classical_integr/cbps0-8_clean.rds")
DefaultAssay(cbps)<-"integrated"
mtd<-data.table(cbps@meta.data,keep.rownames = "bc")

tfg_lga<-Reduce(rbind,lapply(44:50, function(i){
  print(i)
  seed<-i
  set.seed(seed)
  mtd[group=="lga",cells_choosed:=bc%in%sample(bc,50,replace = T),by=.(sample_hto,cell_type_hmap)] #50 cells max by sample_hto and cell_type
  mtd[cells_choosed==T,cells_choosed:=bc%in%sample(bc,5000)] #to have exactly 5k cells

  nrow(mtd[cells_choosed==T] )
  cbpsf<-cbps[,colnames(cbps)%in%mtd[cells_choosed==T]$bc]
  scenics<-RunScenic(cbpsf,
                     dir = fp(out,paste0("run_",i)),
                     run.steps = 1,
                     seed = i,
                   cisTarget_dir = "/disks/DATATMP/PhD_AlexandrePelletier/marcos/ref/cisTarget_databases/human/")

  tfg<-data.table(readRDS(fp(out,paste0("run_",i),"int/1.4_GENIE3_linkList.Rds")))
  return(tfg[,run:=i])

  
  }))

fwrite(tfg_lga,fp(out,"tf_targets_interactions_50_lga_runs_5kcells.csv.gz"))
