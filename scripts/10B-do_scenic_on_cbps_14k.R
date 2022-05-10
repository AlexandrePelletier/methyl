
library(Seurat)
library(SCENIC)
library(here)
source("../methyl/scripts/utils/new_utils.R")

working_dir<-here("outputs/10-SCENIC/cbps_14k")
dir.create(working_dir)
dbDir <- "/disks/DATATMP/PhD_AlexandrePelletier/marcos/ref/cisTarget_databases/human/"
organism <- "hgnc" 
myDatasetTitle <- "SCENIC on cbps 0-8 14k cells"

dir.create(working_dir)
setwd(working_dir)


 # RcisTarget databases location

dir.create("output")
dir.create("int")

cbps<-readRDS("outputs/10A-classical_integr/cbps0-8_clean.rds")

cbps[["sample_hto"]]<-paste(cbps$sample,cbps$hto,sep="_")

mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd[,cells_choosed:=bc%in%sample(bc,50,replace = T),by=.(sample_hto,cell_type)] #50 cells max by sample_hto and cell_type

mtd[,cells_choosed:=bc%in%sample(bc,50,replace = T),by=.(sample_hto,cell_type)] #50 cells max by sample_hto and cell_type

mtd[cells_choosed==T] #14k cells

cbps<-cbps[,colnames(cbps)%in%mtd[cells_choosed==T]$bc]

exprMat<-as.matrix(cbps@assays$integrated@data)
dim(exprMat) #3000 47376

cellInfo<-cbps@meta.data[,c("nFeature_RNA","nCount_RNA","lineage","sample_hto")]

saveRDS(cellInfo, file="int/cellInfo.Rds")
colVars <- list(lineage1=c("HSC"="forestgreen", 
                           "Erythroid"="darkorange", 
                           "Lymphoid"="magenta4", 
                           "Myeloid"="hotpink", 
                           "MPP"="orange",
                           "B cell"="blue2",
                           "T cell"="green3",
                           "unknown-1"="red",
                           "unknown-2"="grey"))
colVars$lineage <- colVars$lineage[intersect(names(colVars$lineage), cellInfo$lineage)]
saveRDS(colVars, file="int/colVars.Rds")



 # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[organism]]
scenicOptions <- initializeScenic(org=organism, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

#CO-EXPRESSION NETWORK
#infer potential transcription factor targets based on the expression data
#input: expression matrix (filtered), and a list of transcription factors
#Output :  co-expression modules based on output of GENIE3 and correlation matrix
# i.e  list of potential targets for each TF

#1) filtering : 
#keeps only genes 1) with at least 6 UMI counts across all samples, and 2) detected in at least 1% of the cells are kept 

genesKept <- geneFiltering(as.matrix(cbps@assays$RNA@counts), scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01) #(Adjust minimum values according to your dataset)


exprMat_filtered <- exprMat[rownames(exprMat)%in%genesKept, ]
message(nrow(exprMat_filtered)," genes kept after filtering")

saveRDS(exprMat_filtered,"int/expr_mat_filtered.rds")

#2) correlation (can be parralelized with GENIE3)
#to distinguish potential activation from repression, split targets into positive/negative correlated target
#=> Spearman correlation between the TF and the potential target
runCorrelation(exprMat_filtered, scenicOptions)

#3) GENIE3 [take ++ time]

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions) #run script run_genie3.R in a job because take long time 

#4) Build and score the GRN
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
#exprMat_log <- log2(exprMat+1)


scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 20
scenicOptions@settings$seed <- 123
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# run: 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status





