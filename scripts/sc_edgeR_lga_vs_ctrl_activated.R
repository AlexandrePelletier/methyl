#scDEGs ctrl vs lga stressed by lineage
library(Seurat)
library(edgeR)
library(parallel)
source("../methyl/scripts/utils/new_utils.R")

out<-"outputs/09-LGA_vs_Ctrl_Activated"
dir.create(out,recursive=T)

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

cbps<-subset(cbps,hto==T)



run_edgeRQLF <- function(counts,mtd.data,design,contrast) {
  dge <- DGEList(counts, samples = mtd.data)
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit,contrast=contrast)
  tt <- topTags(qlf, n = Inf)
  return( data.table(gene = rownames(tt$table),
                         p_val = tt$table$PValue,
                         p_val_adj = tt$table$FDR,
                         avg_logFC= tt$table$logFC)
         )
 
}


Idents(cbps)<-"lineage_hmap"

res_lin<-Reduce(rbind,lapply(levels(cbps),function(lin){
  cbps_sub<-subset(cbps,idents=lin)
  message("lin")
  #get and filter counts
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  counts<-counts[rowSums(counts>0)>=0.20*ncol(counts),]
  
  
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mtd$cdr <- scale(colMeans(counts > 0))
  design <- model.matrix(~cdr+sample,data = mtd)
  
  
  # calculate the vector of coefficient weights in lga
  lga <- colMeans(design[mtd$group == "lga", ])
  
  # calculate the vector of coefficient weights in ctrl
  ctrl <- colMeans(design[mtd$group == "ctrl", ])
  
  
  res<-run_edgeRQLF(counts=counts,
                    data.frame(mtd,row.names = "bc"),
                    design,
                    contrast=lga - ctrl)
  return(res[,lineage_hmap:=lin])

  }))

fwrite(res_lin,fp(out,"res_scEdgeR_by_lineage.csv.gz"),sep = ";")
