out<-"outputs/23-multiGRN"

dir.create(out)

source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)



#I) InferCRE: ####
#function which return CREs, intensity and signifiance of assocation between gene expression and neighbourood peaks.
#1) the function
InferCRE<-function(object,
                   peak.assay="peaks",
                   expression.assay="SCT",
                   expression.slot="data",
                   gene.coords=NULL,
                   distance=5e+05,
                   min.distance=NULL,
                   min.cells=10,
                   method="pearson",
                   genes.use=NULL,
                   n_sample=200,
                   pvalue_cutoff=0,
                   score_cutoff=0,
                   gene.id=FALSE,
                   return.object=FALSE){
  ##function which used LinkPeaks Signac function and 
  #return data.table (default) or updated seurat object 
  #presenting peaks gene association with intensity and signifiance of the assocation.
  DefaultAssay(object)<-peak.assay
    object <- LinkPeaks(
    object = object,
    peak.assay = peak.assay,
    expression.assay = expression.assay,
    expression.slot=expression.slot,
    gene.coords=gene.coords,
    distance=distance,
    min.distance=min.distance,
    min.cells=min.cells,
    method=method,
    genes.use=genes.use,
    n_sample=n_sample,
    score_cutoff = score_cutoff,
    pvalue_cutoff = pvalue_cutoff,
    gene.id=gene.id,
    )
  if(return.object){
    return(object)
  }else{
    peak_links_dt<-data.table(as.data.frame(Links(object)))
    return(peak_links_dt)
      }
}

#2) the test
#do not run with a lot of genes because take a while but test with 3 genes,

cbl12<-readRDS("outputs/21-CREs_inference_with_sc_multiomics_cbls/cbl12.rds")

some_lineage_markers<-c("MPO","VPREB1","LTB","GATA1","GATA2","EGR1","CD99","SELL","CDK6",
                        "HBD","FOS","CEBPA","AVP","DUSP2","ID1")
peak_genes_links<-InferCRE(cbl12,genes.use =some_lineage_markers )


peak_genes_links
table(table(peak_genes_links[order(peak,gene)][peak%in%peak[duplicated(peak)]]$peak))
#   2   3   4   5   6 
# 287  47  10   1   1 
#1 peak are associated to 6 gene expression.

peak_genes_linksf<-peak_genes_links[gene%in%some_lineage_markers]#because cbl12 already process for all genes
peak_genes_linksf[order(peak,gene)]

#II) AnnotCRE()
#function which annotate CREs, with TF motif annotation 

#need get TF motifs on each CREs
#1) the function

#get motif

peaks_motifs<-GetMotif(cbl12,peaks = peak_genes_linksf$peak)
peaks_motifs

AnnotCRE<-function(peak_genes_links){
  #function which annotate CREs, with TF motif annotation 
  #return TF motif if located in the peak, and score of TF role/presence
  #based on fold enrichment and pvalue of the enrichment for this TF motif on CREs compared to comparable non-CRE peaks.
  #Plus, add a supplemental columns validating or not TF binding based on TF footprint on this peak.
  require("Signac")
  require("data.table")
  #get motif
  

  
  #calculate TF motifs enrichment in CREs 
  
  #validate TF binding in peak based on TF footprint
  
  #return TF-genes data.table
  
  return(tf_peak_genes_links)
  
}

#2) the test
