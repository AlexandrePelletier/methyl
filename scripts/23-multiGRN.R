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

#calculate TF motifs enrichment in CREs 
#use the signac motif enrichment function, 
length(unique(peak_genes_links$peak)) #4170
motif_enrich<-FindMotifs(cbl12,features =peak_genes_links$peak,background = 40000 ) #Matching GC.percent distribution
motif_enrich<-data.table(motif_enrich)
motif_enrich
peaks_motifs[,motif.name:=as.character(motif.name)]
pme<-merge(peaks_motifs,motif_enrich,all.x=T)
pme[fold.enrichment>1.25&pvalue<0.001] #764/2331 peak-TF motif (32%) have TF motif enriched in the CREs
#number of TF enriched / number of TF tot :
length(unique(pme[fold.enrichment>1.25&pvalue<0.001]$motif))
length(unique(pme$motif))#141/584 (24%)

#validate TF binding in peak based on TF footprint
some_lineage_markers<-c("MPO","VPREB1","LTB","GATA1","GATA2","EGR1","CD99","SELL","CDK6",
                        "HBD","FOS","CEBPA","AVP","DUSP2","ID1")

library(BSgenome.Hsapiens.UCSC.hg38)
cbl12 <- Footprint(
  object = cbl12,
  motif.name = c("GATA2", "CEBPA", "EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in_peaks=TRUE
)

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Erythro-Mas","Lymphoid","expected")

p <- PlotFootprint(cbl12, features = c("GATA2", "CEBPA", "EBF1"),idents = lins)
p + patchwork::plot_layout(ncol = 1)
foot_dt<-data.table(GetFootprintData(cbl12,features = c("GATA2", "CEBPA", "EBF1") ))
foot_dtf<-foot_dt[abs(position)<100]
foot_dtf[,avg.norm.count:=mean(count),by=.(group,class,feature)]
fs<-unique(foot_dtf,by=c("group","class","feature"))
fs[is.na(group),group:="expected"]
ggplot(fs[group%in%lins])+geom_col(aes(x=feature,y=avg.norm.count,fill=group),position = "dodge")
?Footprint

#for some lineage TFs

some_lineage_tfs<-c("SPI1","CEBPA",
                   "TCF3","STAT1",
                   "GATA1","GATA2","TAL1",
                   "HOXA9","HOXB4",
                   "EGR1","FOS","JUNB","GATA3","KLF2")
some_lineage_tfsf<-intersect(unique(pme$motif.name),some_lineage_tfs)
cbl12
cbl12 <- Footprint(
  object = cbl12,
  motif.name = some_lineage_tfsf,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in_peaks=TRUE,upstream = 100,downstream = 100
)

foot_dt<-data.table(GetFootprintData(cbl12,features = some_lineage_tfsf ))
foot_dt[is.na(group),group:="expected"]
foot_dt[,avg.norm.count:=mean(norm.value),by=.(group,feature)]
foot_dt[,min.close.count:=min(norm.value[abs(position)<6]),by=.(group,feature)]
table(foot_dt[abs(position)<6]$group,foot_dt[abs(position)<6]$feature)

fwrite(foot_dt,fp(out,"footprinting_some_lineage_tf_cbl12.csv.gz"))

fs<-unique(foot_dt,by=c("group","feature"))
ggplot(fs[group%in%lins])+
  geom_col(aes(x=feature,y=avg.norm.count,fill=group),position = "dodge")+coord_cartesian(ylim = c(0.8,1.2))



fs[,dt_count:=avg.norm.count-min.close.count,by=.(group,feature)]
ggplot(fs[group%in%lins])+
  geom_col(aes(x=feature,y=dt_count,fill=group),position = "dodge")

 PlotFootprint(cbl12, features = c("EGR1", "KLF2", "CEBPA"),idents = lins)
fs[,dt_count_expected:=dt_count[group=="expected"][1],by="feature"]

fs[,dt_dt_count:=dt_count-dt_count_expected,by=.(group,feature)]
fs[,deux_dt_dt_count:=2^dt_dt_count]



ggplot(fs[group%in%lins])+
  geom_col(aes(x=feature,y=avg.norm.count,fill=group),position = "dodge")+coord_cartesian(ylim = c(0.8,1.2))

ggplot(fs[group%in%lins])+geom_col(aes(x=feature,y=dt_dt_count,fill=group),position = "dodge")

#a continuer : 1) comment est calculer l'expected ? 2) le gap a +/-5pb du TF est il biologically signif ?
#=> 3) include in annot TF CREs
#Can we have a measure of Tn5 insertion /enrichment at peak level ??
#=> use tiles plot to see



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
