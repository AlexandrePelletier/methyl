
#linkPeaks for every closest genes of background used for enrichment (peaks test for DA pct.det>0.05)
out<-"outputs/27-multiomics_cbl_stim"
source("scripts/utils/new_utils.R")

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)

atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs$group<-ifelse(str_detect(atacs$dataset,"1|3"),"ctrl","lga")

atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
Idents(atacs)<-"predicted.id"
fcs<-FoldChange(atacs,subset.ident = "HSC",group.by = "group",ident.1 ="lga",ident.2 = "ctrl" ,assay = "lin_peaks")
peaks_5pct<-rownames(fcs)[fcs$pct.1>=0.05|fcs$pct.2>=0.05]

peaks_hsc_genes<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")

genes_5pct<-unique(peaks_hsc_genes[query_region%in%peaks_5pct]$gene_name)
length(genes_5pct) #20k

remove(atacs)

stims<-readRDS("outputs/27-multiomics_cbl_stim/cbls_stim_final.rds")
stims
head(stims[[]])

#by batch of 1000
library(parallel)


genes_list <- split(genes_5pct, ceiling(seq_along(genes_5pct)/1000))

pgl<-Reduce(rbind,mclapply(genes_list,function(genes){
  
    stims <- LinkPeaks(
    object = stims,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use = genes)
#13760 genes and 113542 peaks
    return(data.table(as.data.frame(Links(stims))))
    
    },mc.cores = 20))


fwrite(pgl,fp(out,"res_peak_genes_all_linkage.csv.gz"))


message("finished !")

