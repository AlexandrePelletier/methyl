#chromatin_change_LGA_vs_Ctrl
source("scripts/utils/new_utils.R")
library(Seurat)
#renv::install("Signac")
library(Signac)

out<-"outputs/15-chromatin_change_LGA_vs_Ctrl"
dir.create(out)
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
atacs@assays$lin_peaks@motifs<-readRDS("outputs/14-DMCs_atac_integr/atacs_cbps_lin_peaks_motif_object.rds")

DefaultAssay(atacs)<-"lin_peaks"

#1) TF motif enrichment/ cells####
#Computing motif activities using ChromVar
#ChromVAR identifies motifs associated with variability in chromatin accessibility between cells. 
#See the chromVAR paper for a complete description of the method.
#renv::install("GreenleafLab/chromVAR")
library(BSgenome.Hsapiens.UCSC.hg38)
atacs <- RunChromVAR(
  object = atacs,assay = "lin_peaks",
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(atacs@assays$chromvar,fp(out,"cbps_atacs_lin_peaks_chromvar_assay.rds"))
DefaultAssay(atacs) <- 'chromvar'
Idents(atacs) <- 'predicted.id'

# look at the activity of EGR1
motifs<-res_motif_lin[motif.name%in%tfs]$motif
#mot<-res_motif_lin[motif.name%in%c("HES1","EGR3")]$motif

FeaturePlot(
  object = atacs,
  features = motifs,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  label=T,
  reduction="humap"
)
VlnPlot(atacs,motifs)
head(atacs[[]])
atacs$group<-ifelse(str_detect(atacs$dataset,"1|3"),"ctrl","lga")
differential.activity <- FindMarkers(
  object = atacs,
  subset.ident = "HSC",
  group.by = "group",
  ident.1 = 'lga',
  ident.2 = 'ctrl',
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
diffact<-data.table(differential.activity,keep.rownames = "motif")
diffact<-merge(diffact,unique(res_motif_lin[,.(motif,motif.name)]))
fwrite(diffact,fp(out,"differential_motif_activity_hsc_lga_vs_ctrl.csv"))
diffact[p_val_adj<0.001]$motif.name
head(diffact[order(p_val_adj)],100)

MotifPlot(
  object = atacs,
  motifs = head(rownames(differential.activity),10),
  assay = 'lin_peaks'
)

#2) diff accessi ctrl lga####
#Run 15A
peaks_hsc_lga<-fread(fp(out,"differential_peaks_accessibility_lga_vs_ctrl_hsc_logFC0.csv.gz"))

peaks_hsc_lga #rm sexual chromosome

peaks_hsc_lga_xy<-peaks_hsc_lga[!str_detect(peak,"chr[XY]")]

fwrite(peaks_hsc_lga_xy,fp(out,"differential_peaks_accessibility_lga_vs_ctrl_hsc_without_xy.csv.gz"))

ggplot(peaks_hsc_lga_xy,aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.001&abs(avg_log2FC)>0.25))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#Motif enrichment
#background : peaks test for DA (min.pct>= 0.05)
fcs<-FoldChange(atacs,subset.ident = "HSC",group.by = "group",ident.1 ="lga",ident.2 = "ctrl" ,assay = "lin_peaks")
peaks_5pct<-rownames(fcs)[fcs$pct.1>=0.05|fcs$pct.2>=0.05]
sum(fcs$pct.1==0|fcs$pct.2==0) #14k no express in LGA or ctrl
fcs[fcs$pct.1==0|fcs$pct.2==0,]
sum((fcs$pct.1==0|fcs$pct.2==0)&fcs$pct.2>0.05) #2
sum((fcs$pct.1==0|fcs$pct.2==0)&fcs$pct.1>0.05) #3

atacs<-RegionStats(atacs,assay='lin_peaks',genome = BSgenome.Hsapiens.UCSC.hg38)

hsc_lga_tf_up<-FindMotifs(object = atacs,
           assay = "lin_peaks",
           features = peaks_hsc_lga_xy[p_val_adj<0.001&avg_log2FC>0.25]$peak,
           background = peaks_5pct)
hsc_lga_tf_down<-FindMotifs(object = atacs,
           assay = "lin_peaks",
           features = peaks_hsc_lga_xy[p_val_adj<0.001&avg_log2FC<(-0.25)]$peak,
           background = peaks_5pct)
hsc_lga_tf_down
hsc_lga_tf_all<-rbind(data.table(hsc_lga_tf)[,accessibility:="up"],data.table(hsc_lga_tf_down)[,accessibility:="down"])
hsc_lga_tf_all[order(pvalue)][1:20]

fwrite(hsc_lga_tf_all,fp(out,"motif_enrichment_in_peaks_up_and_down_lga_vs_ctrl_hsc.csv.gz"))



MotifPlot(
  object = atacs,
  motifs = head(hsc_lga_tf_all[enriched_in=="ctrl"]$motif,12),
  assay = 'lin_peaks')

#DMCs/DEGs enrichemnt
#annot peaks
annotations<-readRDS("../atac/ref/gene_annotations_hg38_GRanges.rds")
Annotation(atacs) <- annotations

peaks_dt<-fread(fp(out,"cbps_atacs_peaks_by_lineage2.csv.gz"))

peaks_hsc_genes<-data.table(ClosestFeature(atacs,regions =peaks_dt[lineage=="HSC"]$peaks ))
fwrite(peaks_hsc_genes,fp(out,"peaks_hsc_genes_anno.csv.gz"))

#merge with res_degs
res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")
peaks_hsc_genes[,gene:=gene_name]
peaks_expr<-merge(peaks_hsc_genes,res_degs[lineage=="HSC"],by = "gene")
peaks_expr#63959 peaks asso to express genes
length(unique(peaks_expr$gene)) #10157 genes

#merge with res_cpgs
peaks_cpgs<-fread("outputs/14-DMCs_atac_integr/cpgs_in_lin_OCRs.csv.gz")
peaks_meth<-merge(peaks_cpgs,fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz"))
length(unique(peaks_meth$peaks)) #66196 peaks
length(unique(peaks_meth$cpg_id)) #274516 CpGs


#term list : peaks DMCs/DEGs
meth_degs_peaks_list<-list(CpGs=unique(peaks_meth$peaks),
                         DMCs=unique(peaks_meth[P.Value<0.001&abs(logFC)>25]$peaks),
                         expressed=unique(peaks_expr$query_region),
                         DEGs=unique(peaks_expr[padj<0.05&abs(log2FoldChange)>0.5]$query_region),
                         DMCs_DEGs=intersect(unique(peaks_meth[P.Value<0.001&abs(logFC)>25]$peaks),
                                             unique(peaks_expr[padj<0.05&abs(log2FoldChange)>0.5]$query_region)))
lapply(meth_degs_peaks_list, length)
lapply(meth_degs_peaks_list, head)

#query list : peaks up dn
da_peaks_list<-list(da_peaks=unique(peaks_hsc_lga_xy[p_val_adj<0.05&abs(avg_log2FC)>0.25]$peak),
                down_peaks=unique(peaks_hsc_lga_xy[p_val_adj<0.05&avg_log2FC<(-0.25)]$peak),
                up_peaks=unique(peaks_hsc_lga_xy[p_val_adj<0.05&avg_log2FC>0.25]$peak))
lapply(da_peaks_list, length)
lapply(da_peaks_list, head)

#background : peaks test for DA (min.pct>= 0.05)
fcs<-FoldChange(atacs,subset.ident = "HSC",group.by = "group",ident.1 ="lga",ident.2 = "ctrl" ,assay = "lin_peaks")
peaks_5pct<-rownames(fcs)[fcs$pct.1>=0.05|fcs$pct.2>=0.05]
length(peaks_5pct)#80k
background<-peaks_5pct
length(background) #79640/215117 peaks (37%)
head(background)

res_or_da_peaks_dmcs_degs<-OR3(da_peaks_list,
                    meth_degs_peaks_list,
                    background = background,
                    overlap_column = F) 

fwrite(res_or_da_peaks_dmcs_degs,fp(out,"res_da_peaks_enrichment_in_dmcs_degs_hsc_peaks.csv"))


#GRN 2.0


