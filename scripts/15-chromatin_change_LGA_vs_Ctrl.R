#chromatin_change_LGA_vs_Ctrl

out<-"outputs/15-chromatin_change_LGA_vs_Ctrl"
dir.create(out)
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
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
DefaultAssay(atacs)<-"lin_peaks"

peaks_hsc_lga #rm sexual chromosome

peaks_hsc_lga_xy<-peaks_hsc_lga[!str_detect(peak,"chr[XY]")]
atacs@assays$lin_peaks

fwrite(peaks_hsc_lga_xy,fp(out,"differential_peaks_accessibility_lga_vs_ctrl_hsc_without_xy.csv.gz"))

ggplot(peaks_hsc_lga_xy,aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.001&abs(avg_log2FC)>0.25))+
  geom_point()+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#Motif enrichment
atacs@assays$lin_peaks@motifs<-readRDS("outputs/14-DMCs_atac_integr/atacs_cbps_lin_peaks_motif_object.rds")

atacs<-RegionStats(atacs,assay='lin_peaks',genome = BSgenome.Hsapiens.UCSC.hg38)

hsc_lga_tf<-FindMotifs(object = atacs,
           assay = "lin_peaks",
           features = peaks_hsc_lga_xy[p_val_adj<0.001&avg_log2FC>0.25]$peak )
hsc_lga_tf_down<-FindMotifs(object = atacs,
           assay = "lin_peaks",
           features = peaks_hsc_lga_xy[p_val_adj<0.001&avg_log2FC<(-0.25)]$peak )
hsc_lga_tf_down
hsc_lga_tf_all<-rbind(data.table(hsc_lga_tf)[,accessibility:="up"],data.table(hsc_lga_tf_down)[,accessibility:="down"])
hsc_lga_tf_all[order(pvalue)][1:20]

fwrite(hsc_lga_tf_all,fp(out,"motif_enrichment_in_peaks_up_and_down_lga_vs_ctrl_hsc.csv.gz"))



MotifPlot(
  object = atacs,
  motifs = head(hsc_lga_tf_all[enriched_in=="ctrl"]$motif,12),
  assay = 'lin_peaks')

#DMCs/DEGs enrichemnt





