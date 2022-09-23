out<-"outputs/24-JOBIM_22_poster"
dir.create(out)
source("scripts/utils/new_utils.R")

#n DA OCRs by ct
res_da<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_by_lin_logFC0.csv.gz")
res_da<-res_da[!str_detect(peak,"Y|X")]

table(res_da[p_val_adj<0.001&abs(avg_log2FC)>0.25]$lineage)
# Erythro-Mas         HSC    Lymphoid    MPP/LMPP 
#          83         250          43         120 

res_da[,features:=factor(ifelse(avg_log2FC>0,"up","down"))]
cts<-c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")
res_da[,lineage:=factor(lineage,levels = cts)]
p1<-ggplot(res_da[p_val_adj<0.001&abs(avg_log2FC)>0.25])+
  geom_bar(aes(x=lineage,fill=features),position = position_dodge2(width = 0.9, preserve = "single"))+
  scale_fill_manual(values = c("darkblue","cyan"))+
  scale_x_discrete(limits=cts,guide = guide_axis(angle=60))+
  scale_y_continuous(expand=c(0,0))+
  theme_minimal()+
  #theme(legend.position="bottom")+
  labs(x=element_blank())+ggtitle("Differential Accessibility")
#ggsave(fp(out,"n_DA_OCRs_by_lineage.pdf"))

#n DEGs by ct
res_de<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")

table(res_de[padj<0.01&abs(log2FoldChange)>0.25]$lineage)
     # B cell          DC Erythro-Mas         HSC    Lymphoid    MPP/LMPP     Myeloid      T cell 
     #     16           3          18         173           2          89          10           1 

res_de[,features:=factor(ifelse(log2FoldChange>0,"up","down"))]
cts<-c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")
res_de[,lineage:=factor(lineage,levels = cts)]
p2<-ggplot(res_de[padj<0.01&abs(log2FoldChange)>0.25&lineage%in%cts])+
  geom_bar(aes(x=lineage,fill=features),position = position_dodge2(width = 0.9, preserve = "single"))+
  scale_fill_manual(values = c("darkblue","cyan"))+
  scale_x_discrete(limits=cts,guide = guide_axis(angle=60))+
  scale_y_continuous(expand=c(0,0))+
  theme_minimal()+
  #theme(legend.position="bottom")+
  labs(x=element_blank())+ggtitle("Differential Expression")
ps<-p1+p2+plot_layout(guides = 'collect')
ggsave(fp(out,"n_DAR_DEGs_by_lineage.pdf"))

#diff methylation
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

p0<-ggplot(res_meth)+
  geom_point(aes(x=logFC,y=-log10(P.Value),col=P.Value<0.001&abs(logFC)>25),size=0.1)+
  scale_color_manual(values = c("grey","red"))+
  scale_x_continuous(limits = c(-100,100))+theme_minimal()+
  theme(legend.position="none")+
  labs(x="methylation difference (%)")+
  ggtitle("Differential Methylation")
p0+p1+p2+plot_layout(guides = 'collect',widths = c(3,2,2))
ggsave(fp(out,"DMCs_and_n_DAR_DEGs_by_lineage.pdf"))

#UMAPs scGANs batch effect removal
library(Seurat)
cbps<-readRDS("outputs/scGANs_batch_correction/cbps_filtered_with_scGANs_red_dims.rds")
#cbpsf<-subset(cbps,hto==T)
p1<-DimPlot(cbps,reduction = "umap_null",group.by = c("batch"),raster = T)
p2<-DimPlot(cbps,reduction = "zumap",group.by = "batch",raster = T)
pa<-p1|p2+plot_layout(guides = "collect")

p3<-DimPlot(cbps,reduction = "umap_null",group.by = c("lineage_hmap"),raster = T)
p4<-DimPlot(cbps,reduction = "zumap",group.by = "lineage_hmap",raster = T)
pb<-p3|p4+plot_layout(guides = "collect")

(p1|p2)/(p3|p4)+plot_layout(guides = "collect")&NoAxes()
ggsave(fp(out,"umaps_scGAN_batch_correction.pdf"),height = 7)
#scGANs RNA & ATACs integration


#CRE identif with scGANs


