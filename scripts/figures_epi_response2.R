library(data.table)
source("scripts/utils/new_utils.R")
out<-"outputs/figures_epi_response2"
dir.create(out)

degs<-fread("outputs/08-HTO_signature/by_lineage/res_pseudobulk_DESeq2_3replicates.csv.gz")
degs<-degs[padj<0.05&abs(log2FoldChange)>0.5]
table(degs$lineage)
degs[,ndegs:=.N,by="lineage"]
degs_pro<-degs[ndegs>60& lineage!="DC"]
degs_pro[,common_degs:=.N==5,by="gene"]
unique(degs_pro[common_degs==T]$gene)

res_kegg_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_kegg_dn.csv")
res_kegg_dn<-res_kegg_dn[p.adjust<0.05]

res_go_dn<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_go_bp_dn.csv.gz")
res_go_dn[order(p.adjust)][p.adjust<0.05]$Description
res_go_dn<-res_go_dn[p.adjust<0.05]
terms_of_int<-rbind(res_go_dn[,source:="GO:BP"],res_kegg_dn[,source:="KEGG"])

terms_of_int[,intersection_size:=Count]
terms_of_int[,GeneRatio:=as.numeric(str_extract(GeneRatio,"[0-9]+"))/as.numeric(str_extract(GeneRatio,"[0-9]+$"))]

terms_of_int[,top20:=rank(p.adjust)<=20,by="source"]
terms_of_int$Description <- reorder(terms_of_int$Description, terms_of_int$p.adjust)
ggplot(terms_of_int[top20==T])+
  geom_point(aes(x=Description,y=-log10(p.adjust),size=Count,col=GeneRatio))+
  scale_x_discrete(guide = guide_axis(angle = 80))+facet_grid(~source,scales = "free",space = "free")

ggsave(fp(out,"dotplot_go_kegg_enrichment_HSC_LGA_down_top20_padj0.05.pdf"),width = 12,height = 8)



#figure ATAC
#
#a) figure UMAP predicted id
library(Seurat)
library(Signac)
inp<-"outputs/14-DMCs_atac_integr/"
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")

DimPlot(atacs, group.by = "predicted.id", label = TRUE,reduction = "humap")
ggsave(fp(out,"Aa-cbps_atac_umap_predicted_lineage.pdf"))


#b) Lineage spe TF
res_motif_lin<-fread(fp(inp,"enriched_motifs_by_lineage_specific_peaks.csv"))
#rm redundant motif
res_motif_lin[,top15:=rank(pvalue)<=15,by="lineage"]
res_motif_lin[,pvalue:=pvalue+10^-316]

res_motif_lin[,]
ggplot(res_motif_lin[pvalue<10^-50&lineage!="MPP/LMPP"])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  facet_grid(~lineage,scales = "free_x",space="free_x")+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))
ggsave(fp(out,"Ab-lineage_spe_tf_mlog10pval50.pdf"))

ggplot(res_motif_lin[pvalue<10^-50&lineage!="MPP/LMPP"])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  facet_wrap("lineage",scales = "free_x",ncol = 2)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60))+
  theme(axis.text=element_text(size=8))
ggsave(fp(out,"Ab-lineage_spe_tf_mlog10pval50_v2.pdf"))
fwrite(res_motif_lin,fp(out,"enriched_motifs_by_lineage_specific_peaks.csv"))

#c) DMCs enrichment in HSC peaks
res_dmcs_peaks_enr<-fread(fp(inp,"res_dmcs_peaks_enrichment_for_lineage_specific_peaks.csv"))
res_dmcs_peaks_enr[,percent.observed:=precision*100]
ggplot(res_dmcs_peaks_enr)+
  geom_point(aes(x=query,size=n.overlap,col=-log10(padj),y=percent.observed))+
  scale_color_gradient(low = "grey",high = "red")+
    scale_y_continuous(expand = c(0,1))+scale_x_discrete(guide = guide_axis(angle = 0))+
    theme(axis.title.x =element_blank())

ggsave(fp(out,"Ac-dmcs_enrichment_in_lineage_spe_peaks.pdf"))


#d) TF DMCs HSC peaks vs CpG peaks
res_motif_hsc_dmcs<-fread(fp(inp,"res_tf_motif_enrichment_in_DMCs_containing_vs_CpGs_containing_hsc_peaks.csv"))
res_motif_hsc_dmcs[,pvalue:=pvalue+10^-316]
ggplot(res_motif_hsc_dmcs[1:15])+
  geom_point(aes(x=motif.name,size=fold.enrichment,col=-log10(pvalue),y=percent.observed))+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100),expand = c(0,0))+scale_x_discrete(guide = guide_axis(angle = 60),
                                                                       limits=res_motif_hsc_dmcs[1:15][order(pvalue,-fold.enrichment)]$motif.name)+
  scale_size(limits = c(1,2.7))+
  theme(axis.text=element_text(size=9))

ggsave(fp(out,"Ad-dotplot_tf_motif_dmcs_vs_cpg_hsc_peaks_top15.pdf"))


#e)valid trs alteration méca
res_or_degs_dmcs_tf<-fread(fp(inp,"res_degs_enrichment_in_dmcs_tf_hsc_peaks.csv"))

res_or_degs_dmcs_tf[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_degs_dmcs_tf[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_degs_dmcs_tf[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_degs_dmcs_tf[,degs_peaks:=paste0(query,"\n(n=",n.query,")")]

ggplot(res_or_degs_dmcs_tf)+geom_col(aes(y=log2(fold.enrichment),x=candidat_peaks,fill=padj_b),position = "dodge")+
  facet_wrap("degs_peaks")+
  scale_fill_manual(values=c("grey","orange","red","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))
ggsave(fp(out,"Ae-barplot_dmcs_tf_peaks_enrichment_in_degs_peaks.pdf"))

#SUPP
#distrib predic lineage
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

ggplot(atacs@meta.data)+geom_bar(aes(x=predicted.id))

#can remove badly predict lineage ?
ggplot(atacs@meta.data)+geom_bar(aes(x=predicted.id))


