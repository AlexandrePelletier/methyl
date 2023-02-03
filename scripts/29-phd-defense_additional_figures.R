out<-"outputs/29-phd_defense_additional_figures"
dir.create(out)

source("scripts/utils/new_utils.R")


#DMCs  enrichment in peak up / down
inp<-"outputs/15-chromatin_change_LGA_vs_Ctrl/"
res_or_da_peaks_dmcs_degs<-fread(fp(inp,"res_da_peaks_enrichment_in_dmcs_degs_hsc_peaks.csv"))

res_or_da_peaks_dmcs_degs[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
 res_or_da_peaks_dmcs_degs[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_da_peaks_dmcs_degs[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_da_peaks_dmcs_degs[,da_peaks:=paste0(query,"\n(n=",n.query,")")]
res_or_da_peaks_dmcs_degs[,candidat_peaks:=factor(candidat_peaks,levels =res_or_da_peaks_dmcs_degs[query=="down_peaks"][order(fold.enrichment)]$candidat_peaks )]

ggplot(res_or_da_peaks_dmcs_degs[term%in%c("CpGs","DMCs")&query!="da_peaks"])+geom_col(aes(y=log2(fold.enrichment),x=da_peaks,fill=padj_b),position = "dodge")+
  facet_wrap("candidat_peaks")+
  scale_fill_manual(values=c("grey","orange","red","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 8))
ggsave(fp(out,"Ah-barplot_dmcs_degs_peaks_enrichment_in_hsc_da_peaks.pdf"))


#DEGs enrichment in DMCs peaks, down peaks
res_or_deg_peaks_dmcs_dap<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/res_degs_peaks_enrichment_in_dmcs_da_hsc_peaks.csv")

res_or_deg_peaks_dmcs_dap[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_deg_peaks_dmcs_dap[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_deg_peaks_dmcs_dap[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_deg_peaks_dmcs_dap[,degs_peaks:=paste0(query,"\n(n=",n.query,")")]

res_or_deg_peaks_dmcs_dap<-res_or_deg_peaks_dmcs_dap[!term%in%c('peaks','CpGs',"DA_peaks")]
res_or_deg_peaks_dmcs_dap[,candidat_peaks:=factor(candidat_peaks,levels =res_or_deg_peaks_dmcs_dap[query=="downreg"][order(fold.enrichment)]$candidat_peaks[c(2,1,3,4)] )]

ggplot(res_or_deg_peaks_dmcs_dap[query.!="DEGs"])+geom_col(aes(y=log2(fold.enrichment),x=candidat_peaks,fill=padj_b),position = "dodge")+
  facet_wrap("degs_peaks")+
  scale_fill_manual(values=c("grey","orange","red","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))+theme_bw()


ggplot(res_or_deg_peaks_dmcs_dap[query.!="DEGs"])+
  geom_point(aes(y=-log10(pval),x=candidat_peaks,size=n.overlap ,col=fold.enrichment))+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  facet_wrap("degs_peaks")+
  scale_color_gradient(low = "grey",high="red")+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))+theme_minimal()

ggsave(fp(out,"Ah-barplot_dmcs_degs_peaks_enrichment_in_hsc_da_peaks.pdf"))


#DMCs/DEGs enrichment in DA CREs
inp<-"outputs/27-multiomics_cbl_stim/"
res_or_da_cres_dmcs_degs<-fread(fp(inp,"res_da_cres_enrichment_for_dmcs_degs2.csv"))

res_or_da_cres_dmcs_degs[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_da_cres_dmcs_degs[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_da_cres_dmcs_degs[,candidat_CREs:=paste0(term,"\n(n=",term.size,")")]
res_or_da_cres_dmcs_degs
res_or_da_cres_dmcs_degs[,da_CREs:=paste0(query,"\n(n=",n.query,")")]
res_or_da_cres_dmcs_degs[,candidat_CREs:=factor(candidat_CREs,levels =res_or_da_cres_dmcs_degs[query=="down_CREs"][order(fold.enrichment)]$candidat_CREs )]

ggplot(res_or_da_cres_dmcs_degs[query.=="down_CREs"])+geom_col(aes(y=log2(fold.enrichment),x=candidat_CREs,fill=padj_b),position = "dodge")+
  facet_wrap("da_CREs")+
  scale_fill_manual(values=c("grey","red","red4","darkred"))+
  scale_x_discrete(guide =guide_axis(angle = 66))+theme_minimal()+
  theme(axis.text.x = element_text(size = 8))
#ggsave(fp(out,"Ah-barplot_dmcs_degs_peaks_enrichment_in_hsc_da_peaks.pdf"))


#DMCs/ DA CREs enrichment in DEGs
res_or_degs_cres_dmcs_da<-fread("outputs/27-multiomics_cbl_stim/res_degs_cres_enrichment_for_dmcs_da_cres.csv")
res_or_degs_cres_dmcs_da[,padj_b:=ifelse(padj<0.001,"***",ifelse(padj<0.01,"**",ifelse(padj<0.05,'*','ns')))]
res_or_degs_cres_dmcs_da[,padj_b:=factor(padj_b,levels=c("ns","*","**","***"))]
res_or_degs_cres_dmcs_da[,candidat_peaks:=paste0(term,"\n(n=",term.size,")")]
res_or_degs_cres_dmcs_da[,degs_peaks:=paste0(query,"\n(n=",n.query,")")]

ggplot(res_or_degs_cres_dmcs_da[query.!="DEGs"])+
  geom_point(aes(y=-log10(pval),x=candidat_peaks,size=n.overlap ,col=fold.enrichment))+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  facet_wrap("degs_peaks")+
  scale_color_gradient(low = "grey",high="red")+
  scale_x_discrete(guide =guide_axis(angle = 66))+
  theme(axis.text.x = element_text(size = 8))+theme_minimal()


#LGA effect validation

#res degs/dap
res_degs2<-fread("outputs/27-multiomics_cbl_stim/res_hsc_stim2h_lga_vs_ctrl.csv")
res_degs0<-fread("outputs/27-multiomics_cbl_stim/res_hsc_no_stim_lga_vs_ctrl.csv")
res_degs0[p_val<0.001] #5


res_dap2<-fread("outputs/27-multiomics_cbl_stim/res_da_peaks_lga_vs_ctrl_hsc_2h_stim.csv")
res_dap0<-fread("outputs/27-multiomics_cbl_stim/res_da_peaks_lga_vs_ctrl_hsc_no_stim.csv")

res_dap2[p_val<0.0001] #55
res_dap0[p_val<0.0001] #55

res_degs0<-res_degs0[,feature:=gene][,-"gene"]
res_degs2<-res_degs2[,feature:=gene][,-"gene"]
res_dap0<-res_dap0[,feature:=peak][,-"peak"]
res_dap2<-res_dap2[,feature:=peak][,-"peak"]

res_dap0<-res_dap0[,closest_gene:=gene][,-"gene"]
res_dap2<-res_dap2[,closest_gene:=gene][,-"gene"]


res_all<-Reduce(function(...)rbind(...,fill=T),list(res_degs0[,condition:='no_stim'][,test:='DE'],
                                  res_degs2[,condition:='2h_stim'][,test:='DE'],
                                  res_dap0[,condition:='no_stim'][,test:='DA'],
                                  res_dap2[,condition:='2h_stim'][,test:='DA']))

egrn<-fread("outputs/13-GRN_integr/egr1_network_tf_target_interactions.csv")

genes_of_int<-union(egrn$target,egrn$tf)
ggplot(res_all,aes(y=-log10(p_val),x=avg_log2FC))+
  geom_point(aes(col=p_val<0.001))+
  facet_grid(condition~test)+
  scale_color_manual(values = c("grey","red"))+
  geom_text_repel(aes(label=ifelse(p_val<0.001&feature%in%genes_of_int,feature,ifelse(p_val<0.001&closest_gene%in%genes_of_int,closest_gene,""))),
                  
                  max.overlaps=1000)


#with sc multiome interactions
res_motifpeak_genes<-fread("outputs/27-multiomics_cbl_stim/res_peaks_genes_links_EGR1_KLF2_KLF4_motifs_anno.csv")

res_motifpeak_genesf<-res_motifpeak_genes[!is.na(motif)]
tfs<-c("EGR1","KLF2")
genes_egrncres<-unique(res_motifpeak_genesf[motif.name%in%tfs][,.(gene,motif.name,peak)])
genes_of_int<-unique(genes_egrncres$gene)
ggplot(res_all,aes(y=-log10(p_val),x=avg_log2FC))+
  geom_point(aes(col=p_val<0.001))+
  facet_grid(condition~test)+
  scale_color_manual(values = c("grey","red"))+
  geom_text_repel(aes(label=ifelse(p_val<0.001&feature%in%genes_of_int,feature,ifelse(p_val<0.001&closest_gene%in%genes_of_int,closest_gene,"")),
                      ),
                  max.overlaps=1000)



