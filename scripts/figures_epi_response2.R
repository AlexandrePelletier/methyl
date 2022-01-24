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
