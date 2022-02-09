
#TF signature enriched in DMC  ?
source("../methyl/scripts/utils/new_utils.R")
out<-here("../methyl/outputs/tf_signature_enrichment_in_DMC")
dir.create(out)
regulons_list<-readRDS("../singlecell/outputs/06-SCENIC/cbps0-8_clean/regulons_list.rds")

#in GeneScore
res_meth<-fread("../methyl/outputs/model14_without_iugr/2020-10-08_all_res_with_perm.csv")

res_meth[compa=="CF.LF"&gene=="HES1"]

#GSEA
library(fgsea)
?fgsea
res_methg<-unique(res_meth[order(compa,gene,pval)],by=c("compa","gene"))
res_methg[GeneScore>0,gs_rank:=rank(GeneScore),by="compa"]
res_methg[GeneScore<=0,gs_rank:=-rank(abs(GeneScore)),by="compa"]


genes_rank<-res_methg[compa=="CF.LF"]$gs_rank
names(genes_rank)<-res_methg[compa=="CF.LF"]$gene
res_gsea<-fgsea(pathways=regulons_list,
      stats=genes_rank,eps=0)

fwrite(res_gsea,f)
res_gsea[padj<0.001]
res_gsea[order(padj)]$pathway[1:10]
res_gsea[pathway=="STAT3e"]
res_gsea[pathway=="STAT3"]

fwrite(res_gsea,fp(out,"res_cflf_gsea_genescore_ranked.csv.gz"),sep=";")

#for each compa
#run tf1-

res_gsea_all<-fread("../methyl/outputs/tf_signature_enrichment_in_DMC/res_gsea_genescore_ranked.csv.gz")

res_gsea_all[,top10:=padj<=sort(padj)[10],by="compa"]

res_gsea_all[top10==T,.(pathway,compa,padj,NES)]
res_gsea_all[pathway=="STAT3",.(pathway,compa,padj,NES)]

#C vs L
#top 20 regulons by gsea in C vs L
res_gsea<-res_gsea_all[compa=="C.L"][,regulon:=pathway][,-"pathway"]
regulons_ord<-res_gsea[padj<=sort(padj)[15]][order(padj)]$regulon

ggplot(res_gsea[padj<=sort(padj)[15]])+geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+scale_x_discrete(limits=regulons_ord)

regulons_ord<-res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]][order(padj)]$regulon
ggplot(res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]])+
  geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+
  scale_x_discrete(limits=regulons_ord)

#=> TF epigen alteration by lineage
res_gsea<-res_gsea_all[compa=="C.L"][,regulon:=pathway][,-"pathway"]

regulons_lineage<-fread("../singlecell/outputs/06-SCENIC/cbps0-8_clean/tf_markers_lineage.csv.gz")

regulons_lineage_f<-regulons_lineage[p_val_adj<10^-30&avg_log2FC>0]
table(regulons_lineage_f$lineage)
res_gsea_lineage<-merge(res_gsea,regulons_lineage_f[,.(regulon,lineage,p_val_adj)])
res_gsea_lineage<-res_gsea_lineage[padj<0.001&p_val_adj<0.001]

ggplot(res_gsea_lineage[!str_detect(regulon,'e$')&lineage%in%c("HSC","MPP","Erythroid","Lymphoid","Myeloid")])+
  geom_col(aes(x=regulon,y=-log10(padj),fill=-log10(p_val_adj)))+facet_grid( .~ lineage, space = "free_x",scales="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 9),
        axis.title=element_text(size=12,face="bold"))
table(res_gsea_lineage[padj<0.001]$lineage)
res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')]

ggplot(res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')])+geom_bar(aes(x=lineage))
ggplot(res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')])+geom_boxplot(aes(x=lineage,y=-log10(padj)))



#CF vs LF
#top 20 regulons by gsea in C vs L
res_gsea<-res_gsea_all[compa=="CF.LF"][,regulon:=pathway][,-"pathway"]
regulons_ord<-res_gsea[padj<=sort(padj)[15]][order(padj)]$regulon

ggplot(res_gsea[padj<=sort(padj)[15]])+geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+scale_x_discrete(limits=regulons_ord)

regulons_ord<-res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]][order(padj)]$regulon
ggplot(res_gsea[!str_detect(regulon,'e$')][padj<=sort(padj)[15]])+
  geom_col(aes(x=regulon,y=-log10(padj),fill=NES))+
  scale_x_discrete(limits=regulons_ord)


#=> TF epigen alteration by lineage
res_gsea<-res_gsea_all[compa=="CF.LF"][,regulon:=pathway][,-"pathway"]

res_gsea_lineage<-merge(res_gsea,regulons_lineage_f[,.(regulon,lineage,p_val_adj)])
res_gsea_lineage<-res_gsea_lineage[padj<0.001&p_val_adj<0.001]

res_gsea_lineage[,n.enriched:=length(tr(leadingEdge,sep = "|")),by="regulon"]
res_gsea_lineage[,pct.enriched:=n.enriched/size]

ggplot(res_gsea_lineage[!str_detect(regulon,'e$')&lineage%in%c("HSC","MPP","Erythroid","Lymphoid","Myeloid")])+
  geom_point(aes(x=regulon,y=pct.enriched,size=size,col=-log10(p_val_adj+1e-300)))+facet_grid( .~ lineage, space = "free_x",scales="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 9),
        axis.title=element_text(size=12,face="bold"))

table(res_gsea_lineage[padj<0.001]$lineage)
res_gsea_lineage[padj<0.001&!str_detect(regulon,'e$')]

res_gsea[regulon=="JUN"]
VlnPlot(cbps,c("JUN"),group.by = "hto",split.by ="group" ,idents = "HSC",pt.size = 0)
dt<-data.table(JUN=cbps@assays$TF_AUC@counts["JUN",],
           cbps@meta.data)

plot(density(dt$JUN))
abline(v=0.5)
dt[,pct.JUN:=sum(JUN>0.5)/.N,by=c("lineage","sample","hto")]
ggplot(unique(dt,by=c("sample","hto")))+geom_boxplot(aes(x=hto,y=pct.JUN,fill=group))+facet_wrap("lineage")


#over representation test
res_meth<-fread("../methyl/outputs/model14_without_iugr/2020-10-08_all_res_with_perm.csv")
res_meth[,gene_score_scaled:=scale(GeneScore,center = F),by="compa"]
res_methg<-unique(res_meth[order(compa,gene,pval)],by=c("compa","gene"))
ggplot(res_methg)+geom_density(aes(x=gene_score_scaled,fill=compa))+facet_wrap("compa")
ggplot(res_methg)+geom_density(aes(x=gene_score_scaled,fill=compa))+facet_wrap("compa")+geom_vline(xintercept = 1.25)

ggplot(res_methg[abs(gene_score_scaled)>1.25])+geom_bar(aes(x=compa,fill=compa))
res_methg[,gene_altered:=abs(gene_score_scaled)>1.25]

#CFLF
genes_altered<-res_methg[gene_altered==T&compa=='CF.LF']$gene
size_universe<-length(res_methg[compa=='CF.LF']$gene)
res_or<-data.table(regulon=names(regulons_list))
res_or[,pval:=sapply(regulon,function(tf){
  return(over_repr_test(genes_altered,regulons_list[[tf]],size_universe =size_universe ))
  })]

res_or[,padj:=p.adjust(pval,method = 'BH')]
res_or[padj<0.05]

#for all

res_or_all<-Reduce(rbind,lapply(unique(res_methg$compa),function(comp){
  print(comp)
  res_or<-data.table(regulon=names(regulons_list),regulon.size=sapply(regulons_list,length),compa=comp)
  genes_altered<-res_methg[gene_altered==T&compa==comp]$gene
  res_or[,n.genes.altered:=length(genes_altered)]
  res_or[,n.enriched:=sum(genes_altered%in%regulons_list[[regulon]]),by="regulon"]
  res_or[,genes.enriched:=paste(genes_altered[genes_altered%in%regulons_list[[regulon]]],collapse="|"),by="regulon"]
  res_or[,pct.enriched:=n.enriched/regulon.size]
  
  size_universe<-length(res_methg[compa==comp]$gene)
  
  res_or[,pval:=phyper(q=n.enriched-1, 
       m=n.genes.altered, 
       n=size_universe-n.genes.altered, 
       k=regulon.size, 
       lower.tail=FALSE),
       by="regulon"]
  
  res_or[,padj:=p.adjust(pval,method = 'BH')]
  
  return(res_or)
  })
  )

res_or_all[padj<0.05]
ggplot(res_or_all[padj<0.05])+geom_bar(aes(x=compa))
ggplot(res_or_all[padj<1e-10])+geom_bar(aes(x=compa))

res_or_all[padj<1e-10&compa=="CF.LF"]

ggplot(res_or_all[padj<0.001&!str_detect(regulon,"e$")])+geom_bar(aes(x=compa))
res_or_all[padj<0.001&compa=="CF.LF"&!str_detect(regulon,"e$")]$regulon


ggplot(res_or_all[padj<0.001&compa=="CF.LF"&!str_detect(regulon,"e$")])+
  geom_point(aes(x=regulon,y=pct.enriched,size=regulon.size,col=-log10(padj)))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 9),
        axis.title=element_text(size=12,face="bold"))
ggplot(res_or_all[padj<0.001&compa=="C.L"&!str_detect(regulon,"e$")])+
  geom_point(aes(x=regulon,y=pct.enriched,size=regulon.size,col=-log10(padj)))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 9),
        axis.title=element_text(size=12,face="bold"))

fwrite(res_or_all,fp(out,"res_over_repr_by_compa_genescore_scaled_thr1.25.csv.gz"),sep=";")
res_or_all<-fread(fp(out,"res_over_repr_by_compa_genescore_scaled_thr1.25.csv.gz"),sep=";")

#pca on the cpgs determining TF

genes_determining_tf<-Reduce(union,regulons_list[unique(res_or_all[padj<0.001&!str_detect(regulon,"e$")]$regulon)])

res_meth[,tf_determining_genes:=gene%in%genes_determining_tf]

res_meth[,tf_determining_cpgs:=tf_determining_genes&pval<0.001&abs(meth.change)>25]

cpgs_determining_tf<-unique(res_meth[tf_determining_cpgs==T]$locisID)
length(cpgs_determining_tf)
meth<-fread("../methyl/datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
meth_deter<-meth[locisID%in%cpgs_determining_tf] #1413 cpgs

pca<-prcomp(t(as.matrix(meth_deter[,.SD,.SDcols=mtd[Group_name%in%c("C","L")]$sample])))

pca$x
vars_pcs<-GetVarPCs(pca)
round(vars_pcs*100)

res_pca_deter<-as.data.table(pca$x,keep.rownames = "sample")

mtd<-fread("../methyl/datasets/cd34/cleaned_batch_CD34_250121.csv")

res_pca_deter<-merge(res_pca_deter,mtd,by="sample")

ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC2,col=Group_name))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC2,col=Group_Sex))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC3,col=Group_Sex))
ggplot(res_pca_deter)+geom_point(aes(x=PC1,y=PC3,col=Group_Sex))

ggplot(res_pca_deter)+geom_boxplot(aes(x=Group_Sex,y=PC1,fill=Group_Sex))
