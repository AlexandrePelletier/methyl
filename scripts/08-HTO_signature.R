library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("scripts/utils/new_utils.R")
out<-"outputs/08-HTO_signature"
dir.create(out)

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

#cellcycle activation ####
mtd<-data.table(cbps@meta.data,keep.rownames = "bc")
mtd$Phase<-factor(mtd$Phase,levels = c("G1","S","G2M"))
mtd$lineage_hmap<-factor(mtd$lineage_hmap,levels = c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid"))

#G2M score
#for dupli
mtdf<-mtd[lineage_hmap%in%c("HSC",'MPP/LMPP')]
ggplot(mtdf[sample%in%c("ctrlM555","ctrlM518","ctrlM537")])+geom_boxplot(aes(x=sample,y=G2M.Score,fill=hto),outlier.shape = NA)+
  coord_cartesian(ylim = c(-0.2,0.12))

#for all
p<-ggplot(mtdf)+geom_boxplot(aes(x=hto,y=G2M.Score,fill=hto,group=sample_hto),outlier.shape = NA)+
  coord_cartesian(ylim = c(-0.2,0.12))
p+labs( x="Sample", y = "G2/M Score",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

#pct G2M
#calc pct global
mtd[,n.sample:=.N,.(sample_hto)]
mtd[,pct.phase:=.N/n.sample,by=.(sample_hto,Phase)]
mts<-unique(mtd,by=c("sample_hto","Phase"))

#calc pct by lineage
mtd[,n.sample.lin:=.N,.(sample_hto,lineage_hmap)]
mtd[,pct.phase.lin:=.N/n.sample.lin,by=.(sample_hto,lineage_hmap,Phase)]
mtsl<-unique(mtd,by=c("sample_hto","lineage_hmap","Phase"))

#plot 
#duplicates barplot
#global
mtsd<-mts[sample%in%c("ctrlM555","ctrlM518","ctrlM537")]
mtsd[,mean.pct:=mean(pct.phase),by=.(Phase,hto)]
ggplot(mtsd)+geom_col(aes(x=Phase,y=mean.pct,fill=hto),position = "dodge")

#global by sample for only G2M
p1<-ggplot(mtsd[Phase=="G2M"])+geom_col(aes(x=sample,y=pct.phase*100,fill=hto),position = "dodge")
p1<-p1+labs( x="Sample", y = "HSPC in G2/M Phase (%)",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

pa<-ggplot(mtsd[Phase=="G2M"])+geom_col(aes(y=pct.phase*100,x=hto,fill=hto),position = "dodge")
pa<-pa+labs( x="Sample", y = "HSPC in G2/M Phase (%)",x="HTO labelling",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

#only in HSC/MPP cells for only G2M
mtdf<-mtd[lineage_hmap%in%c("HSC",'MPP/LMPP')]
mtdf[,n.sample:=.N,.(sample_hto)]
mtdf[,pct.phase:=.N/n.sample,by=.(sample_hto,Phase)]
mtsf<-unique(mtdf,by=c("sample_hto","Phase"))
mtsfd<-mtsf[sample%in%c("ctrlM555","ctrlM518","ctrlM537")]
p2<-ggplot(mtsfd[Phase=="G2M"])+geom_col(aes(x=sample,y=pct.phase*100,fill=hto),position = "dodge")
p2<-p2+labs( x="Sample", y = "HSC/MPP cells in G2/M Phase (%)",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

pb<-ggplot(mtsfd[Phase=="G2M"])+geom_col(aes(x=hto,y=pct.phase*100,fill=hto),position = "dodge")
pb<-pb+labs( x="HTO labelling", y = "HSC/MPP cells in G2/M Phase (%)",fill="HTO labelling")+
   theme_classic()  +
   scale_fill_manual(values=c('#999999','#E69F00'))


p1+p2
ggsave(fp(out,"cellcycle_activation.pdf"))

pa+pb
ggsave(fp(out,"cellcycle_activation2.pdf"))

ggplot(mtd[sample%in%c("ctrlM555","ctrlM518","ctrlM537")])+
  geom_bar(aes(x=hto,fill=Phase),position = "fill")+
  facet_wrap('sample')

px<-ggplot(mtsfd[Phase!="G1"])+
  geom_col(aes(x=hto,y=pct.phase*100,fill=Phase),position = "stack")+
  facet_wrap('sample')
px<-px+labs( x="HTO labelling", y = "HSC/MPP in Cell Cycle Phase (%)",fill="Cell Cycle Phase")+
   theme_minimal() 

ggsave(fp(out,"cellcycle_activation3.pdf"))

mtdf2<-mtd[lineage_hmap%in%c("HSC",'MPP/LMPP')&sample%in%c("ctrlM555","ctrlM518","ctrlM537")]
mtdf2[,n.cells:=.N,.(hto)]
mtdf2[,pct.phase:=.N/n.cells,by=.(Phase,hto)]
mtsf2<-unique(mtdf2,by=c("Phase",'hto'))

py<-ggplot(mtsf2[Phase!="G1"])+
  geom_col(aes(x=hto,y=pct.phase*100,fill=Phase),position = "stack")
py+labs( x="HTO labelling", y = "HSC/MPP in Cell Cycle Phase (%)",fill="Cell Cycle Phase")+
   theme_minimal() 

ggsave(fp(out,"cellcycle_activation4.pdf"))


#all samples boxplot
#global by sample for only G2M
p<-ggplot(mts[Phase=="G2M"])+geom_boxplot(aes(x=hto,y=pct.phase*100,fill=hto))
p+labs( x="Sample", y = "Cells in G2/M Phase (%)",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))


#only in HSC/MPP cells for only G2M
mtdf<-mtd[lineage_hmap%in%c("HSC",'MPP/LMPP')]
mtdf[,n.sample:=.N,.(sample_hto)]
mtdf[,pct.phase:=.N/n.sample,by=.(sample_hto,Phase)]
mtsf<-unique(mtdf,by=c("sample_hto","Phase"))
p<-ggplot(mtsf[Phase=="G2M"])+geom_boxplot(aes(x=hto,y=pct.phase*100,fill=hto))
p+labs( x="Sample", y = "HSC/MPP cells in G2/M Phase (%)",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

#by group
p<-ggplot(mtsf[Phase=="G2M"])+geom_boxplot(aes(x=hto,y=pct.phase*100,fill=group))
p+labs( x="HTO labelling", y = "HSC/MPP cells in G2/M Phase (%)",fill="Group")+
   theme_classic()

#by lineage
mtsld<-rbind(mtsld,data.table(Phase="G2M",differentiated=F,sample="ctrlM537",hto=F,lineage_hmap='HSC',pct.phase.lin=0),fill=T)
p<-ggplot(mtsld[Phase=="G2M"&differentiated==F&lineage_hmap!="LT-HSC"],
          aes(x=sample,y=pct.phase.lin*100,fill=hto))+
  geom_col(position = "dodge")+
  facet_wrap("lineage_hmap")+scale_fill_discrete(drop=FALSE)

p+labs( x="Sample", y = "cells in G2/M Phase (%)",fill="HTO labelling")+
   theme_classic() +
   scale_fill_manual(values=c('#999999','#E69F00'))

p<- ggplot(mtsld[Phase=="G2M"&differentiated==F&lineage_hmap!="LT-HSC"],
           aes(x=lineage_hmap, y=mean.pct, fill=hto)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean.pct-se.pct, ymax=mean.pct+se.pct), width=.2,
                 position=position_dodge(.9)) 

p


#all samples boxplot
ggplot(mts)+geom_boxplot(aes(x=Phase,y=pct.phase,fill=hto))

ggplot(mtsl[differentiated==F])+geom_boxplot(aes(x=Phase,y=pct.phase,fill=hto))+
  facet_wrap("lineage_hmap")


#0) signature : same sample hto vs not####
#a) signature all cbps
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(DESeq2)
source("../methyl/scripts/utils/new_utils.R")


cbps_dup<-subset(cbps,sample%in%c("ctrlM555","ctrlM518","ctrlM537"))
table(cbps_dup$hto,cbps_dup$sample)
  #      ctrlM518 ctrlM537 ctrlM555
  # FALSE     2005     2795     1976
  # TRUE       672      467      610

#get mtd of interest
mtd<-data.table(cbps_dup@meta.data,keep.rownames = "bc")
#cell cycle activation
ggplot(mtd[differentiated==F])+geom_bar(aes(x=hto,fill=Phase),position = "fill")+
  facet_grid(sample~lineage_hmap)

ggplot(mtd[lineage_hmap%in%c("HSC","MPP/LMPP")])+geom_bar(aes(x=hto,fill=Phase),position = "fill")+
  facet_grid(lineage_hmap~sample)+
  coord_cartesian(ylim = c(0,0.5))+
  theme_minimal()+scale_fill_manual(values = c("white","orange","blue"))

mtd[,pct.cc:=sum(Phase%in%c("G2M","S"))/.N,.(lineage_hmap,sample,hto)]

ggplot(unique(mtd[lineage_hmap%in%c("HSC","MPP/LMPP","Lymphoid")],by=c("lineage_hmap","sample","hto")))+
  geom_boxplot(aes(x=hto,y=pct.cc,fill=hto),position = "dodge")+
  facet_wrap("lineage_hmap",scales = "free")

mts<-unique(mtd,by=c("sample","orig.ident"))


#DEGs
#get counts and filter genes lowly express
counts<-as.matrix(cbps_dup@assays$RNA@counts)
dim(counts) 

counts <- counts[rowSums(counts > 0) >= 100, ] 
dim(counts) #13203 genes

# Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample_hto, fun = "sum")) 


#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ hto )

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds,name =  "htoTRUE",
                alpha = 0.05)


res_dt<-data.table(as.data.frame(res),keep.rownames="gene")

res_dt[is.na(padj),padj:=1]
res_dt[padj<0.05&abs(log2FoldChange)>0.5] #1518
fwrite(res_dt,fp(out,"res_pseudobulk_DESeq2_3replicates.csv"),sep=";")


#pathway enrichment of the signature
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

res_dt<-fread("outputs/08-HTO_signature/res_pseudobulk_DESeq2_3replicates.csv")
res_dt[padj<0.05&abs(log2FoldChange)>0.5] #1518
res_dt[padj<0.05&log2FoldChange>0.5] #1075
1518-1075 #443


#gseGO
res_dt[,deg_score:=log2FoldChange*-log10(padj)]
#res_dt[,deg_score_rk:=rank(deg_score)]#du plus surexpr au plus downreg

res_dt1<-merge(res_dt,
              data.table(bitr(res_dt$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F))[,gene:=SYMBOL])[!is.na(ENTREZID)]
genelist<-res_dt1[order(-deg_score)]$deg_score
names(genelist)<-res_dt1[order(-deg_score)]$ENTREZID

res_gsea_go<- gseGO(geneList     =genelist , 
                    ont="BP",
                    exponent = 1,
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
saveRDS(res_gsea_go,fp(out,"res_gsea_go_bp.rds"))
res_gsea_go<-readRDS(fp(out,"res_gsea_go_bp.rds"))

dotplot(res_gsea_go,showCategory=50)

gsea_go_dt<-data.table(as.data.frame(res_gsea_go))
gsea_go_dt[,sens:=ifelse(NES>0,"up","dn")]
gsea_go_dt[,n.enriched:=length(tr(core_enrichment)),"ID"]
gsea_go_dt[p.adjust<0.05]
gsea_go_dt[p.adjust<0.05&NES>0]$Description

fwrite(gsea_go_dt,fp(out,"res_gsea_go_bp.csv.gz"))
gsea_go_dt<-fread(fp(out,"res_gsea_go_bp.csv.gz"))

#with rank degscore
res_dt[log2FoldChange>0,deg_score:=rank(log2FoldChange*-log10(padj))]
res_dt[log2FoldChange<=0,deg_score:=-rank(abs(log2FoldChange)*-log10(padj))]

#res_dt[,deg_score_rk:=rank(deg_score)]#du plus surexpr au plus downreg

res_dt1<-merge(res_dt,
              data.table(bitr(res_dt$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F))[,gene:=SYMBOL])[!is.na(ENTREZID)]
genelist<-res_dt1[order(-deg_score)]$deg_score
names(genelist)<-res_dt1[order(-deg_score)]$ENTREZID

res_gsea_go<- gseGO(geneList     =genelist , 
                    ont="BP",
                    exponent = 1,
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
saveRDS(res_gsea_go,fp(out,"res_gsea_go_bp_rank.rds"))
dotplot(res_gsea_go,showCategory=50)

gsea_go_dt<-data.table(as.data.frame(res_gsea_go))
gsea_go_dt[,sens:=ifelse(NES>0,"up","dn")]
gsea_go_dt[,n.enriched:=length(tr(core_enrichment)),"ID"]
gsea_go_dt[p.adjust<0.05]
gsea_go_dt[p.adjust<0.05&NES>0]$Description

fwrite(gsea_go_dt,fp(out,"res_gsea_go_bp_rank.csv.gz"))

#kegg
res_gsea_kegg<- gseKEGG(geneList     =genelist , 
                    organism = "hsa",
                    exponent = 1,
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0)

dotplot(res_gsea_kegg,showCategory=50)

res_gsea_kegg_dt<-data.table(as.data.frame(res_gsea_kegg))
res_gsea_kegg_dt[,sens:=ifelse(NES>0,"up","dn")]
res_gsea_kegg_dt[,n.enriched:=length(tr(core_enrichment)),"ID"]
res_gsea_kegg_dt[p.adjust<0.05]
res_gsea_kegg_dt[p.adjust<0.05&NES>0]$Description

saveRDS(res_gsea_kegg,fp(out,"res_gsea_kegg.rds"))
fwrite(res_gsea_kegg_dt,fp(out,"res_gsea_kegg.csv.gz"))


#enrich
possible_genes<-rownames(readRDS('outputs/06-integr_singlecell_cbps/cbps_light.rds'))
# -kegg
res_kegg<-enrichKEGG(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
saveRDS(res_kegg,fp(out,"res_hto_signature_kegg.rds"))
res_kegg_dt<-data.table(as.data.frame(res_kegg))
res_kegg_dt[p.adjust<0.1]#39 : MAPK, TNF, microrna in cancer
fwrite(res_kegg_dt,fp(out,"res_hto_signature_kegg.csv"))


res_kegg_up<-enrichKEGG(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
saveRDS(res_kegg_up,fp(out,"res_hto_signature_kegg_up.rds"))
res_kegg_up<-readRDS("outputs/08-HTO_signature/res_hto_signature_kegg_up.rds")

res_kegg_dt<-data.table(as.data.frame(res_kegg_up))
res_kegg_dt[p.adjust<0.1]# 38 : MAPK, TNF, IL17, FOXO, TGFB, Hippo..
fwrite(res_kegg_dt,fp(out,"res_hto_signature_kegg_up.csv"))


res_kegg_dn<-enrichKEGG(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
saveRDS(res_kegg_dn,fp(out,"res_hto_signature_kegg_dn.rds"))
res_kegg_dn<-readRDS("outputs/08-HTO_signature/res_hto_signature_kegg_dn.rds")

res_kegg_dt<-data.table(as.data.frame(res_kegg_dn))
res_kegg_dt[p.adjust<0.1]#11 : Biosynthesis of unsaturated fatty acids,Oxidative phosphorylation 
fwrite(res_kegg_dt,fp(out,"res_hto_signature_kegg_dn.csv"))

#gsekegg
fcs<-res_dt$log2FoldChange*-log10(res_dt$pvalue)
names(fcs)<-possible_genes
fcs<-sort(fcs,decreasing = T)
head(fcs)

genes.df<-bitr(names(fcs),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
head(genes.df)

fcs<-fcs[genes.df$SYMBOL]
names(fcs)<-genes.df$ENTREZID
head(fcs)
res_gsek<-gseKEGG(geneList     = fcs,
                      exponent=0,
                  eps=0,
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 0.05,
                        verbose = FALSE)

saveRDS(res_gsek,fp(out,"res_hto_signature_gsea_kegg.rds"))
res_gsek<-readRDS(fp(out,"res_hto_signature_gsea_kegg.rds"))

res_kegg_dt<-data.table(as.data.frame(res_gsek))
res_kegg_dt
fwrite(res_kegg_dt,fp(out,"res_hto_signature_kegg.csv"))


# -go
res_go_mf<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "MF",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)


saveRDS(res_go_mf,fp(out,"res_hto_signature_go_mf.rds"))
res_go_mf_dt<-data.table(as.data.frame(res_go_mf))
fwrite(res_go_mf_dt,fp(out,"res_hto_signature_go_mf.csv"))

res_go_bp<-enrichGO(bitr(res_dt[padj<0.05&abs(log2FoldChange)>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

saveRDS(res_go_bp,fp(out,"res_hto_signature_go_bp.rds"))
res_go_bp_dt<-data.table(as.data.frame(res_go_bp))
fwrite(res_go_bp_dt,fp(out,"res_hto_signature_go_bp.csv"))

#go bp up
res_go_bp<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange>0.5]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
saveRDS(res_go_bp,fp(out,"res_hto_signature_go_bp_up.rds"))
res_go_bp<-readRDS(fp(out,"res_hto_signature_go_bp_up.rds"))

res_go_bp_dt<-data.table(as.data.frame(res_go_bp))
res_go_bp_dt[p.adjust<0.1]#312
fwrite(res_go_bp_dt,fp(out,"res_hto_signature_go_bp_up.csv"))

#go bp dn
res_go_bp<-enrichGO(bitr(res_dt[padj<0.05&log2FoldChange<(-0.5)]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,
                 universe =bitr(possible_genes,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)
saveRDS(res_go_bp,fp(out,"res_hto_signature_go_bp_dn.rds"))
res_go_bp_dt<-data.table(as.data.frame(res_go_bp))
res_go_bp_dt[p.adjust<0.1]#103
fwrite(res_go_bp_dt,fp(out,"res_hto_signature_go_bp_dn.csv"))

#b) signature by lineage####
out1<-fp(out,"by_lineage")
dir.create(out1)
Idents(cbps_dup)<-"lineage_hmap"
res_lin_dup<-Reduce(rbind,lapply(levels(cbps_dup),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_dup,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample","orig.ident"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample_hto, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                                 design = ~ hto)
  
  dds <- DESeq(dds)
  
  res <- results(dds,name = "htoTRUE",alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

fwrite(res_lin_dup,fp(out1,"res_pseudobulk_DESeq2_3replicates.csv.gz"))



#sign enrich by lineage #[to redo]
res_lin_dup<-fread("outputs/08-HTO_signature/by_lineage/res_pseudobulk_DESeq2_3replicates.csv.gz")

# -kegg
res_kegg_lin<-lapply(unique(res_lin_dup$lineage), function(lin)enrichKEGG(bitr(res_lin_dup[padj<0.05&abs(log2FoldChange)>0.6&lineage==lin]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05))
names(res_kegg_lin)<-unique(res_lin_dup$lineage)
saveRDS(res_kegg_lin,fp(out1,"res_hto_signature_kegg_by_lineage.rds"))

res_kegg_lin_dt<-Reduce(rbind,lapply(names(res_kegg_lin)[!sapply(res_kegg_lin, is.null)],function(lin)data.table(as.data.frame(res_kegg_lin[[lin]]))[,lineage:=lin]))
res_kegg_lin_dt#yes !
fwrite(res_kegg_lin_dt,fp(out1,"res_hto_signature_kegg_by_lineage.csv"))


# -go
res_go_mf_lin<-lapply(unique(res_lin_dup$lineage), function(lin)enrichGO(bitr(res_lin_dup[padj<0.05&abs(log2FoldChange)>0.6&lineage==lin]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "MF",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05))
names(res_go_mf_lin)<-unique(res_lin_dup$lineage)

saveRDS(res_go_mf_lin,fp(out1,"res_hto_signature_go_mf_by_lineage.rds"))
res_go_mf_lin_dt<-Reduce(rbind,lapply(names(res_go_mf_lin)[!sapply(res_go_mf_lin, is.null)],function(lin)data.table(as.data.frame(res_go_mf_lin[[lin]]))[,lineage:=lin]))
res_go_mf_lin_dt
table(res_go_mf_lin_dt$lineage)
  # B cell          DC Erythro-Mas         HSC      LT-HSC    Lymphoid    MPP/LMPP     Myeloid 
  #        12          25           6          61          13          10          15          14 
fwrite(res_go_mf_lin_dt,fp(out1,"res_hto_signature_go_mf_by_lineage.csv"))


res_go_bp_lin<-lapply(unique(res_lin_dup$lineage), function(lin)enrichGO(bitr(res_lin_dup[padj<0.05&abs(log2FoldChange)>0.6&lineage==lin]$gene,fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 0.05))
names(res_go_bp_lin)<-unique(res_lin_dup$lineage)

saveRDS(res_go_bp_lin,fp(out1,"res_hto_signature_go_bp_by_lineage.rds"))
res_go_bp_lin_dt<-Reduce(rbind,lapply(names(res_go_bp_lin)[!sapply(res_go_bp_lin, is.null)],function(lin)data.table(as.data.frame(res_go_bp_lin[[lin]]))[,lineage:=lin]))
res_go_bp_lin_dt
table(res_go_bp_lin_dt$lineage)
   # B cell          DC Erythro-Mas         HSC    Lymphoid    MPP/LMPP     Myeloid 
   #       57         210          98         443         129          46         155 
fwrite(res_go_bp_lin_dt,fp(out1,"res_hto_signature_go_bp_by_lineage.csv"))


#I) ctrl vs ctrl hto####
source("scripts/utils/new_utils.R")
cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
cbps_c<-subset(cbps,group=="ctrl")


out<-"outputs/08-HTO_signature/pseudobulk_DESeq2_ctrl_hto"
dir.create(out,recursive=T)

#get mtd of interest
mtd<-data.table(cbps_c@meta.data,keep.rownames = "bc")
table(mtd$hto)
# FALSE  TRUE 
# 18521  5824 

mts<-unique(mtd,by=c("sample","hto"))
table(mts$hto)
# FALSE  TRUE 
#     7     8 

#get counts and filter genes lowly express
counts<-as.matrix(cbps_c@assays$RNA@counts)
dim(counts) #34889 24345

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
nrow(counts) # 14836 genes

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample_hto, fun = "sum"))

#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ orig.ident+sex)

dds <- DESeq(dds)

# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
# calculate the vector of coefficient weights in hto
hto <- colMeans(mod_mat[dds$hto == T, ])
hto

# calculate the vector of coefficient weights in non hto
basal <- colMeans(mod_mat[dds$hto == F, ])
basal

# The contrast we are interested in is the difference between hto and not
hto - basal

# get the results for this contrast
res <- results(dds, contrast = hto - basal)
res_dt<-data.table(as.data.frame(res),keep.rownames="gene")
res_dt[padj<0.05&log2FoldChange>0.6] #1537 DEGs
res_dt[,lineage:="all_cbps"]

fwrite(res_dt,fp(out,"res_all_cbps_de_analysis.csv"),sep=";")

#by lineage

Idents(cbps_c)<-"lineage_hmap"
res_lin<-Reduce(rbind,lapply(unique(cbps_c$lineage_hmap[!cbps_c$differentiated]),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_c,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                                 design = ~ orig.ident+sex)
  
  dds <- DESeq(dds)
  
  
  mod_mat <- model.matrix(design(dds), colData(dds))
  hto <- colMeans(mod_mat[dds$hto == T, ])
  basal <- colMeans(mod_mat[dds$hto == F, ])


  res <- results(dds,contrast = hto-basal,alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

fwrite(res_lin,fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"))



#II) same in lga####
cbps_l<-subset(cbps,group=="lga")

out<-"outputs/08-HTO_signature/pseudobulk_DESeq2_lga_hto"
dir.create(out,recursive=T)

#get mtd of interest

mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
table(mtd$hto)
# FALSE  TRUE 
# 16791  6861  
mts<-unique(mtd,by=c("sample","orig.ident"))
table(mts$hto)
# FALSE  TRUE 
#     6     6 
#get counts and filter genes lowly express
counts<-as.matrix(cbps_sub@assays$RNA@counts)
dim(counts) #34889 23652

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
nrow(counts) # 14801 genes

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample_hto, fun = "sum"))

#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ orig.ident+sex)

dds <- DESeq(dds)
resultsNames(dds)
# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
# calculate the vector of coefficient weights in hto
hto <- colMeans(mod_mat[dds$hto == T, ])
hto

# calculate the vector of coefficient weights in non hto
basal <- colMeans(mod_mat[dds$hto == F, ])
basal

# The contrast we are interested in is the difference between hto and not
hto - basal

# get the results for this contrast
res <- results(dds, contrast = hto - basal)

res_dt<-data.table(as.data.frame(res),keep.rownames="gene")
res_dt[padj<0.05&log2FoldChange>0.6] #1128 DEGs

res_dt[,lineage:="all_cbps"]
fwrite(res_dt,fp(out,"res_all_cbps_de_analysis.csv"),sep=";")

#by lineage
Idents(cbps_l)<-"lineage_hmap"
res_lin<-Reduce(rbind,lapply(unique(cbps_c$lineage_hmap[!(cbps_c$differentiated)&cbps_c$lineage_hmap!="LT-HSC"]),function(lin){
  print(lin)
  cbps_sub<-subset(cbps_l,lineage_hmap==lin)
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 
  if(nrow(counts)>0){
    # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                                 design = ~ orig.ident+sex)
  
  dds <- DESeq(dds)
  
  
  mod_mat <- model.matrix(design(dds), colData(dds))
  hto <- colMeans(mod_mat[dds$hto == T, ])
  basal <- colMeans(mod_mat[dds$hto == F, ])


  res <- results(dds,contrast = hto-basal,alpha = 0.05)
  
  return(data.table(as.data.frame(res),keep.rownames="gene")[,lineage:=lin])
  }else{
    return(data.table())
      }
  

  }))

cbps_sub<-subset(cbps_l,lineage_hmap=="LT-HSC")
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
  message(nrow(counts)," genes kept after filtering") 

# Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                                 design = ~ orig.ident)
  
  dds <- DESeq(dds)
  
  
mod_mat <- model.matrix(design(dds), colData(dds))
hto <- colMeans(mod_mat[dds$hto == T, ])
basal <- colMeans(mod_mat[dds$hto == F, ])


res <- results(dds,contrast = hto-basal,alpha = 0.05)
  
res<-data.table(as.data.frame(res),keep.rownames="gene")[,lineage:="LT-HSC"]
res_lin<-rbind(res_lin,res)
fwrite(res_lin,fp(out,"res_pseudobulkDESeq2_by_lineage.csv.gz"))


#IV) ctrl lga hto####
cbps_sub<-subset(cbps,hto==T&group%in%c("ctrl","lga"))

out<-"outputs/08-HTO_signature/pseudobulk_DESeq2_lga_vs_ctrl_hto"
dir.create(out,recursive=T)

#get mtd of interest

mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
table(mtd$hto)
#   TRUE 
#  12685
mts<-unique(mtd,by=c("sample","orig.ident"))
table(mts$hto)
# TRUE 
#  14
#get counts and filter genes lowly express
counts<-as.matrix(cbps_sub@assays$RNA@counts)

counts <- counts[rowSums(counts > 0) >= 100|rowSums(counts > 0)>=ncol(counts)*0.1, ] 
nrow(counts) # 14140 genes

  # Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample_hto, fun = "sum"))

#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample_hto")[colnames(sample_counts),], 
                               design = ~ group+orig.ident+sex)

dds <- DESeq(dds)
resultsNames(dds)

# get the model matrix
mod_mat <- model.matrix(design(dds), colData(dds))
mod_mat
# calculate the vector of coefficient weights in hto
lga <- colMeans(mod_mat[dds$group == "lga", ])
lga

# calculate the vector of coefficient weights in non hto
ctrl <- colMeans(mod_mat[dds$group == "ctrl", ])
ctrl

# The contrast we are interested in is the difference between hto and not
lga - ctrl

# get the results for this contrast
res <- results(dds, contrast = lga - ctrl)

res_dt<-data.table(as.data.frame(res),keep.rownames="gene")
res_dt[padj<0.05] #478
res_dt[,lineage:="all_cbps"]
fwrite(res_dt,fp(out,"res_all_cbps_de_analysis.csv"),sep=";")



#now that we have validated hsc response to stimulation, we can evaluate the differential response of lga HSC compared to control

#see 09-


