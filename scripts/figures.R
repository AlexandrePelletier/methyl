#figures
library(clusterProfiler)
library(org.Hs.eg.db)


#pathways enrichment in CTRL vs LGA
source("scripts/utils/new_utils.R")
res<-fread("analyses/model14_without_iugr/2020-12-15_all_res_genescore1.csv")

res_cl<-res[compa=="C.L"]

geneList<-unique(res_cl,by="gene")$gene_score
names(geneList)<-unique(res_cl,by="gene")$gene
geneList<-sort(geneList,decreasing = T)
head(geneList,20)
genes.df<-bitr(names(geneList),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
genes.df$GeneScore<-geneList[genes.df$SYMBOL]
geneList.Entrez<-genes.df$GeneScore
names(geneList.Entrez)<-genes.df$ENTREZID

res_kegg<-gseKEGG(geneList= rank(geneList.Entrez),
            organism     = 'hsa', 
            minGSSize    = 50,
            pvalueCutoff = 0.001,
            verbose = FALSE)

res_kegg_dt<-data.table(res_kegg@result)

dotplot(res_kegg,showCategory=20)

emapplot(res_kegg,showCategory =30 )

res_kegg_0.01<-gseKEGG(geneList= rank(geneList.Entrez),
            organism     = 'hsa', 
            minGSSize    = 50,
            pvalueCutoff = 0.01,
            verbose = FALSE)

res_kegg_dt_0.01<-data.table(res_kegg_0.01@result)
emapplot(res_kegg_0.01,showCategory =68 )

#CXCL12/CXL14
source("scripts/utils/methyl_utils.R")
unique(res[gene=="CXCL12"][order(pval)],by=c("gene","compa"))
plotMeth(res[gene=="CXCL12"&pval<0.001]$cpg_id)


unique(res[gene=="CXCL14"][order(pval)],by=c("gene","compa"))
plotMeth(res[gene=="CXCL14"&pval<0.001]$cpg_id)

res_genes<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-09-01_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
res_genes[gene=="CXCL12"]


