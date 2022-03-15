#Gene score calculation and validation
source("scripts/utils/new_utils.R")
source("scripts/utils/methyl_utils.R")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

out<-"outputs/03-pathway_analysis"
dir.create(out)


resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")

#OR KEGG
plot(density(resg$gene_score_add))
abline(v=300)
resg[gene_score_add>300]
res_or_kegg<-enrichKEGG(bitr(resg[gene_score_add>300]$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,
                 organism = "hsa",pvalueCutoff = 0.05)
as.data.frame(res_or_kegg)


emapplot(pairwise_termsim(res_or_kegg))

#GSE KEGG
gene_scores<-resg$gene_score_add
names(gene_scores)<-resg$gene
gene_scores<-sort(gene_scores,decreasing = T)
head(gene_scores)

genes.df<-bitr(names(gene_scores),
               fromType = 'SYMBOL',
               toType = 'ENTREZID',
               OrgDb = org.Hs.eg.db)
head(genes.df)
head(gene_scores)
gene_scores<-gene_scores[genes.df$SYMBOL]
names(gene_scores)<-genes.df$ENTREZID
res_gsea_kegg<- gseKEGG(geneList     = rank(gene_scores),
                        organism     = 'hsa', 
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        verbose = FALSE)

dotplot(res_gsea_kegg,showCategory=63)
gsea_kegg_dt<-data.table(as.data.frame(res_gsea_kegg))
gsea_kegg_dt[p.adjust<0.001]
gsea_kegg_dt[p.adjust<0.001,gene_score.avg:=mean(resg$gene_score_add[resg$gene %in% tr(core_enrichment,tradEntrezInSymbol = T)],na.rm=T),.(ID)]

dotplot(res_gsea_kegg,x=gsea_kegg$gene_score.avg,showCategory=63)

emapplot(pairwise_termsim(res_gsea_kegg,showCategory = 63))
saveRDS(res_gsea_kegg,fp(out,"res_gsea_kegg.rds"))
fwrite(gsea_kegg_dt,fp(out,"res_gsea_kegg.csv"))

#GSE GO

res_gsea_go<- gseGO(geneList     = rank(gene_scores), 
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
nrow(as.data.frame(res_gsea_go))#5986

dotplot(res_gsea_go,showCategory=20)
gsea_go_dt<-data.table(as.data.frame(res_gsea_go))
gsea_go_dt[Description=="regulation of growth"]

saveRDS(res_gsea_go,fp(out,"res_gsea_go.rds"))
fwrite(gsea_go_dt[order(p.adjust)],fp(out,"res_gsea_go.csv"))

dotplot(res_gsea_go,x=gsea_go[order(p.adjust)]$gene_score.avg[1:40],showCategory=40)

emapplot(pairwise_termsim(res_gsea_go,showCategory = 40),showCategory = 40)
res_gsea_go<- gseGO(geneList     = rank(gene_scores), 
                    ont="BP",
                        minGSSize    = 10,maxGSSize = 600,
                        pvalueCutoff = 1,
                        eps = 0,
                        OrgDb = org.Hs.eg.db)
nrow(as.data.frame(res_gsea_go))#5986

dotplot(res_gsea_go,showCategory=20)
gsea_go_dt<-data.table(as.data.frame(res_gsea_go))
gsea_go_dt[Description=="regulation of growth"]

saveRDS(res_gsea_go,fp(out,"res_gsea_go_bp_all.rds"))
fwrite(gsea_go_dt[order(p.adjust)],fp(out,"res_gsea_go_bp_all.csv"))

#GSE GWAS

source("scripts/utils/new_utils.R")
library(clusterProfiler)
library(enrichplot)
out<-"outputs/03-pathway_analysis"

resg<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
gwas_genes<-fread("../methyl/ref/gwas/list_ref_GWAS_042721.csv",
                  select = c(2,3,4),
                  col.names = c("reported_gene","mapped_genes","disease_trait"),
                  skip = 1)

#☺1) with reported genes
gwas_genes[reported_gene=="NR",reported_gene:=NA]
gwas_genes_rep<-unique(gwas_genes[!is.na(reported_gene)],by=c("reported_gene","disease_trait"))
gwas_genes_rep #92k gene-trait associations
length(unique(gwas_genes_rep$disease_trait)) #4k5 disease traits
gwas_genes_rep[,n.gene.trait:=.N,by="disease_trait"]
summary(unique(gwas_genes_rep[,.(disease_trait,n.gene.trait)])$n.gene.trait)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    2.00    5.00   20.38   11.00 1532.00 

# remove traits if < 10 genes 
gwas_genes_rep10<-gwas_genes_rep[n.gene.trait>=10]
gwas_genes_rep10 #81k gene-trait associations
fwrite(gwas_genes_rep10,"../methyl/ref/gwas/reported_gene_traits_GWAS_10genes.csv")
gwas_genes_rep10<-fread("../methyl/ref/gwas/reported_gene_traits_GWAS_10genes.csv")

length(unique(gwas_genes_rep10$disease_trait)) #1299 disease traits
gwas_genes_rep10[,n.gene.trait:=.N,by="disease_trait"]
summary(unique(gwas_genes_rep10[,.(disease_trait,n.gene.trait)])$n.gene.trait)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.00   14.00   23.00   62.75   50.00 1532.00 

gene_scores<-resg$gene_score_add
names(gene_scores)<-resg$gene
gene_scores<-sort(gene_scores,decreasing = T)

res_gsea_gwas<-GSEA(geneList = rank(gene_scores),
                    TERM2GENE = gwas_genes_rep10[,.(disease_trait,reported_gene)],
                    maxGSSize    = 500,
                    eps = 0,
                    pvalueCutoff = 0.01
                    )
nrow(as.data.frame(res_gsea_gwas))#86/1300
dotplot(res_gsea_gwas,showCategory=20)
emapplot(pairwise_termsim(res_gsea_gwas,showCategory = 86),showCategory = 86)

gsea_gwas<-data.table(as.data.frame(res_gsea_gwas))
gsea_gwas[,gene_score.avg:=mean(resg$gene_score_add[resg$gene %in% tr(core_enrichment,tradEntrezInSymbol = F)],na.rm=T),.(ID)]

saveRDS(res_gsea_gwas,fp(out,"res_gsea_gwas.rds"))
fwrite(gsea_gwas[order(p.adjust)],fp(out,"res_gsea_gwas.csv"))


#2) if not satisfying, with mapped gene #satisfying

#PERMUT Patway====
res_perm<-fread("outputs/02-gene_score_calculation_and_validation/res_1000perm_genescore_add.csv.gz")
 #pathway
library(clusterProfiler)
library(org.Hs.eg.db)

#kegg
res_k_perm<-Reduce(rbind,mclapply(1:1000,function(i){
  resg<-res_perm[perm==i]
  gene_scores<-resg$gene_score_add
  names(gene_scores)<-resg$gene
  gene_scores<-sort(gene_scores,decreasing = T)
  
  genes.df<-bitr(names(gene_scores),
                 fromType = 'SYMBOL',
                 toType = 'ENTREZID',
                 OrgDb = org.Hs.eg.db)
  gene_scores<-gene_scores[genes.df$SYMBOL]
  names(gene_scores)<-genes.df$ENTREZID
  res_gsea_kegg<- data.table(as.data.frame(gseKEGG(geneList     = rank(gene_scores),
                        organism     = 'hsa', 
                        minGSSize    = 50,
                        pvalueCutoff = 1,
                        verbose = FALSE)))
  
  return(res_gsea_kegg[,perm:=i][,.(ID,p.adjust,perm)])}))

fwrite(res_k_perm,"outputs/03-pathway_analysis/res_1000perm_kegg.csv.gz")

#check res
source("scripts/utils/new_utils.R")
res_k_perm<-fread("outputs/03-pathway_analysis/res_1000perm_kegg.csv.gz")
res_k_perm
res_k<-fread("outputs/03-pathway_analysis/res_gsea_kegg.csv")
res_kp<-merge(res_k,res_k_perm[ID%in%res_k$ID],all=T)
res_kp[,p.perm:=sum(p.adjust[is.na(perm)]>=p.adjust[!is.na(perm)])/1000,by="ID"]

res_kpf<-res_kp[is.na(perm)][,-"perm"]
res_kpf[p.perm<0.05]# 32/63 are signif

fwrite(res_kpf,"outputs/03-pathway_analysis/res_gsea_kegg_perm.csv")
#go bp
system("nohup Rscript scripts/03A-go_bp_perm.r > logs/03A-go_bp_perm.log &")
#run go_bp_perm.r


#gwas
system("nohup Rscript scripts/03A-gwas_perm.r > logs/03A.log &")
#run gwas_perm.r
res_gw_perm<-fread("outputs/03-pathway_analysis/res_gsea_gwas_perm.csv")

res_gw_perm[p.perm<0.05] #36/87
