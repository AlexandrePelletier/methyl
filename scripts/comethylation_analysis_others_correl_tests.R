#use Hmisc because calculate correlation value and significance in one 
renv::install("Hmisc")
library(Hmisc)
#need a matrix obs/var
methg0_mat<-as.matrix(dcast(methg0,sample~gene,value.var = "gene_meth_scaled"),rownames = "sample")
dim(methg0_mat) #108 7401
7401^2 #54M de tests de correl
head(methg0_mat[,1:10])
res<-rcorr(methg0_mat,type="spearman")
res$P.adj<-p.adjust(res$P,method = "BH")
sum(res$P.adj<0.05,na.rm = T) #3,7M of signif correl

saveRDS(res,fp(out,"gene_correl_spearman.rds"))

108^2 #12k tests de correl
res_sample<-rcorr(t(methg0_mat),type="spearman")
res_sample$P.adj<-p.adjust(res_sample$P,method = "BH")
sum(res_sample$P.adj<0.05,na.rm = T) #12k of signif correl
saveRDS(res_sample,fp(out,"sample_correl_spearman.rds"))

      #c) clusteriser la matrice de correlation pour identifier des modules
res<-readRDS(fp(out,"gene_correl_spearman.rds"))
res_sample<-readRDS(fp(out,"sample_correl_spearman.rds"))

#1) test with pheatmap : is there module of genes ?
renv::install("pheatmap")
mtd<-fread("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")
library(pheatmap)
pheatmap(
  t(methg0_mat), 
  annotation_col = data.frame(mtd,row.names = "sample")[colnames(methg0_mat),c("Group_Sex","Group_Complexity_Fac")],
  clustering_distance_cols = as.dist(1 - res_sample$r),
  clustering_distance_rows = as.dist(1 - res$r),
  show_rownames = F,show_colnames = F
  ) #doesnt work
