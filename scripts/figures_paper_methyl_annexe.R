##annexe paper methyl
source("../methyl/scripts/utils/new_utils.R")
library(limma)
library(pheatmap)
out<-"../methyl/outputs/figures_paper_meth"
dir.create(out)


meth<-fread("../methyl/datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
mtd<-fread("../methyl/datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")

#genescore heatmap
meth_df<-melt(meth,id.vars = 1,
              variable.name = "sample",
              value.name = "methylation")
cpgs_weight<-fread("../methyl/ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_weight[,cpg_weight:=RegWeight+LinksWeight]
meth_df<-merge(meth_df,cpgs_weight[,.(locisID,gene,cpg_weight)],by="locisID")
meth_df[,gene_meth:=sum(methylation*cpg_weight)/sum(cpg_weight),by=c("sample","gene")]
fwrite(meth_df,fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")
meth_df<-fread(fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")

meth_genes<-unique(meth_df,by=c("sample","gene"))
plot(density(meth_genes$gene_meth))
meth_genes[,gene_meth_scaled:=scale(gene_meth),by=c("sample")]

meth_genes_mat<-dcast(meth_genes[,-c("locisID","hi.conf.genes","methylation","gene_meth")],sample~gene)
meth_genes_mat[1:10,1:10]

meth_genes_mat<-as.matrix(data.frame(meth_genes_mat,row.names = "sample"))
meth_genes_mat<-t(meth_genes_mat)
meth_genes_mat[1:10,1:10]
dim(meth_genes_mat)


pheatmap(meth_genes_mat,
         annotation_col = data.frame(mtd,row.names = "sample")[colnames(meth_genes_mat),c("Group_Sex","Group_Complexity_Fac")],
         show_rownames = F,show_colnames = F
        ) #doesn work

#top1000 sd Ctrl LGA samples
meth_genes<-merge(meth_genes,mtd,by="sample")
meth_genes_cl<-meth_genes[Group_name%in%c("C","L")]
meth_genes_cl[,sd:=sd(gene_meth_scaled),by="gene"]
meth_genes_cl[,top1000:=sd>=sort(sd,decreasing = T)[1000],by="sample"]
pheatmap(meth_genes_mat[rownames(meth_genes_mat)%in%unique(meth_genes_cl[top1000==T]$gene),unique(meth_genes_cl$sample)],
         annotation_col = data.frame(mtd,row.names = "sample")[unique(meth_genes_cl$sample),c("Group_name","Gender","Library_Complexity")],
         show_rownames = F,show_colnames = F
        )

#hi sd in complex samples
meth_genes_cl[Group_Complexity_Fac==4,sd.complex:=sd(gene_meth_scaled),by="gene"]
meth_genes_cl[,top1000.complex:=sd.complex>=sort(sd.complex,decreasing = T)[1000],by="sample"]

pheatmap(meth_genes_mat[rownames(meth_genes_mat)%in%unique(meth_genes_cl[top1000.complex==T]$gene),unique(meth_genes_cl$sample)],
         annotation_col = data.frame(mtd,row.names = "sample")[unique(meth_genes_cl$sample),c("Group_name","Gender","Library_Complexity")],
         show_rownames = F,show_colnames = F
        )


#hi conf genes
meth_df[,hi.conf.genes:=sum(cpg_weight>=2)>5,by=c("sample","gene")]

meth_genes_hi<-unique(meth_df,by=c("sample","gene"))[hi.conf.genes==T]
meth_genes_hi[,gene_meth_scaled:=scale(gene_meth),by=c("sample")]

meth_genes_mat_hi<-dcast(meth_genes_hi,sample~gene,value.var="gene_meth_scaled")
meth_genes_mat_hi[1:10,1:10]

meth_genes_mat_hi<-as.matrix(data.frame(meth_genes_mat_hi,row.names = "sample"))
meth_genes_mat_hi<-t(meth_genes_mat_hi)
meth_genes_mat_hi[1:10,1:10]
dim(meth_genes_mat_hi)


pheatmap(meth_genes_mat_hi[rownames(meth_genes_mat_hi)%in%unique(meth_genes_cl[top1000.complex==T]$gene),unique(meth_genes_cl$sample)],
         annotation_col = data.frame(mtd,row.names = "sample")[unique(meth_genes_cl$sample),c("Group_name","Gender","Library_Complexity")],
         show_rownames = F,show_colnames = F
        )

#see key genes

ggplot(meth_genes[gene%in%c("SOCS3","HES1","JUN")&Group_name%in%c("C","L")])+
geom_boxplot(aes(x=Group_Sex,y=gene_meth))+facet_wrap("gene")



#validation genescore
res_de<-fread("../singlecell/outputs/07-DEGs_LGA_stress/pseudobulk_deseq2_all_cbps/res_de_analysis_all_genes.csv")

res_cl_merge<-merge(res_de,unique(res_cl[order(gene,pval)],by="gene"))

ggplot(res_cl_merge)+geom_boxplot(aes(x=padj<0.05,y=GeneScore))          

wilcox.test(res_cl_merge[padj<=0.05]$GeneScore,res_cl_merge[padj>0.05]$GeneScore)

res_de_sc<-rbind(fread("../singlecell/outputs/08-DEGs_LGA_no_stress/sc_edger_deseq2_all_cbps/res_de_analysis_all_genes.csv")[,condition:="basale"],
              fread("../singlecell/outputs/07-DEGs_LGA_stress/sc_edger_deseq2_all_cbps/res_de_analysis_all_genes.csv")[,condition:="stress"])

res_cl_sc_merge<-merge(res_de_sc,unique(res_cl,by="gene"),allow.cartesian=TRUE)

ggplot(res_cl_sc_merge)+geom_boxplot(aes(x=p_val_adj<0.01&abs(avg_logFC)>0.6,y=GeneScore))+facet_wrap("condition")          

#doexnt work.. pb of degs ?
genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")
res_de[gene%in%genes_of_interest]
res_de_sc[gene%in%genes_of_interest]

#not really, try correl in hsc : 

res_de_hsc<-rbind(fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_by_lin/res_de_analysis_all_genes.csv")[lineage=="HSC"][,condition:="basale"],
              fread("../singlecell/outputs/07-DEGs_LGA_stress/pseudobulk_deseq2_by_lin/res_de_analysis_all_genes.csv")[lineage=="HSC"][,condition:="stress"])

res_cl_hsc_merge<-merge(unique(res_de_hsc,by=c("gene","condition")),unique(res_cl,by="gene"),allow.cartesian=TRUE)

ggplot(res_cl_hsc_merge)+geom_boxplot(aes(x=padj<0.2,y=GeneScore))+facet_wrap("condition")          

res_de_hsc[padj<0.05]
res_de_hsc[gene%in%genes_of_interest]


#doesnt validate genescore, try with res_cflf
res_clf<-fread("../methyl/outputs/model14_without_iugr/2020-09-16_res_CF.LF_with_GeneScore_and_permut.csv")

#stress
res_de_clf_hto<-fread("../singlecell/outputs/07-DEGs_LGA_stress/pseudobulk_deseq2_all_progens_female/res_de_analysis_all_genes.csv")
res_de_clf_hto[padj<0.2]
res_clf_hto_merge<-merge(res_de_clf_hto,unique(res_clf[order(gene,pval)],by="gene"))

p1<-ggplot(res_clf_hto_merge)+geom_boxplot(aes(x=padj<0.2,y=GeneScore),outlier.shape = NA) +scale_y_continuous(limits = c(-50,75))
p2<-ggplot(res_clf_hto_merge)+geom_boxplot(aes(x=padj<0.2,y=meth.change),outlier.shape = NA)+scale_y_continuous(limits = c(-20,60))
p3<-ggplot(res_clf_hto_merge)+geom_boxplot(aes(x=padj<0.2,y=-log10(pval)),outlier.shape = NA)+scale_y_continuous(limits = c(0,4))     
p1+p2+p3

#no stress
res_de_clf<-fread("../singlecell/outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_progens_female/res_de_analysis_all_genes.csv")
res_de_clf[padj<0.2]
#stress
res_clf_merge<-merge(res_de_clf,unique(res_clf[order(gene,pval)],by="gene"))

p1<-ggplot(res_clf_merge)+geom_boxplot(aes(x=padj<0.2,y=GeneScore),outlier.shape = NA) +scale_y_continuous(limits = c(-50,75))
p2<-ggplot(res_clf_merge)+geom_boxplot(aes(x=padj<0.2,y=meth.change),outlier.shape = NA)+scale_y_continuous(limits = c(-20,60))
p3<-ggplot(res_clf_merge)+geom_boxplot(aes(x=padj<0.2,y=-log10(pval)),outlier.shape = NA)+scale_y_continuous(limits = c(0,4))     
p1+p2+p3
ggplot(res_clf_merge)+geom_boxplot(aes(x=pvalue<0.01,y=GeneScore)) 
#ccl : un peu mieux mais pas ouf quand meme

#validation with degs ctrl vs lga basic (with hto in covariate) #see 04-script, but spoiler : validate pas



#Over repre test GO BP and KEGG pathway
#basale
#kegg
pathways<-read.gmt("../methyl/ref/msigdb/c2.cp.kegg.v7.1.symbols.gmt")
pathways_list<-split(pathways$gene,f = pathways$term)

res_kegg<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Upreg"]$gene ,
                     terms_list = pathways_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_kegg[padj<0.5]#0

#with hypermet-downreg
res_kegg2<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Downreg"]$gene ,
                     terms_list = pathways_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_kegg2[padj<0.5]#0


#GO
go_term<-read.gmt("../methyl/ref/msigdb/c5.go.bp.v7.3.symbols.gmt")
go_list<-split(go_term$gene,f = go_term$term)

res_go<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Upreg"]$gene ,
                     terms_list = go_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_go[padj<0.5]#0
#with hypermet-downreg
res_go<-over_repr_test_multi(genes_of_interest =res_cl_merge[correl=="Hypermet-Downreg"]$gene ,
                     terms_list = go_list ,
                     size_universe = nrow(res_cl_merge)
                    )

res_go[padj<0.5]#0

#in stress
res_clh_merge[,gs_scaled:=scale(GeneScore,center = F)]
ggplot(res_clh_merge)+geom_density(aes(x=gs_scaled))
ggplot(res_clh_merge)+geom_density(aes(x=avg_logFC))

res_clh_merge[,correl:=ifelse(gs_scaled>1&avg_logFC>1,"Hypermet-Upreg",
                             ifelse(gs_scaled>1&avg_logFC<(-1),"Hypermet-Downreg",
                              ifelse(gs_scaled<(-1)&avg_logFC<(-1),"Hypomet-Downreg",
                               ifelse(gs_scaled<(-1)&avg_logFC>1,"Hypomet-Upreg","no change"))))]
ggplot(res_clh_merge) + geom_point(aes(x=avg_logFC,y=gs_scaled,col=correl))

#kegg
pathways<-read.gmt("../methyl/ref/msigdb/c2.cp.kegg.v7.1.symbols.gmt")
pathways_list<-split(pathways$gene,f = pathways$term)

res_kegg<-over_repr_test_multi(genes_of_interest =res_clh_merge[correl=="Hypermet-Upreg"]$gene ,
                     terms_list = pathways_list ,
                     size_universe = nrow(res_clh_merge)
                    )

res_kegg[padj<0.5]#0

#with hypermet-downreg
res_kegg2<-over_repr_test_multi(genes_of_interest =res_clh_merge[correl=="Hypermet-Downreg"]$gene ,
                     terms_list = pathways_list ,
                     size_universe = nrow(res_clh_merge)
                    )

res_kegg2[padj<0.5]#0


#GO
go_term<-read.gmt("../methyl/ref/msigdb/c5.go.bp.v7.3.symbols.gmt")
go_list<-split(go_term$gene,f = go_term$term)

res_go<-over_repr_test_multi(genes_of_interest =res_clh_merge[correl=="Hypermet-Upreg"]$gene ,
                     terms_list = go_list ,
                     size_universe = nrow(res_clh_merge)
                    )

res_go[padj<0.5]#0
#with hypermet-downreg
res_go<-over_repr_test_multi(genes_of_interest =res_clh_merge[correl=="Hypermet-Downreg"]$gene ,
                     terms_list = go_list ,
                     size_universe = nrow(res_clh_merge)
                    )

res_go[padj<0.5]#0


#GSEA [to do]
#score correl by gene = distance fold change obs with fold change calc ~ linear regression
summary(lm(res_cl_merge$avg_logFC~res_cl_merge$GeneScore))
#a = -0.0005153
#b=  -0.0162490

PredFoldChange<-function(gene_score,a=-0.0005153,b=-0.0162490)a*gene_score+b
PredFoldChange(300)
res_cl_merge[,avg_logFC_pred:=PredFoldChange(GeneScore)]
ggplot(res_cl_merge)+geom_point(aes(x=GeneScore,y=avg_logFC_pred))
res_cl_merge[,avg_logFC_dist:=abs(avg_logFC_pred-avg_logFC)]
res_cl_merge[order(avg_logFC_dist)]
#doesnt worj because correl based on data






#degs cf lf
out<-"outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_all_progens_female"
dir.create(out,recursive=T)

cbps_cl<-subset(cbps,hto==F&group%in%c("ctrl","lga")&ambigous==F&orig.ident!="cd34_hto1_0C1I1L")

#all progens
lineages<-c("HSC","MPP","Erythroid","Lymphoid","Myeloid")
cbps_sub<-subset(cbps_cl,lineage%in%lineages&sex=="F")
#get mtd of interest
mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
mts<-unique(mtd,by=c("sample","orig.ident"))
table(mts$group)
#get counts and filter genes lowly express
counts<-as.matrix(cbps_sub@assays$RNA@counts)
dim(counts) 

counts <- counts[rowSums(counts > 0) >= 100, ] 
message(nrow(counts)," genes kept after filtering") 

# Aggregate across cluster-sample groups
sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                     groupings = mtd$sample, fun = "sum"))
#DEseq2_analysis
dds <- DESeqDataSetFromMatrix(sample_counts, 
                               colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                               design = ~ group)

dds <- DESeq(dds)

res <- results(dds,
             contrast = c("group","lga","ctrl"),
              alpha = 0.05)
res<-data.table(as.data.frame(res),keep.rownames="gene")[is.na(padj),padj:=1][,lineage:='all_progens']
fwrite(res,fp(out,"res_de_analysis_all_genes.csv"),sep=";")

#by lin
out<-"outputs/08-DEGs_LGA_no_stress/pseudobulk_deseq2_female_by_lineage"
dir.create(out,recursive=T)

res_lin<-Reduce(rbind,lapply(lineages,function(lin){
  print(lin)
  cbps_sub<-subset(cbps_cl,lineage==lin&sex=="F")
  #get mtd of interest
  mtd<-data.table(cbps_sub@meta.data,keep.rownames = "bc")
  mts<-unique(mtd,by=c("sample","orig.ident"))
  #get counts and filter genes lowly express
  counts<-as.matrix(cbps_sub@assays$RNA@counts)
  dim(counts) 

  counts <- counts[rowSums(counts > 0) >= 100, ] 
  message(nrow(counts)," genes kept after filtering") 
  
  # Aggregate across cluster-sample groups
  sample_counts <- t(aggregate.Matrix(t(counts[,mtd$bc]), 
                       groupings = mtd$sample, fun = "sum"))
  #DEseq2_analysis
  dds <- DESeqDataSetFromMatrix(sample_counts, 
                                 colData = data.frame(mts,row.names="sample")[colnames(sample_counts),], 
                                 design = ~ group)
  
  dds <- DESeq(dds)
  
  res <- results(dds,
               contrast = c("group","lga","ctrl"),
                alpha = 0.05)
  return(data.table(as.data.frame(res),keep.rownames="gene")[is.na(padj),padj:=1][,lineage:=lin])

  }))


fwrite(res_lin,fp(out,"res_de_analysis_all_genes.csv"),sep=";")
res_lin[padj<0.2]
res_lin[gene%in%genes_of_interest&padj<0.2]
res_lin[gene%in%genes_of_interest&padj<0.5]

