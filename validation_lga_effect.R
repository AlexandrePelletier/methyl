library("Seurat")
source("scripts/utils/new_utils.R")

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")

degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")

degs[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5]$gene
genes<-c('SESN2','SOCS3','EGR1','ID1',"KLF2","IER2","IER3","DUSP2","PLK3","JUNB","ZFP36")
gene_exp<-merge(data.table(cbps@meta.data,keep.rownames = "bc"),
           melt(data.table(data.frame(t(as.matrix(cbps@assays$SCT@data[genes,]))),keep.rownames = "bc"),
                variable.name = "gene",value.name = "expression",id.vars = "bc"))

ggplot(gene_exp[lineage_hmap=="HSC"])+geom_boxplot(aes(x=hto,y=expression,fill=group))+facet_wrap("gene")
ggplot(gene_exp[lineage_hmap=="HSC"])+geom_boxplot(aes(x=hto,y=expression,fill=group,group=sample_hto))+facet_wrap("gene",scales="free")

ggplot(gene_exp[lineage_hmap=="HSC"&hto==T])+geom_boxplot(aes(x=group,y=expression,fill=sex,group=sample_hto))+facet_wrap("gene",scales="free")

gene_exp[,avg_expr:=median(expression),by=.(sample_hto,lineage_hmap,gene)]
gene_exp_s<-unique(gene_exp,by = c("sample_hto","lineage_hmap","gene"))
ggplot(gene_exp_s[lineage_hmap=="HSC"])+geom_boxplot(aes(x=hto,y=expression,fill=group))+facet_wrap("gene",scales = "free")


degs_sc<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_scEdgeR_by_lineage.csv.gz")
degs_sc[p_val_adj<0.05&abs(avg_logFC)>0.5&lineage_hmap=='HSC']$gene
