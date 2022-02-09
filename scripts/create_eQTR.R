library(data.table)
library(stringr)
library(biomartr)
eqtl<-fread("ref/GTEx_Analysis_v8.metasoft.txt.gz")
#I) recover col gene, chr, pos 
#1) in hg38 :
eqtl[,chr:=str_extract(RSID,"chr[0-9XY]+")]
eqtl[,pos:=str_extract(RSID,"_[0-9]+")]
eqtl[,pos:=as.numeric(str_extract(pos,"[0-9]+"))]
eqtl[,gene_id:=str_extract(RSID,"ENSG[0-9]+")]
eqtl

#trans in symbol :  with hsapiens_gene_ensembl when biomartr server will work
genes_translator <- biomart( genes      = unique(eqtl$gene_id), # genes that we wanted info
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "hsapiens_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = "hgnc_symbol", # attributes were selected with biomartr::getAttributes()
                                filters    = "ensembl_gene_id")

genes_translator<-data.table(genes_translator)
genes_translator<-genes_translator[,gene_id:=ensembl_gene_id][,-"ensembl_gene_id"]
genes_translator<-genes_translator[,gene:=hgnc_symbol][,-"hgnc_symbol"]
eqtl<-merge(eqtl,genes_translator,all.x=T,by="gene_id")


#2) trans in hg38 :
translator<-fread("ref/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
                  select = c(1,8))
translator

eqtl[,variant_id:=str_remove(RSID,",ENSG[0-9]+\\.[0-9]+")]


eqtl<-merge(eqtl,translator,all.x=T,by="variant_id")


eqtl[,chr.hg19:=paste0("chr",str_extract(variant_id_b37,"^[0-9XY]{1,2}"))]

eqtl[,pos.hg19:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2))]
eqtl

fwrite(eqtl,"ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz")


#II) Create eQTR
#+/- 500pb
eqtl<-fread("ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz")
eqtl[,start.eQTR:=pos-500]
eqtl[,end.eQTR:=pos+500]

eqtl[,start.eQTR.hg19:=pos-500]
eqtl[,end.eQTR.hg19:=pos+500]

#anno for tissue_wide eQTL : 
# - eQTL with low Heterogeneity : Qvalue <0.1 ;
# - eQTL in >30/49 tissue =>  mval>0.8 & pval <10^-4 in >30 tissue
# - PVALUE_FE < 10^-30
# -at least 20 tissues with mval > 0.8
eqtl[,n.tissue_sig:=rowSums(eqtl[,.SD,.SDcols=str_detect(colnames(eqtl),"mval")]>0.8 & eqtl[,.SD,.SDcols=str_detect(colnames(eqtl),"pval")]>0.0001,na.rm = T)]
eqtl[,tissue_wide:=PVALUE_Q>0.1&n.tissue_sig>30&PVALUE_FE<10e-30]
nrow(eqtl[tissue_wide==T]) #552921

#anno for liver
eqtl[,liver:=mval_Liver>0.8&pval_Liver>10^-4]
nrow(eqtl[liver==T]) #3578315
fwrite(eqtl,"ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz")

#save eQTR in bed format
eqtl<-fread("ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz")

fwrite(eqtl[!is.na(gene)&gene!=""&liver==T][,chr:=chr.hg19][,start:=start.eQTR.hg19][,end:=end.eQTR.hg19][,.(chr,start,end,gene)][order(chr,start)],"ref/eQTL/liver_eQTR_hgnc_gene_hg19.bed",sep="\t")
fwrite(eqtl[!is.na(gene)&gene!="" &tissue_wide==T][,chr:=chr.hg19][,start:=start.eQTR.hg19][,end:=end.eQTR.hg19][,.(chr,start,end,gene)][order(chr,start)],"ref/eQTL/tissue_wide_eQTR_hgnc_gene_hg19.bed",sep="\t")


#anno for whole blood and save eQTR in bed format
eqtl<-fread("ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz")
eqtl[,blood:=mval_Whole_Blood>0.8&pval_Whole_Blood>10^-4]

nrow(eqtl[blood==T]) #2151314
fwrite(eqtl,"ref/eQTL/GTEx_Analysis_v8.metasoft_annotated.csv.gz")

fwrite(eqtl[!is.na(gene)&gene!=""&(blood==T|tissue_wide==T)][,chr:=chr.hg19][,start:=start.eQTR.hg19][,end:=end.eQTR.hg19][,.(chr,start,end,gene)][order(chr,start)],
       "ref/eQTL/whole_blood_and_tissue_wide_eQTR_hgnc_gene_hg19.bed",sep="\t")

