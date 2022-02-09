#DMR 
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(reticulate)
library("org.Hs.eg.db")
source("scripts/utils/new_utils.R")

out<-"outputs/01B-DMR_comb-p"
dir.create(out)

res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz",sep="\t")
meth<-fread("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",
            select = c("id","chr","start"),
            col.names = c("cpg_id","chr","pos"))

res<-merge(res,meth[,.(cpg_id,chr,pos)],all.x = T)
fwrite(res[,start:=pos][,end:=pos+1][order(chr,start)][!is.na(start)][,.(chr,start,end,P.Value)],fp(out,"res_cpgs.bed"),sep="\t")


reticulate::py_install("toolshed",pip=T,envname = "/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/")
reticulate::use_python("/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin/python")


old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path,"/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin", sep = ":"))
Sys.getenv("PATH")

cmd<-paste("python tools/comb-p/combined-pvalues/cpv/pipeline.py -p outputs/01B-DMR_comb-p/ctrl_lga --region-filter-p 0.1 --anno hg19 --seed 0.01 --dist 1000",
           fp(out,"res_cpgs.bed"))
system(cmd)
#
res_dmr<-fread("outputs/02-gene_score_calculation_and_validation/ctrl_lga.regions-p.bed.gz")
res_dmr<-res_dmr[,chr:=`#chrom`][,-1][,.(chr,start,end,min_p,n_probes,z_p,z_sidak_p)]
#merge with gene

dmr<-makeGRangesFromDataFrame(res_dmr[,1:3])
anno_dmr <- annotatePeak(dmr, tssRegion=c(-2000, 2000),
                       TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                       annoDb="org.Hs.eg.db")
anno_dmr<-as.data.frame(anno_dmr)
head(anno_dmr,10)
fwrite(anno_dmr,fp(out,"chipseaker_anno_dmr.csv.gz"),sep=";")
anno_dmr<-fread(fp(out,"chipseaker_anno_dmr.csv.gz"),sep=";")

dmr_tss<-data.table(anno_dmr)[,chr:=seqnames][,tss_dist:=distanceToTSS][,gene:=SYMBOL][,.(chr,start,end,tss_dist,gene)]

res_dmr<-merge(res_dmr,dmr_tss,all.x=T)
res_dmr[z_sidak_p<0.05&abs(tss_dist)<2000]$gene
res_dmr[gene=="SOCS3"]
fwrite(res_dmr,fp(out,"res_dmr_anno.csv.gz"),sep=";")
