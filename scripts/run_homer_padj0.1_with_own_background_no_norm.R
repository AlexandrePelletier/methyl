source("scripts/utils/new_utils.R")
set.seed(1234)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/apelletier/HSPC_EpiStress/Data/tools/homer/bin/", sep = ":"))

out<-"outputs/03B-motif_analysis_padj0.1_no_norm"
dir.create(out)
bed_path<-fp(out,"40bp_win_dmc_padj0.01_meth.change25")
background_path<-fp(out,"50k_background_cpgs.bed")


res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
res[adj.P.Val<0.1&abs(logFC)>25] #1255
cpgs<-fread("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",
            select = c("id","chr","start"),
            col.names = c("cpg_id","chr","pos"))
res<-merge(res,cpgs)
res[,start:=pos-20][start<1,start:=1][,end:=pos+20]
res<-res[!is.na(chr)&!is.na(pos)]


random_cpgs<-res[cpg_id%in%sample(cpg_id[!(adj.P.Val<0.1&abs(logFC)>25)],50000),.(chr,start,end)][order(chr,start)]


fwrite(random_cpgs,
       background_path,sep="\t",
       col.names=F)

n_cpgs<-nrow(res[adj.P.Val<0.1&abs(logFC)>25]) 

message("analyzing motif enrichment in ",n_cpgs,"DMCs versus 50k random cpgs background")


fwrite(res[adj.P.Val<0.1&abs(logFC)>25][,.(chr,start,end)][order(chr,start)],
     bed_path,sep="\t",
     col.names=F)

cmd<-paste("findMotifsGenome.pl",bed_path,"hg19",out,"-size 41 -p 4 -noweight -nlen 0 -bg",background_path)
system(cmd)

message("Success !")
