source("scripts/utils/new_utils.R")
set.seed(1234)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/apelletier/HSPC_EpiStress/Data/tools/homer/bin/", sep = ":"))

out<-"outputs/03B-motif_analysis/with_auto_background"
dir.create(out,recursive = T)

bed_path<-fp(out,"40bp_win_dmc_pvalnom0.001_meth.change30")
background_path<-fp(out,"15k_randoms_cpgs_40bp.bed")


res<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
  
cpgs<-fread("datasets/cd34/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",
            select = c("id","chr","start"),
            col.names = c("cpg_id","chr","pos"))
res<-merge(res,cpgs)
res[,start:=pos-20][start<1,start:=1][,end:=pos+20]
res<-res[!is.na(chr)&!is.na(pos)]


random_cpgs<-res[cpg_id%in%sample(cpg_id,15000),.(chr,start,end)][order(chr,start)]

fwrite(random_cpgs,
       background_path,sep="\t",
       col.names=F)

n_cpgs<-nrow(res[P.Value<0.001&abs(logFC)>30]) #4774

message("analyzing motif enrichment in ",n_cpgs,"DMCs versus autobackground")


fwrite(res[P.Value<0.001&abs(logFC)>30][,.(chr,start,end)][order(chr,start)],
     bed_path,sep="\t",
     col.names=F)

cmd<-paste("findMotifsGenome.pl",bed_path,"hg19",out,"-size 40 -p 4")
system(cmd)

