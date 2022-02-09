

source("scripts/utils/new_utils.R")
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/apelletier/HSPC_EpiStress/Data/tools/homer/bin/", sep = ":"))

out0<-"analyses/motif_analysis_dmc/"
cond<-"pval0.001_homer_with_background"

for(compa in c("C.L","MC.ML","FC.FL")){
  print(compa)
  res<-fread(paste0("analyses/model14_without_iugr/2020-04-16_res_locis_in_",compa,"_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv"),
           select = c(6,7,1,2:4),col.names = c("chr","pos","cpg_id","meth.change","pval","pval_adj"),dec = ",")
  
  out<-fp(out0,paste0(compa,cond))
  message("saving in",out)
  
  
  res[,start:=pos-20][start<1,start:=1][,end:=pos+20]
  res<-res[!is.na(chr)&!is.na(pos)]
  
  n_cpgs<-nrow(res[pval<0.001&abs(meth.change)>30])
  
  message("analyzing motif enrichment in ",n_cpgs,"DMCs versus 15k random cpgs background")
  
  bed_path<-fp(out0,paste0(compa,"40bp_win_dmc_pvalnom0.001_meth.change30"))
  background_path<-fp(out0,"background_cpgs_40bp_win_random15k.bed")
  
  fwrite(res[pval<0.001&abs(meth.change)>30][,.(chr,start,end)][order(chr,start)],
       bed_path,sep="\t",
       col.names=F)
  
  cmd<-paste("findMotifsGenome.pl",bed_path,"hg19",out,"-size 40 -p 4 -bg",background_path)
  system(cmd)
}



