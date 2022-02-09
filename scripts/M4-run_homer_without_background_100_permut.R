

source("scripts/utils/new_utils.R")
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/home/apelletier/HSPC_EpiStress/Data/tools/homer/bin/", sep = ":"))

n_perm<-100
compa<-"C.L"

out0<-ps("../methyl/analyses/motif_analysis_dmc/permut_",n_perm,compa)
message("motif analysis permutation in ",compa," (", n_perm,' times)')
dir.create(out0)

res_cpgs<-fread(ps("../methyl/analyses/model14_without_iugr/2020-04-16_res_locis_in_",compa,"_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv"),
select = c(6,7,1,2:4),col.names = c("chr","pos","cpg_id","meth.change","pval","pval_adj"),dec = ",")
  


cond<-"pval0.001_homer_without_background"

res_cpgs[,start:=pos-20][start<1,start:=1][,end:=pos+20]
res_cpgs<-res_cpgs[!is.na(chr)&!is.na(pos)]
n_cpgs<-nrow(res_cpgs[pval<0.001&abs(meth.change)>30])

message("analyzing motif enrichment in ",n_cpgs," permuted cpgs")
  

for(i in 1:n_perm){
  message("perm ",i,"/",n_perm)
  out<-fp(out0,paste0(i,cond))
  dir.create(out)
  
  bed_path<-fp(out,paste0("random",i,".bed"))
  
  random_cpgs<-res_cpgs[cpg_id%in%sample(cpg_id,n_cpgs),.(chr,start,end)][order(chr,start)]

  fwrite(random_cpgs[,.(chr,start,end)][order(chr,start)],
       bed_path,sep="\t",
       col.names=F)
  
  cmd<-paste("findMotifsGenome.pl",bed_path,"hg19",out,"-size 40 -p 4")
  system(cmd)
  
  res<-fread(fp(out,"knownResults.txt"),
             select = c(1,3),col.names = c("motif_name","pval"))
  
  res[,permut:=i]
  
  if(i==1){
    res_all<-copy(res)
  }else{
    res_all<-rbind(res_all,res)
  }
  
  if(i>5){
    message("removing res")
    unlink(out,recursive = T)
  }
}

message("calculating pval permut")

res_obs<-fread("../methyl/analyses/motif_analysis_dmc/C.Lpval0.001_homer_without_background/knownResults.txt",
             select = c(1,3),col.names = c("motif_name","pval"))

res_obs[,permut:=0]
res_all<-rbind(res_all,res_obs)
res_all[,obs:=permut==0]
res_all[,pval.perm:=sum(pval[obs==T]>=pval[obs==F])/n_perm,by="motif_name"]

fwrite(res_all,fp(out0,"res_permut.tsv"),sep='\t')

fwrite(res_all[obs==T,.(motif_name,pval,pval.perm)],fp(out0,"res_permut.tsv"),sep='\t')

