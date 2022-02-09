#motif analysis in in DMC
source("scripts/utils/new_utils.R")
out<-"analyses/motif_analysis_dmc/"
dir.create(out)

#which tool choosed ?

#1) meme-chip [options] [-db <motif database>]*
#docs : http://gensoft.pasteur.fr/docs/meme/5.0.4/meme-chip.html

dir.create(out)

res_cpgs<-fread("analyses/model14_without_iugr/2020-04-16_res_locis_in_C.L_pval_1_locisF.msp1.NA.fullMethyl.confScore.nbMethylNonZeros_model_14_.csv",
           select = c(6,7,1,2:4),col.names = c("chr","pos","cpg_id","meth.change","pval","pval_adj"),dec = ",")

res_cpgs[pval_adj<0.1&abs(meth.change)>30]
res_cpgs[,start:=pos-20][start<1,start:=1][,end:=pos+20]
res_cpgs<-res_cpgs[!is.na(chr)&!is.na(pos)]
fwrite(res_cpgs[pval_adj<0.1&abs(meth.change)>30][,.(chr,start,end)][order(chr,start)],
       file.path(out,"40bp_win_dmc_lga_vs_ctrl_padj0.1_meth.change30.bed"),sep="\t",
       col.names=F)


res_cpgs[pval<0.001&abs(meth.change)>30]
fwrite(res_cpgs[pval<0.001&abs(meth.change)>30][,.(chr,start,end)][order(chr,start)],
       file.path(out,"40bp_win_dmc_lga_vs_ctrl_pvalnom0.001_meth.change30.bed"),sep="\t",
       col.names=F)

random_cpgs<-res_cpgs[cpg_id%in%sample(cpg_id,15000),.(chr,start,end)][order(chr,start)]

fwrite(random_cpgs,
       file.path(out,"background_cpgs_40bp_win_random15k.bed"),sep="\t",
       col.names=F)

#get fasta
system("bedtools getfasta  -fi ref/hg19/ucsc.hg19.fasta -bed analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_padj0.1_meth.change30.bed -fo analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_padj0.1_meth.change30.fasta")
system("bedtools getfasta  -fi ref/hg19/ucsc.hg19.fasta -bed analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_pvalnom0.001_meth.change30.bed -fo analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_pvalnom0.001_meth.change30.fasta")
system("bedtools getfasta  -fi ref/hg19/ucsc.hg19.fasta -bed analyses/motif_analysis_dmc/40bp_win_cpgs_background_random5000.bed -fo analyses/motif_analysis_dmc/40bp_win_cpgs_background_random5000.fasta")


#simple motif analysis on methyl : HOMER puis FIMO
# Finding Enriched Motifs in Genomic Regions 
#use: findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
#By default this will perform de novo motif discovery as well as check the enrichment of known motifs
#all parameter in http://homer.ucsd.edu/homer/ngs/peakMotifs.html

#[in shell]
#find motif in pval<0.001 DMC region  ctrl vs lga with auto background : 
findMotifsGenome.pl analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_pvalnom0.001_meth.change30.bed hg19 analyses/motif_analysis_dmc/ctrl_vs_lga_pval0.001_homer/ -size 40

#with background : 
findMotifsGenome.pl analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_pvalnom0.001_meth.change30.bed hg19 analyses/motif_analysis_dmc/ctrl_vs_lga_pval0.001_homer_with_bg/ -size 40 -bg analyses/motif_analysis_dmc/40bp_win_cpgs_background_random10k.bed -p 4

#with background and mask
findMotifsGenome.pl analyses/motif_analysis_dmc/40bp_win_dmc_lga_vs_ctrl_pvalnom0.001_meth.change30.bed hg19 analyses/motif_analysis_dmc/ctrl_vs_lga_pval0.001_homer_with_bg_and_mask/ -size 40 -bg analyses/motif_analysis_dmc/40bp_win_cpgs_background_random10k.bed -mask -p 4

#CTRL F vs LGA F, and CTRL M vs LGA M
#run M1


#methyl TF, permut 3000 random CpGs => Stat3, ETS.. really Stress spe ?
#run M4
source("scripts/utils/new_utils.R")
out<-"analyses/motif_analysis_dmc/permut_100"

res<-fread("analyses/motif_analysis_dmc/permut_100C.L/res_permut.tsv")

res[pval.perm<0.1][order(pval.perm)]
res[,full_motif_name:=motif_name]
res[,motif_name:=sapply(full_motif_name,function(x)strsplit(x,"/")[[1]][1])]
res

ggplot(res[pval.perm<=0.1&pval<=0.05][order(pval.perm)])+geom_col(aes(x=motif_name,y=-log10(pval),fill=pval.perm))+
  scale_x_discrete(limits=res[pval.perm<=0.1&pval<=0.05][order(pval.perm)]$motif_name)+
  theme(axis.text.x = element_text(angle = 60, hjust=1))

res[str_detect(motif_name,"Stat")] #not sig



