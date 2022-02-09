#Concept : Find comethylated genes across samples 

#postulats : 
  #des TFs peuvent réguler l'activité de la DNMT/TET autour des sites/motifs sur lesquelles ils se fixent (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6923324/)

#question biologique : y a til des groupes (modules) de genes cométhylés a travers les échantillons, qui seraient la signature de l'action d'un TF ?
#question subsidiaire : ces modules sont ils associés avec un trait phénotypique particulier ?
# ==> question technique : Comment identifier des modules de genes co-méthylés ?

#PLAN DANALYSE
  #main steps (Verifier/Valider chaque étape avant de passer à la suivante)
    #1) aggréger la methylation des CpGs à l'échelle du gène (genescore)
      #a) lier les CpGs au genes
      #b) assigner un poids pour chaque lien CpG-gene
      #c) pour chaque gène, faire la moyenne pondéré des valeurs de méthylation des CpGs liés à lui

    #2) cross correlation de ce genescore à travers les samples
      #a) caractérisation du score pour choisir le bon test de correlation
      #b) analyse de correlation
      #c) clusteriser la matrice de correlation pour identifier des modules

    #3) caractérisations/validation des modules 
      #a) pathway/Biological process/TF signature enrichment analysis in these modules
      #b) TF binding motifs analysis 
      #c) correlation avec des traits phénotypiques (Birth weight, maternal age..)
      #d) Integration avec les modules de co-expression identifies en scRNA-seq




#ANALYSIS
source("scripts/utils/new_utils.R")
out<-"../methyl/outputs/04-comethylation_analysis"
dir.create(out)

    #1) aggréger la methylation des CpGs à l'échelle du gène (genescore)
      #a) lier les CpGs au genes 
#see 02A 
      #b) assigner un poids pour chaque lien CpG-gene
#see 02B
cpgs_weight<-fread("../methyl/outputs/02A-CpGs_annotations/cpgs_genes_annot_and_weight.csv.gz")
cpgs_weight[,cpg_weight:=regul_weight+links_weight]

      #c) pour chaque gène, faire la moyenne pondéré des valeurs de méthylation des CpGs liés à lui
mtd<-fread("datasets/cd34/metadata_cl_pcs_040521.csv")
meth<-fread("datasets/cd34/meth_data_filtered.csv.gz",
            select = c("cpg_id",mtd$sample))
meth_df<-melt(meth,id.vars = 1,
              variable.name = "sample",
              value.name = "methylation")
meth_df<-merge(meth_df,cpgs_weight[,.(cpg_id,gene,cpg_weight)],by="cpg_id")
meth_df[,gene_meth:=sum(methylation*cpg_weight)/sum(cpg_weight),by=c("sample","gene")]
fwrite(meth_df,fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")
meth_df<-fread(fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")
rm(meth)
#caracterisation score : 
methg<-unique(meth_df,by=c("sample","gene"))

rm(meth_df)
plot(density(methg$gene_meth))

#tecnic bias by sample ?
ggplot(methg)+geom_density(aes(gene_meth))+facet_wrap("sample") #yes

#QC : rm genes with 0 values
methg0<-methg[!(gene%in%gene[gene_meth==0])]
length(unique(methg0$gene)) #6224

plot(density(methg0$gene_meth))
#tecnic bias by sample ?
ggplot(methg0)+geom_density(aes(gene_meth))+facet_wrap("sample") #yes

methg0[,gene_meth_scaled:=scale(gene_meth),by=c("sample")]
ggplot(methg0)+geom_density(aes(gene_meth_scaled))+facet_wrap("sample") #ok

    #2) cross correlation de ce genescore à travers les samples
      #a) caractérisation du score pour choisir le bon test de correlation
#between 0-100 (-4 and 8 for the scaled version), 
#moyenne pondéré of the percentage of HpA2 (methylation sensitive) cleaving counts compared to Msp1 (reference, unsensitive) cleaving counts
#=> based on https://statsandr.com/blog/correlation-coefficient-and-correlation-test-in-r/, use spearman correl because :
#"Spearman correlation (which is actually similar to Pearson but based on the ranked values for each variable rather than on the raw data) is often used to evaluate relationships involving qualitative ordinal variables or quantitative variables if the link is partially linear"

      #b) analyse de correlation with spearman
# meth based on neurons jc paper
#use MEGENA ()
#-test for planarity with PMFG 
renv::install("songw01/MEGENA")

library(MEGENA)
# input parameters
n.cores <- 10; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "spearman" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 100; # number of permutations for calculating connectivity significance p-value. 
# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
###########

methg0_mat<-as.matrix(dcast(methg0,sample~gene,value.var = "gene_meth_scaled"),rownames = "sample")
dim(methg0_mat) #70 6224
head(t(methg0_mat))
ijw <- calculate.correlation(t(methg0_mat),
                             doPerm = cor.perm,
                             output.corTable = FALSE,
                             output.permFDR = FALSE,
                             method=method)

head(ijw)

# calculate PFN
#In this step, Planar Filtered Network (PFN) is calculated by taking significant correlation pairs, ijw. In the case of utilizing a different similarity measure, one can independently format the results into 3-column data frame with column names c("row","col","weight"), and make sure the weight column ranges within 0 to 1. Using this as an input to calculate.PFN() will work just as fine. 
#### register multiple cores if needed: note that set.parallel.backend() is deprecated. 
run.par = doPar & (getDoParWorkers() == 1) 
if (run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  # check how many workers are there
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
}
##### calculate PFN
el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores,keep.track = FALSE)
g <- graph.data.frame(el,directed = FALSE) #make the graph
vcount(g) #6224 nodes
fwrite(el,fp(out,"graph.csv.gz"))
saveRDS(g,fp(out,"graph.rds"))

# perform clustering
#MCA clustering is performed to identify multiscale clustering analysis. "MEGENA.output" is the core output to be used in the down-stream analyses for summarization and plotting.
##### perform MCA clustering.
megena_out <- do.MEGENA(g,
 mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE, #for hub analysis within significant modules
 min.size = 10,max.size = vcount(g)/2, #max module size = 3112
 doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
 save.output = FALSE)
###### unregister cores as these are not needed anymore.
if (getDoParWorkers() > 1)
{
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
saveRDS(megena_out,fp(out,"megena_spearman_pfn_mca_outputs.rds"))

summary_out <- MEGENA.ModuleSummary(megena_out,
	mod.pvalue = module.pval,hub.pvalue = hub.pval,
	min.size = 10,max.size = vcount(g)/2,
	annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
	output.sig = TRUE)
saveRDS(summary_out,fp(out,"summary_megena_outputs.rds"))

#some signif modules ?
modules<-megena_out$module.output$modules
modules_pval<-megena_out$module.output$module.pvalue
modules_sig<-modules[modules_pval<0.05]
length(modules_sig) #449

#heatmap module 
source("scripts/utils/new_utils.R")
out<-"../methyl/outputs/04-comethylation_analysis"

meth_df<-fread(fp(out,"methylation_sum_by_gene_datatable.tsv.gz"),sep="\t")

methg0<-methg[!(gene%in%gene[gene_meth==0])]

methg0_mat<-as.matrix(dcast(methg0,sample~gene,value.var = "gene_meth_scaled"),rownames = "sample")



    #3) caractérisations/validation des modules 
source("scripts/utils/new_utils.R")
out<-"../methyl/outputs/04-comethylation_analysis"
networks<-readRDS(fp(out,"megena_spearman_pfn_mca_outputs.rds"))


modules<-networks$module.output$modules
modules_pval<-networks$module.output$module.pvalue
modules_sig<-modules[modules_pval<0.05]
modules_sig_10<-modules_sig[sapply(modules_sig,function(x)length(x)>10)]
length(modules_sig_10) #152

      #a) pathway/Biological process/TF signature enrichment analysis in these modules
#kegg
modules_kegg<-lapply(modules_sig_10,function(x){
  
  print(length(x))
  genes_id<-bitr(x,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
  return(as.data.frame(enrichKEGG(genes_id,
                 organism = "hsa",pvalueCutoff = 0.05)))
  })

modules_kegg_dt<-Reduce(rbind,lapply(names(modules_kegg),function(mod_name)data.table(modules_kegg[[mod_name]])[,module:=mod_name][,module.size:=length(modules[[mod_name]])]))
table(modules_kegg_dt$module)
# c1_100  c1_11 c1_120 c1_129 c1_132 c1_168 c1_186 c1_209  c1_21  c1_22 c1_255 
#     26      1      2      4      1      1      4      7      5      8      2 
# c1_258  c1_27 c1_281 c1_283 c1_299 c1_311 c1_321 c1_341 c1_354 c1_373  c1_39 
#      6      5      1      1      2      1      1      3      1      1      1 
#   c1_4 c1_425  c1_63  c1_64  c1_66  c1_74  c1_91 
#      1      6      1      1      2      1      6 

summary(unique(modules_kegg_dt,by="module")$module.size)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 11.00   12.00   18.00   57.79   25.00 1099.00 

modules_kegg_dt[module=="c1_258"]

modules_kegg_dt[,n.term.enriched:=.N,by="module"]
modules_kegg_dt[,top5:=p.adjust<=sort(p.adjust)[5],by="module"]
split(modules_kegg_dt[top5==T,.(module,Description,Count)],by="module")
fwrite(modules_kegg_dt,fp(out,"res_kegg_on_modules_sig_10genes.csv.gz"))


#go
modules_go<-lapply(modules_sig_10,function(x){
  
  print(length(x))
  genes_id<-bitr(x,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID
  return(as.data.frame(enrichGO(genes_id,
                                OrgDb = org.Hs.eg.db,
                                pvalueCutoff = 0.05)))
  })

modules_go_dt<-Reduce(rbind,lapply(names(modules_go),function(mod_name)data.table(modules_go[[mod_name]])[,module:=mod_name][,module.size:=length(modules[[mod_name]])]))
table(modules_go_dt$module)
# c1_10 c1_100  c1_11 c1_129 c1_130 c1_132 c1_140 c1_156 c1_158 c1_184 c1_209 c1_229 c1_256 
#     14     28     12      7     22      1     16      2      7      1     13      1     16 
# c1_258 c1_264 c1_266 c1_281 c1_283  c1_29 c1_294 c1_310 c1_311 c1_313 c1_333 c1_334 c1_342 
#      1      8      1      2     21      4      3     17     11     17      3     30      2 
# c1_355 c1_359 c1_368 c1_373 c1_391   c1_4 c1_425  c1_43 c1_457  c1_48  c1_63  c1_64  c1_74 
#      9     20     14      7      1      1     12      1      1      4      8      1      6 
#  c1_81  c1_91 
#      4      3 

modules_go_dt[,n.term.enriched:=.N,by="module"]
modules_go_dt[Count>=5] #only 1

modules_go_dt[,top5:=p.adjust<=sort(p.adjust)[5],by="module"]
split(modules_go_dt[top5==T,.(module,Description,Count,geneID)],by="module")
fwrite(modules_go_dt,fp(out,"res_go_on_modules_sig_10genes.csv.gz"))


#TF
TFsign<-read.gmt("ref/tf_signatures/TRRUST_Transcription_Factors_2019.txt")


modules_tf<-lapply(modules_sig_10,function(x){
  
  print(length(x))
  return(as.data.frame(enricher(x,TERM2GENE= TFsign,
                               pvalueCutoff = 0.05)))
  })
modules_tf_dt<-Reduce(rbind,lapply(names(modules_tf)[sapply(modules_tf, function(x)ncol(x)>0)],
                                         function(mod_name)data.table(modules_tf[[mod_name]])[,module:=mod_name][,module.size:=length(modules[[mod_name]])]))

summary(unique(modules_tf_dt,by="module")$module.size)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #  11.00   13.00   19.00   23.28   27.00   66.00 
modules_tf_dt[,n.term.enriched:=.N,by="module"]
modules_tf_dt[,top5:=p.adjust<=sort(p.adjust)[5],by="module"]
split(modules_tf_dt[top5==T,.(module,Description,Count,geneID)],by="module")

fwrite(modules_tf_dt,fp(out,"res_tf_on_modules_sig_10genes.csv.gz"))



#merge the res

      #b) TF binding motifs analysis 
      #c) correlation avec des traits phénotypiques (Birth weight, maternal age..)
      #d) Integration avec les modules de co-expression identifies en scRNA-seq
