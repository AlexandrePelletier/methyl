#2020-12-08 test cbps with new genescore
library(data.table)
library(stringr)
library(ggplot2)
source("scripts/utils/new_utils.R")

out<-"analyses/2020-12-08_cbps_new_genescore"
dir.create(out)

#[recalc limma res] with same than before (because not save all limma res with locisID.. )
library(limma)
methyl<-fread("datasets/cd34/2020-05-25_methyl_data_before_limma.csv")
methyl<-methyl[,cpgID:=locisID][,.SD,.SDcols=c(ncol(methyl),1:(ncol(methyl)-1))][,-"locisID"]

meta<-fread("datasets/cd34/cleaned_batch_CD34_library_date_220620.csv")
meta
samples<-colnames(methyl)[colnames(methyl)%in%meta[Group_name%in%c("C","L")]$sample]
methyl<-methyl[,.SD,.SDcols=c("cpgID",samples)]
dim(methyl) #70 samples
meta<-meta[sample%in%colnames(methyl)]
nrow(meta) #70





var_fac_to_model<-c("Group_Sex","latino","Group_Complexity_Fac","batch")
var_to_model<-c(var_fac_to_model,"Mat.Age")
meta[,(var_fac_to_model):=lapply(.SD, as.factor),.SDcols=var_fac_to_model]

metaF<-data.frame(meta[rowSums(is.na(meta[,.SD,.SDcols=var_to_model]))==0],row.names = "sample")
head(metaF)
formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac +batch

design<-model.matrix(formule,data = metaF)

fit <- lmFit(data.frame(methyl,row.names = "cpgID")[,row.names(design)], design)

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
results <- decideTests(fit2)
sum(abs(results)) #8>7
colSums(abs(results))
sum(apply(fit2$p.value<0.001,1,any)) #11561
compas<-colnames(fit2$p.value)

res<-data.table(topTable(fit2,coef = "CF.LF",n = Inf),keep.rownames = "cpgID")
fwrite(res,"analyses/model14_without_iugr/2020-12-15_res_limma_CF.LF.csv")

#check if concordant with res_old GS
res_o<-fread("analyses/model14_without_iugr/2020-09-16_res_CF.LF_with_GeneScore_and_permut.csv")
res[cpgID==552]
res_o[cpgID==552]

res[,cpgID:=as.integer(cpgID)]
res_m<-merge(res,res_o,by="cpgID")
res_m[logFC>10&meth.change<(-10)] 

plot(res_m$logFC,res_m$meth.change)
cor(res_m$logFC,res_m$meth.change)
plot(-log10(res_m$P.Value),-log10(res_m$pval))
cor(-log10(res_m$P.Value),-log10(res_m$pval))

summary(res_m[P.Value<0.001]$pval) #small correlation
summary(res_m$pval) 

#this big pvalue and methcange change, is due to limma or to a script error ?
#=> redo limma analysis
#10 times for CF.CL

#prep data
var_fac_to_model<-c("Group_Sex","latino","Group_Complexity_Fac","batch")
meta[,(var_fac_to_model):=lapply(.SD, as.factor),.SDcols=var_fac_to_model]#factorize factor

meta<-meta[rowSums(is.na(meta[,.SD,.SDcols=var_to_model]))==0] #rm samples with NA
meta_df<-data.frame(meta,row.names = "sample")
head(meta_df)

formule<- ~0 + Group_Sex  + latino + Mat.Age + Group_Complexity_Fac +batch
design<-model.matrix(formule,data = meta_df)

methyl_df<-data.frame(methyl,row.names = "cpgID")[,row.names(design)]#prepare matrix
head(methyl_df)

for(i in 1:10){
  print(i)
  fit <- lmFit(methyl_df,design )

  cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                             F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                             CF.CM="Group_SexC_F-Group_SexC_M",
                             CF.LF="Group_SexC_F-Group_SexL_F",
                             CF.LM="Group_SexC_F-Group_SexL_M",
                             CM.LM="Group_SexC_M-Group_SexL_M",
                             CM.LF="Group_SexC_M-Group_SexL_F",
                             LM.LF="Group_SexL_M-Group_SexL_F",
                             levels=design)


  fit2  <- contrasts.fit(fit, cont.matrix)
  fit2  <- eBayes(fit2)
  res<-data.table(topTable(fit2,coef = "CF.LF",n = Inf),keep.rownames = "cpg_id")
  res[,run_id:=i]
  if(i==1){
    res_all<-copy(res)
  }else{
    res_all<-rbind(res_all,res)
  }
}
res_all[,pval.mean:=mean(P.Value),by='cpg_id']
res_all[,pval.sd:=sd(P.Value),by='cpg_id']

res_all[cpg_id==552] #same each time
res_all[pval.sd!=0] #no
summary(res_all$pval.sd)

#maybe because not in the good order, but anyway, need redo all : 
meta_df<-data.frame(meta,row.names = "sample")
design<-model.matrix(formule,data = meta_df)
methyl_df<-data.frame(methyl,row.names = "cpgID")[,row.names(design)]#prepare matrix
fit <- lmFit(methyl_df,design )

cont.matrix <- makeContrasts(C.L = "(Group_SexC_F+Group_SexC_M)-(Group_SexL_F+Group_SexL_M)",
                           F.M="(Group_SexC_F+Group_SexL_F)-(Group_SexC_M+Group_SexL_M)",
                           CF.CM="Group_SexC_F-Group_SexC_M",
                           CF.LF="Group_SexC_F-Group_SexL_F",
                           CF.LM="Group_SexC_F-Group_SexL_M",
                           CM.LM="Group_SexC_M-Group_SexL_M",
                           CM.LF="Group_SexC_M-Group_SexL_F",
                           LM.LF="Group_SexL_M-Group_SexL_F",
                           levels=design)


fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)
compas<-colnames(fit2$p.value)
res_list<-lapply(compas,function(x)data.table(topTable(fit2,coef = x,n = Inf,sort.by = "none"),keep.rownames = "cpg_id"))
names(res_list)<-compas
saveRDS(list(limma_obj=fit2,res=res_list),"analyses/model14_without_iugr/2020-15-12_limma_obj_and_all_res.rds")

#recalc old GeneScore
cpgs_anno<-fread("ref/2020-06-29_All_CpG-Gene_links.csv")
cpgs_anno<-cpgs_anno[,cpg_id:=locisID][,-"locisID"]
res_anno_list<-lapply(res_list,function(x){
  x[,cpg_id:=as.integer(cpg_id)]
  x[,meth.change:=logFC][,pval:=P.Value][,pval_adj:=adj.P.Val]
  return(merge(x[,.(cpg_id,meth.change,pval,pval_adj)],cpgs_anno,by="cpg_id"))
})
res_anno_list
res_anno_list<-lapply(res_anno_list,function(x){
  return(x[,cpg_score:=(-log10(pval)/4*meth.change)*RegWeight*LinksWeight])
})

res_anno_list<-lapply(res_anno_list,function(x){
  x[,n_cpg.gene:=.N,by=.(gene)]
  x[,n_cpg_sig.gene:=sum(pval<0.001),by=.(gene)]
  x[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
  x[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]
  
  
})

res_anno_list<-lapply(names(res_anno_list),function(comp)res_anno_list[[comp]][,compa:=comp])
res_all<-Reduce(rbind,res_anno_list)
fwrite(res_all,"analyses/model14_without_iugr/2020-12-15_all_res_genescore1.csv",sep=";")
res_all[compa=="CF.LF"&cpg_id==552]
#[end recalc limma res]

res<-res_all[compa=="CF.LF"]
res[,chromatin_feature:=sapply(chromatin_state,function(x)ifelse(x==6,"Promoter",
                                                                   ifelse(x==4,"Active Enhancer",
                                                                          ifelse(x==5,"Inactive Enhancer",
                                                                                 ifelse(x%in%1:3,"Gene Body","Heterochromatin")))))]


res

#ANNOT CPGs
cpgs<-res[,.(chr,pos,cpgID)][!is.na(chr)]

#1)cpgs_genes_links
#by distance to tss
#4 first genes within +/-200kb


genes_in_200kb<-bed_inter(cpgs[,start:=pos-200e3][start<0,start:=0][,end:=pos+200e3][,.(chr,start,end,cpgID)][order(chr,start)],
          "ref/hg19/tss_genes.bed",
          select = c(5,6,8,9,4),col.names = c("chr","tss_pos","gene","strand","cpgID"))

genes_in_200kb
cpgs_genes_tss<-merge(cpgs,genes_in_200kb,all.x=T)[,.(chr,pos,cpgID,gene,tss_pos,strand)]
cpgs_genes_tss[strand=="+",tss_dist:=pos-tss_pos][strand=="-",tss_dist:=tss_pos-pos]
cpgs_genes_tss[,ngene.cpg:=length(unique(gene)),by="cpgID"]
summary(cpgs_genes_tss$ngene.cpg)
summary(abs(cpgs_genes_tss$tss_dist))
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #   1.00    6.00   11.00   13.06   18.00  100.00 

cpgs_genes_4<-cpgs_genes_tss[order(cpgID,abs(tss_dist))][,top4:=gene%in%unique(gene)[1:4],by="cpgID"][top4==T][,-"top4"]
summary(abs(cpgs_genes_4$tss_dist))
#for >1 cpg_gene link (because multiple TSS), keep most closest cpg_gene link
cpgs_genes_4<-unique(cpgs_genes_4[order(cpgID,gene,abs(tss_dist))],by=c("cpgID","gene")) 
summary(abs(cpgs_genes_4$tss_dist))
cpgs_genes_4

rm(cpgs_genes_tss)

fwrite(cpgs_genes_4,fp(out,"cpgs_4closest_gene_tss_linked_within_200kb_around.tsv"),sep="\t")

#by presence in eQTL region
cpgs_eQTR<-bed_inter(a=cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpgID)][order(chr,start)],
          b="ref/eQTL/whole_blood_and_tissue_wide_eQTR_hgnc_gene_hg19.bed",
          select = c(4,8),col.names = c("cpgID","gene"))

cpgs_eQTR #887k  cpgs - snp+/-500pb match
unique(cpgs_eQTR) #554k cpg-gene match
unique(cpgs_eQTR,by="cpgID")#246k cpgs linked to a gene
unique(cpgs_eQTR,by="gene")#13k genes linked to a cpg
cpgs_eQTR[,in_eQTR:=T]

#merge the 2 links
cpgs_genes_4[,in_eQTR:=F]
cpgs_genes<-merge(cpgs_genes_4[,.(cpgID,gene,in_eQTR,tss_dist)],cpgs_eQTR,all=T)
cpgs_genes[in_eQTR==T]

#merge with coord
cpgs_genes<-merge(cpgs[,.(cpgID,chr,pos)],cpgs_genes,all.x = T,by="cpgID")

fwrite(cpgs_genes,fp(out,"cpgs_4closest_gene_tss_linked_within_200kb_around_and_eQTR_linked_genes.tsv"),sep="\t")

#2)cpgs_regulatory region
#ensembl regulatory reg matching


cpgs_ensembl<-bed_inter(a=unique(cpgs,by="cpgID")[,chr:=str_remove(chr,'chr')][,start:=pos][,end:=pos+1][,.(chr,start,end,cpgID)][order(chr,start)],
          b="ref/ensembl_regulatory/ensembl_regulatory_hg19.bed",
          select = c(4,8),col.names = c("cpgID","ensembl_regulatory_domain"))



cpgs_ensembl[,ensembl_regulatory_domain:=paste(ensembl_regulatory_domain,collapse = "/"),by="cpgID"] 
cpgs_ensembl<-unique(cpgs_ensembl)
cpgs_ensembl

cpgs_chromatin<-res[,.(cpgID,chr,pos,chromatin_state,chromatin_feature)]

cpgs_reg<-merge(cpgs_chromatin,cpgs_ensembl,all.x=T,by="cpgID")

fwrite(cpgs_reg,fp(out,"cpgs_annot_ensembl_regulatory_domain_and_chromHMM_chromatin_features.tsv"),sep="\t")

cpgs_anno<-merge(cpgs_reg,cpgs_genes,all=T,by="cpgID")

fwrite(unique(cpgs_anno[order(cpgID,gene)][,.(cpgID,chr,pos,in_eQTR,tss_dist,gene,chromatin_state,chromatin_feature,ensembl_regulatory_domain)]),
       fp(out,"cpgs_annot.tsv"),sep="\t")


# Calculate GeneScore

#1) calculate cpg weight
cpgs_genes<-cpgs_anno[!is.na(gene)&gene!=""]
unique(cpgs_genes,by="cpgID") #741408/~800k cpgs link to a gene

unique(cpgs_genes[in_eQTR==T],by="cpgID") #dont 246k linked thx to eQTR
cpgs_genes[,double.linked:=any(in_eQTR==T)&any(in_eQTR==F),by=c("cpgID","gene")]
unique(cpgs_genes[double.linked==T],by="cpgID")#dont 63504 also gene linked based on tss
#=> + ~180k cpgs gene link thx to eQTR

# a) linksWeight
cpgs_genes[in_eQTR==F,links_score:=sapply(abs(tss_dist),function(x){
  if(x<1000)return(1)
  else if(x<20000)return(0.5+0.5*sqrt(1000/x))
  else return(0.5*sqrt(20000/x))
    })]

cpgs_genes[in_eQTR==T,links_score:=1]

cpgs_genes[,links_weight:=max(links_score),by=c("cpgID","gene")]


#b) Regulatory Weights
cpgs_reg<-unique(cpgs_genes,by=c("cpgID"))
unique(cpgs_reg$chromatin_feature)
cpgs_reg[,chromatin_score:=sapply(chromatin_feature,function(x){
      if(x%in%c("Promoter","Active Enhancer"))return(1)
      else if(x=="Inactive Enhancer")return(0.75)
      else if(x%in%"Gene Body")return(0.5)
      else return(0)
    })]


unique(cpgs_reg$chromatin_feature)
cpgs_reg[,ensembl_reg_score:=sapply(ensembl_regulatory_domain,function(x){
    vecX<-strsplit(as.character(x),"/")[[1]]
    if(any(c("CTCF_binding_site","promoter","enhancer")%in%vecX)){
      score<-0.5
    }else if(any(c("open_chromatin_region","promoter_flanking_region")%in%vecX)){
      score<-0.25
    }else{
      score<-0
    }
    if("TF_binding_site"%in%vecX){
      score<-score+0.5
    }
    return(score)
  })]

cpgs_reg[,regul_weight:=(0.5+1.5*((chromatin_score+ensembl_reg_score)/2))]

cpgs_score<-merge(cpgs_genes,cpgs_reg[,.(cpgID,chromatin_score,ensembl_reg_score,regul_weight)],by="cpgID")
fwrite(cpgs_score,fp(out,"cpgs_genes_annot_and_weight_reference.tsv"),sep="\t")

#2) merge with res and calculate CpG Score
res<-merge(res[,.(cpgID,pval,meth.change)],unique(cpgs_score[order(cpgID,gene,-links_score)],by=c("cpgID","gene")),by=c("cpgID"))
res[,cpg_score:=(-log10(pval)/4*meth.change)*regul_weight*links_weight]

#3) GeneScore

res[,n_cpg.gene:=.N,by=.(gene)]
res[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^(1/4),by="gene"]
res[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

plot(density(unique(res,by="gene")$n_cpg.gene))
summary(unique(res,by="gene")$n_cpg.gene)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 #    1.0    50.0    84.0   108.2   139.0   963.0  


plot(density(unique(res,by="gene")$gene_score))

ggplot(unique(res,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()

res[order(-abs(gene_score),-abs(cpg_score))][pval<0.05][gene_score>800]

plot(density(log10(abs(unique(res,by=c("gene"))$gene_score))))


unique(res[gene_score>400]$gene)

fwrite(res,fp(out,"res_gene_score_lgaF_vs_ctrlF.csv"),sep=";")

#VALIDATION #NOOOP
library(data.table)
library(stringr)
#WHy ??
#loose CpG score interest ?
#correl btw old and new
res_n<-fread("analyses/2020-12-08_cbps_new_genescore/res_gene_score_lgaF_vs_ctrlF.csv")
res_n<-res_n[,cpg_id:=cpgID][,-"cpgID"]

res_o<-fread("analyses/model14_without_iugr/2020-12-15_all_res_genescore1.csv")[compa=="CF.LF"]

res_merge<-merge(res_n,res_o,by=c("cpg_id","gene"))

plot(res_merge$cpg_score.x,res_merge$cpg_score.y)#some difference
cor(res_merge$cpg_score.x,res_merge$cpg_score.y)^2 #0.96

#due to what ?
res_merge[,cpg_score.dt:=abs(cpg_score.x-cpg_score.y)]
res_merge[,cpg_score.increase:=sign(cpg_score.x-cpg_score.y)==sign(cpg_score.y)]
summary(res_merge[,.(cpg_score.dt,cpg_score.increase)])
 #   cpg_score.dt      cpg_score.increase
 # Min.   : 0.00000   Mode :logical     
 # 1st Qu.: 0.00000   FALSE:242315      
 # Median : 0.00022   TRUE :195908      
 # Mean   : 0.32093                     
 # 3rd Qu.: 0.06674                     
 # Max.   :70.94492                     

#globally no significant change / increase

summary(res_merge[abs(cpg_score.x)>20|abs(cpg_score.y)>20][,.(cpg_score.dt,cpg_score.increase)])
 # cpg_score.dt    cpg_score.increase
 # Min.   : 0.000   Mode :logical     
 # 1st Qu.: 0.000   FALSE:7416        
 # Median : 0.000   TRUE :8126        
 # Mean   : 2.429                     
 # 3rd Qu.: 1.396                     
 # Max.   :70.945                     
#=> no big difference

#which are this cpgs that change ?
#that increase >15
res_merge[cpg_score.increase==T&cpg_score.dt>15] #only 368 cpgs, increase becuse in eQTR or close to TSS
#Note : 
# - dans old GS, LinksWeight mal programmé car tss_dist == -1 et LinksWeight pas maximal (car dans eQTR)
# - dans new GS, cpg 16156 pas dans eQTR alors qu'il y est dans old GS 


res_merge[cpg_score.increase==F&cpg_score.dt>15]#only 248

#now, optimize ncpg_weitgh

#that we have before :
summary(res_o$n_cpg.gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   20.00   37.00   54.74   69.00  514.00 
ggplot(unique(res_o,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()

degs<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
degs[is.na(padj),padj:=1]
degs[padj<0.05]
res_o_degs<-merge(res_o,degs[,.(gene,pvalue,padj,log2FoldChange)])
ggplot(res_o_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,140))



summary(res_n$n_cpg.gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.0    89.0   145.0   182.7   235.0   963.0 

res_n[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.65,by="gene"]
res_n[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

ggplot(unique(res_n,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
res_n_degs<-merge(unique(res_n,by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)])
ggplot(res_n_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,50))
#loose correl degs genescore even if  i rm ncpg influence => pas assez stringent.
#need decrease ncpg by gene, by 1) link tss dist, 2) eQTR admission

#TEST OPTIM NEW GS
#=> new genescore sans eqtr ?
res_n<-fread("analyses/2020-12-08_cbps_new_genescore/res_gene_score_lgaF_vs_ctrlF.csv")

cpgs_score<-fread("analyses/2020-12-08_cbps_new_genescore/cpgs_genes_annot_and_weight_reference.tsv")

res_n_sans_eqtr<-merge(unique(res_n[,.(cpgID,pval,meth.change)]),unique(cpgs_score[in_eQTR==F][order(cpgID,gene,-links_score)],by=c("cpgID","gene")),by=c("cpgID"))

res_n_sans_eqtr[,cpg_score:=(-log10(pval)/4*meth.change)*regul_weight*links_weight]

res_n_sans_eqtr[,n_cpg.gene:=.N,by=.(gene)]
res_n_sans_eqtr[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res_n_sans_eqtr[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.65,by="gene"]
res_n_sans_eqtr[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

ggplot(unique(res_n_sans_eqtr,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
degs<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
degs[is.na(padj),padj:=1]
res_n_sans_eqtr_degs<-merge(unique(res_n_sans_eqtr[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)],by="gene")

pnew<-ggplot(res_n_sans_eqtr_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,50))
pnew
#retrouve la correl ==> need to be more stringeant in eQTR def OR take it differently 
#(more susceptible to be enhancer, so influence positive of th methylation)

#compare to pold
res_o<-fread("analyses/model14_without_iugr/2020-12-15_all_res_genescore1.csv")[compa=="CF.LF"]
res_o_degs<-merge(unique(res_o[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)])
pold<-ggplot(res_o_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,100))
pnew+pold

wilcox.test(res_n_sans_eqtr_degs[padj<0.05]$gene_score,res_n_sans_eqtr_degs[padj>=0.05]$gene_score) #0.019
wilcox.test(res_o_degs[padj<0.05]$gene_score,res_o_degs[padj>=0.05]$gene_score) #0.002

#res_old est encore mieux. Car uniquement lié au gene le plus proche ?
library(data.table)
library(ggplot2)
library(patchwork)
res_n<-fread("analyses/2020-12-08_cbps_new_genescore/res_gene_score_lgaF_vs_ctrlF.csv")
cpgs_score<-fread("analyses/2020-12-08_cbps_new_genescore/cpgs_genes_annot_and_weight_reference.tsv")

res_n_sans_4genes<-merge(unique(res_n[,.(cpgID,pval,meth.change)]),
                       unique(cpgs_score[in_eQTR==F][,is.closest:=abs(tss_dist)==min(abs(tss_dist)),
                                                     by=.(cpgID)][is.closest==T][order(cpgID,gene,-links_score)],
                              by=c("cpgID","gene")),by=c("cpgID"))

res_n_sans_4genes[,cpg_score:=(-log10(pval)/4*meth.change)*regul_weight*links_weight]

res_n_sans_4genes[,n_cpg.gene:=.N,by=.(gene)]
res_n_sans_4genes[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res_n_sans_4genes[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.25,by="gene"]
res_n_sans_4genes[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

ggplot(unique(res_n_sans_4genes,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
degs<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
degs[is.na(padj),padj:=1]
res_n_sans_4genes_degs<-merge(unique(res_n_sans_4genes[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)],by="gene")

pnew<-ggplot(res_n_sans_4genes_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,100))
pnew #bonne correl

res_o<-fread("analyses/model14_without_iugr/2020-12-15_all_res_genescore1.csv")[compa=="CF.LF"]
res_o_degs<-merge(unique(res_o[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)])
wilcox.test(res_n_sans_4genes_degs[padj<0.05]$gene_score,res_n_sans_4genes_degs[padj>=0.05]$gene_score) #0.000276
wilcox.test(res_o_degs[padj<0.05]$gene_score,res_o_degs[padj>=0.05]$gene_score) #0.002
#=> new GS sans eqtr et avec closest gene, mieux que old pour pred

#def if 1) 2 genes tts-linked is best, 2) more stringeant eQTR if best ac a) influence comme tss b) influence +/- ~ localisation in enh/prom

#1) 2 genes tss-linked is best?
res_n<-fread("analyses/2020-12-08_cbps_new_genescore/res_gene_score_lgaF_vs_ctrlF.csv")
cpgs_score<-fread("analyses/2020-12-08_cbps_new_genescore/cpgs_genes_annot_and_weight_reference.tsv")

res_n_2genes<-merge(unique(res_n[,.(cpgID,pval,meth.change)]),
                       unique(cpgs_score[in_eQTR==F][,top2.closest:=abs(tss_dist)%in%sort(abs(tss_dist))[1:2],
                                                     by=.(cpgID)][top2.closest==T][order(cpgID,gene,-links_score)],
                              by=c("cpgID","gene")),by=c("cpgID"))

res_n_2genes[,cpg_score:=(-log10(pval)/4*meth.change)*regul_weight*links_weight]

res_n_2genes[,n_cpg.gene:=.N,by=.(gene)]
res_n_2genes[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res_n_2genes[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.25,by="gene"]
res_n_2genes[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]

ggplot(unique(res_n_2genes,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
degs<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
degs[is.na(padj),padj:=1]
res_n_2genes_degs<-merge(unique(res_n_2genes[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)],by="gene")

pnew2<-ggplot(res_n_2genes_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,100))
pnew+pnew2 

wilcox.test(res_n_2genes_degs[padj<0.05]$gene_score,res_n_2genes_degs[padj>=0.05]$gene_score) #0.000276 > 0.0018, donc moins bien 

#meme en affinant le ncpg weight ?
res_n_2genes[,n_cpg.gene:=.N,by=.(gene)]
res_n_2genes[,n_cpg_sig_nom3.gene:=sum(pval<10^-3),by=.(gene)]
res_n_2genes[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.45,by="gene"]
res_n_2genes[,gene_score:=sum(cpg_score)*n_cpg_weight,by="gene"]
ggplot(unique(res_n_2genes,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()
ggplot(unique(res_n_2genes,by=c("gene")))+geom_point(aes(x=n_cpg.gene,y=gene_score))+scale_x_log10()

res_n_2genes_degs<-merge(unique(res_n_2genes[order(gene,pval)],by="gene"),degs[,.(gene,pvalue,padj,log2FoldChange)],by="gene")

pnew3<-ggplot(res_n_2genes_degs)+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score),fill=padj<0.05),outlier.shape = NA)+coord_cartesian(ylim = c(-10,100))
pnew+pnew2+pnew3
wilcox.test(res_n_2genes_degs[padj<0.05]$gene_score,res_n_2genes_degs[padj>=0.05]$gene_score) #0.014 donc oui, moins bien 
#==> garder le closest gene

#2) more stringeant eQTR is best ?
#ac a) influence == tss


#ac b) influence +/- ~ localisation in enh/prom


