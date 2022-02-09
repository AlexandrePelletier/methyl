

#eqtr opti
#I) opti region : top snp +/-1000pb 
#1) liver 
message("creation of regulatory region in liver based on GTEx_Analysis_v8_eQTL signif_variant_gene_pairs")
eqtls_liver<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz")
print(paste(nrow(eqtls_liver),"signif_variant_gene_pairs")) #
print(paste(length(unique(eqtls_liver)$gene_id),"gene")) #
eqtls_liver[,n_eqtl.gene:=.N,by="gene_id"]
print(paste("summary of n snp asso by gene : ",nrow))
print(summary(unique(eqtls_liver,by="gene_id")$n_eqtl.gene))

message("extract chr and pos hg38 of snps")
eqtls_liver[,chr.hg38:=str_extract(variant_id,"^chr[0-9XYM]{1,2}") ]
eqtls_liver<-eqtls_liver[!is.na(chr.hg38),]

eqtls_liver[,pos.hg38:=as.numeric(str_sub(str_extract(variant_id,"_[0-9]+"),2)) ]
eqtls_liver

message("calculate dist snps asso to a same gene")
GetDistToNext<-function(positions){
  return(c(sapply(1:(length(positions)-1), function(i)abs(positions[i]-positions[i+1])),NA))
}

eqtls_liver[,dist_to_next:=GetDistToNext(pos)]
print("summary :")
print(summary(eqtls_liver$dist_to_next))
message("finding min.local")

message("definition of eQTR : 1) search most signif eQTL locally (at +/-10kb)")
eqtls_liver[,minLocal.cand:=pval_nominal<quantile(pval_nominal,0.25),by=c("gene_id")]

is.minLocal<-function(dists,pvals){
  isMins<-sapply(1:length(dists), function(i){
    locisA5kb<-which(dists>(dists[i]-5000)&dists<(dists[i]+5000))
    
    if(all(pvals[locisA5kb]>=pvals[i])){
      return(T)
    }else{
      return(F)
    }
  })
  
  return(isMins)
}

eqtls_liver[minLocal.cand==T,minLocal:=is.minLocal(tss_distance,pval_nominal),by=c("gene_id")]

findeQTR<-function(minLocal,pvals,dists,seed=1500,log10.pval.dt.thr=4,length.max=10000,length.min=1000){
  reg<-0
  length.seed.reg<-1
  eQTRs<-rep(NA,length(pvals))
  iMinsLocals<-which(minLocal==T)
  for(i in iMinsLocals){
    #pour chaque minlocal, si pas deja affectÃ© a une reg : 
    if(!(i%in%which(!is.na(eQTRs)))){
      # => une nouvelle reg, avec un nouveau pval.thr de seed : 
      reg<-reg+1
      
      pvals.thr<-10^(log10(pvals[i])+log10.pval.dt.thr)
      # si eQTL a une dist <seed =>  inclus dans reg
      idansSeedsReg<-which(dists>(dists[i]-seed)&dists<(dists[i]+seed))
      idansReg<-idansSeedsReg
      eQTRs[idansReg]<-reg
      #pour chaque nouvelle reg si length.seed.reg < length.max, regarde si eQTL a une dist <seed, si oui, :
      length.seed.reg<- max(dists[idansSeedsReg])-min(dists[idansSeedsReg])
      
      if(length.seed.reg>0){
        while(length.seed.reg<length.max){
          
          start.seed.reg<-min(dists[idansSeedsReg])
          end.seed.reg<-max(dists[idansSeedsReg])
          candIdansSeedsReg<-setdiff(which(dists>(start.seed.reg-seed)&dists<(end.seed.reg+seed)),idansSeedsReg)
          if(length(candIdansSeedsReg)>0){
            
            #si length.seed.reg <length.min => inclus dans seed.reg
            if(length.seed.reg<length.min){
              idansSeedsReg<-c(idansSeedsReg,candIdansSeedsReg)
              eQTRs[idansSeedsReg]<-reg
              start.seed.reg<-min(dists[idansSeedsReg])
              end.seed.reg<-max(dists[idansSeedsReg])
              length.seed.reg<-end.seed.reg-start.seed.reg
            }
            #si length.seed.reg > length.min => inclus mais stopseed si diff pval min <log10.pval.dt.thr
            else{
              #ajout eQTL dans reg mais pas oblogatoirement dans seeed : 
              idansReg<-c(idansSeedsReg,candIdansSeedsReg)
              
              #inclus dans seed si assez signif :
              sigIdansSeedsReg<-candIdansSeedsReg[which(pvals[candIdansSeedsReg]<pvals.thr)]
              #sil y a des signifs on continue le seed :
              if(length(sigIdansSeedsReg)>0){
                idansSeedsReg<-c(idansSeedsReg,sigIdansSeedsReg)
                eQTRs[idansSeedsReg]<-reg
                start.seed.reg<-min(dists[idansSeedsReg])
                end.seed.reg<-max(dists[idansSeedsReg])
                length.seed.reg<-end.seed.reg-start.seed.reg
              }else{
                #sinon, on inclus les eQTL dans reg et on stop seed
                eQTRs[idansReg]<-reg
                length.seed.reg<-length.max
                
              }
              
            }
            
            
          }else{
            length.seed.reg<-length.max
          }
        }
        #sinon, on inclus les eQTL dans reg et on stop seed
        eQTRs[idansReg]<-reg
        
        
      }
      
    }
    
    
  }
  
  
  
  
  return(eQTRs)
}

message("cluster snps by region")
eqtls_liver[,eQTR:=findeQTR(minLocal,pval_nominal,tss_distance,seed = 2000,log10.pval.dt.thr = 4,length.max=10000,length.min = 1000),by="gene_id"]
message("define eQTR : between snps clustered and at least 1kb around the top snps")
eqtls_liver[!is.na(eQTR),is.start.eQTR:= tss_distance==min(tss_distance),by=.(gene_id,eQTR)]
eqtls_liver[!is.na(eQTR),is.end.eQTR:= tss_distance==max(tss_distance),by=.(gene_id,eQTR)]

eqtls_liver[,start.eQTR.hg38:=min(pos.hg38),by=.(gene_id,eQTR)]
eqtls_liver[,end.eQTR.hg38:=max(pos.hg38),by=.(gene_id,eQTR)]

print("n eQTR by gene :")
eqtls_liver[,n.eQTL:=.N,by=.(gene_id,eQTR)]#[reprendre here]
print(summary(unique(eqtls_liver,by=c("gene_id"))$n.eQTL))

print("n eQTL in eQTR :")
eqtls_liver[,n.eQTL:=.N,by=.(gene_id,eQTR)]
print(summary(unique(eqtls_liver,by=c("gene_id","eQTR"))$n.eQTL))
eqtls_liver[n.eQTL==1,start.eQTR.hg38:=start.eQTR.hg38-500]
eqtls_liver[n.eQTL==1,end.eQTR.hg38:=end.eQTR.hg38+500]

message("translate snp pos in h19 pos")
translator<-fread("ref/eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
                  select = c(1,8))
eqtls_liver<-merge(eqtls_liver,translator,all.x=T,by="variant_id")
eqtls_liver[,chr.hg19:=paste0("chr",str_extract(variant_id_b37,"^[0-9XYM]{1,2}")) ]
eqtls_liver[,pos.hg19:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2)) ]
eqtls_liver
eqtls_liver[,start.eQTR.hg19:=min(pos.hg19),by=.(gene_id,eQTR)]
eqtls_liver[,end.eQTR.hg19:=max(pos.hg19),by=.(gene_id,eQTR)]
eqtls_liver[,n.eQTL:=.N,by=.(gene_id,eQTR)]
eqtls_liver[n.eQTL==1,start.eQTR.hg19:=start.eQTR.hg19-500]
eqtls_liver[n.eQTL==1,end.eQTR.hg19:=end.eQTR.hg19+500]

message("translate gene_id in gene symbol")
trans<-TransEnsembltoSymbol(unique(eqtls_liver)$gene_id)
trans<-trans[,gene_id:=ensembl_gene_id_version][,gene:=hgnc_symbol][,-c("ensembl_gene_id_version","hgnc_symbol")]
merge(eqtls_liver,trans,by="gene_id")

fwrite(eqtls_liver,"ref/eQTL/eQTLs_liver_with_eQTR")


#2) cbps


ggplot(eqtls_liver[qval<0.05])+geom_density(aes(x=tss_distance))
ggplot(eqtls_liver[qval<0.05])+geom_density(aes(x=abs(tss_distance)))+scale_x_log10()
summary(abs(eqtls_liver[qval<0.05]$tss_distance))

eqtls_wb<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz")#22k genes
ggplot(eqtls_wb[qval<0.05])+geom_density(aes(x=tss_distance))
ggplot(eqtls_wb[qval<0.05])+geom_density(aes(x=abs(tss_distance)))+scale_x_log10()
summary(abs(eqtls_wb[qval<0.05]$tss_distance))

eqtls_placenta<-fread("ref/eQTL/placenta_eQTL_fabien.csv")
ggplot(eqtls_placenta[qval<0.05])+geom_density(aes(x=tss_distance))
ggplot(eqtls_placenta[qval<0.05])+geom_density(aes(x=abs(tss_distance)))+scale_x_log10()
summary(abs(eqtls_placenta[qval<0.05]$tss_distance))

#merge
cols<-c("gene_name","pval_nominal","qval","variant_id")
eqtls_placenta<-eqtls_placenta[,gene_name:=hgnc_symbol][,-'hgnc_symbol'][,pval_nominal:=pval]
eqtls<-rbind(list())
#1) add tss_dist
tss_dists<-fread("ref/hg19/hg19_refseq_curated.txt.gz")
tss_dists<-tss_dists[!is.na(tss_pos)&!str_detect(chr,"random|alt|fix|hap|gl")]
eqtls_pos_hg19<-unique(eqtls[,.(chr.hg19,pos.hg19,gene,variant_id)])
eqtls_pos_hg19<-eqtls_pos_hg19[!is.na(gene)&gene!=""]
eqtls_pos_hg19<-merge(eqtls_pos_hg19,tss_dists,all.x=T,allow.cartesian = T)
length(unique(eqtls_pos_hg19[is.na(tss_pos)]$gene))
length(unique(eqtls_pos_hg19[!is.na(tss_pos)]$gene))




#II) more stringeant / sophisticated






#geneScore eqtl +/- et sans apriori sur le gène
source("scripts/utils/new_utils.R")
res_cbps<-fread("analyses/2020-12-08_cbps_new_genescore/res_gene_score_lgaF_vs_ctrlF.csv")
res_cbps[in_eQTR==T,keep:=T]
res_cbps[in_eQTR==F,keep:=abs(tss_dist)==min(abs(tss_dist)),by="cpgID"]
res_cbps<-res_cbps[keep==T]
res_cbps[,tissue:="cbps"]
degs<-fread("../singlecell/analyses/04-DEG_in_LGA/2020-07-08_pseudo_bulk_DEseq2_LgaVsCtrl_CBP1andcbp558_559_samples_excluded_regr_on_batch_and_sex_all_genes.csv")
degs[is.na(padj),padj:=1]


res_cbps<-merge(res_cbps,degs[,.(gene,log2FoldChange,pvalue,padj)])
res_cbps[,tissue:="cbps"]


res_liver<-fread("analyses/liver/res_gene_score_HCCvsCtrl_liver.csv")
res_liver[in_eQTR==T,keep:=T]
res_liver[in_eQTR==F,keep:=abs(tss_dist)==min(abs(tss_dist)),by="cpgID"]
res_liver<-res_liver[keep==T]
res_liver[,tissue:="liver"]
res_liver[,meth.change:=meth.change*100]
degs<-unique(fread("analyses/liver/2021-01-11_res_1gene_by_cpg_bonne_correl_with_degs.csv")[,.(gene,log2FoldChange,pvalue,padj)])

res_liver<-merge(res_liver,degs[,.(gene,log2FoldChange,pvalue,padj)])
res_liver[,tissue:="liver"]


res<-rbind(res_cbps[,.SD,.SDcols=intersect(colnames(res_cbps),colnames(res_liver))],res_liver[,.SD,.SDcols=intersect(colnames(res_cbps),colnames(res_liver))])
res

res[,cpg_score:=(-log10(pval)/4*meth.change)*regul_weight*links_weight]

#geneScore eqtl - :
res_1<-res[in_eQTR==F]
res_1[,method:="closest.cpg"]
res_1[,n_cpg.gene:=.N,by=.(gene,tissue)]
res_1[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.25,by=.(gene,tissue)]
res_1[,gene_score:=sum(cpg_score)*n_cpg_weight,by=.(gene,tissue)]
test<-wilcox.test(res_1[tissue=="liver"&padj<0.01]$gene_score,res_1[tissue=="liver"&padj>=0.01]$gene_score) #0.006575
test
res_1[,p.gs_degs:=wilcox.test(unique(gene_score[padj<0.05]),unique(gene_score[padj>=0.05]))$p.value,by="tissue"]
ggplot(unique(res_1,by=c("gene",'tissue')))+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score)),outlier.shape = NA)+
  facet_wrap("tissue",scales = "free")+scale_y_continuous(limits = c(0,125))

#geneScore eqtl +
res_2<-copy(res)
res_2[,method:="closest.cpg+eQTR_asso_large"]
res_2[,n_cpg.gene:=.N,by=.(gene,tissue)]
res_2[,n_cpg_weight:=(1/sum(1/(abs(cpg_score)+1)))^0.25,by=.(gene,tissue)]
res_2[,gene_score:=sum(cpg_score)*n_cpg_weight,by=.(gene,tissue)]

res_2[,p.gs_degs:=wilcox.test(unique(gene_score[padj<0.05]),unique(gene_score[padj>=0.05]))$p.value,by="tissue"]
ggplot(unique(res_2,by=c("gene",'tissue')))+geom_boxplot(aes(x=padj<0.05,y=abs(gene_score)),outlier.shape = NA)+
  facet_wrap("tissue",scales = "free")+scale_y_continuous(limits = c(0,50))

res_g<-rbind(unique(res_1,by=c("gene",'tissue')),unique(res_2,by=c("gene",'tissue')))
res_g




