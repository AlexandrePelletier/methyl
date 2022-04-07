#02A1-create_eQTR

source("scripts/utils/new_utils.R")
out<-"outputs/02A1-create_eQTR"
dir.create(out)
#I) with whole_blood eQTL
eqtls_wb<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
eqtls_wb #2 414 653 signif_variant_gene_pairs
unique(eqtls_wb$gene_id) #1000 genes
eqtls_wb[,n_eqtl.gene:=.N,by="gene_id"]

summary(unique(eqtls_wb,by="gene_id")$n_eqtl.gene)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   #  1.0    22.0    87.0   195.4   230.0  8571.0

#extract chr and pos hg38 of snps
eqtls_wb[,chr.hg38:=str_extract(variant_id,"^chr[0-9XYM]{1,2}") ]
eqtls_wb<-eqtls_wb[!is.na(chr.hg38),]
eqtls_wb[,pos.hg38:=as.numeric(str_sub(str_extract(variant_id,"_[0-9]+"),2)) ]
eqtls_wb

#trans in hg19 :
translator<-fread("ref/eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz",
                  select = c(1,8))
translator


eqtls_wb<-merge(eqtls_wb,translator,all.x=T,by="variant_id")
rm(translator)

eqtls_wb[,chr.hg19:=paste0("chr",str_extract(variant_id_b37,"^[0-9XY]{1,2}"))]

eqtls_wb[,pos.hg19:=as.numeric(str_sub(str_extract(variant_id_b37,"_[0-9]+"),2))]
eqtls_wb

eqtls_wb[,chr:=chr.hg19]
eqtls_wb[,pos:=pos.hg19]

#calculate dist snps asso to a same gene
GetDistToNext<-function(positions){
  return(c(sapply(1:(length(positions)-1), function(i)abs(positions[i]-positions[i+1])),NA))
}

eqtls_wb[,dist_to_next:=GetDistToNext(pos),by="gene_id"]
summary(eqtls_wb$dist_to_next)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's
#    0       114       397     12418      1158 180679628     39827

#definition of eQTR : 1 to 10kb region around most signif eQTL locally

# 1) find min.local, i.e eqtl with pval in top25% of eQTL of the gene, and first eQTL at +/-5kb

eqtls_wb[,minLocal.cand:=pval_nominal<=quantile(pval_nominal,0.25),by=c("gene_id")]

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

eqtls_wb[minLocal.cand==T,minLocal:=is.minLocal(pos,pval_nominal),by=c("gene_id")]
eqtls_wb[minLocal==T] #189470 min local

eqtls_wb[,n.local.gene:=sum(minLocal,na.rm = T),by="gene_id"]
summary(unique(eqtls_wb,by="gene_id")$n.local.gene)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 1.00    3.00    8.00   15.33   18.00  966.00


# 2) find eQTR, i.e regions at +/-500pb around min local iteractively expanded if eQTLs at +/-2000pb and log10(pval) difference < 4
findeQTR<-function(minLocal,pvals,dists,seed=2000,log10.pval.dt.thr=4,length.max=10000,length.min=1000){
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

#cluster snps by region
eqtls_wb[,eQTR:=findeQTR(minLocal,pval_nominal,pos,seed = 2000,log10.pval.dt.thr = 4,length.max=10000,length.min = 1000),by="gene_id"]
eqtls_wb[!is.na(eQTR)] #1.2M / 2.4M eQTL on eQTR


eqtls_wb[!is.na(eQTR),start.eQTR:=min(pos),by=.(gene_id,eQTR)]
eqtls_wb[!is.na(eQTR),end.eQTR:=max(pos),by=.(gene_id,eQTR)]

#n eQTL in eQTR
eqtls_wb[!is.na(eQTR),n.eQTL:=.N,by=.(gene_id,eQTR)]
summary(unique(eqtls_wb[!is.na(eQTR)],by=c("gene_id","eQTR"))$n.eQTL)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 #   1.00    2.00    6.00   10.73   13.00  689.00

#for eQTR with only 1 eQTL, increase regions at +/-500pb from the eQTL
eqtls_wb[n.eQTL==1,start.eQTR:=start.eQTR-500]
eqtls_wb[n.eQTL==1,end.eQTR:=end.eQTR+500]



#translate gene_id in gene symbol
trans<-fread("ref/eQTL/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz",select = 1:2,col.names = c("gene_id","gene"))
trans[gene=="",gene:=NA]
eqtls_wb<-merge(eqtls_wb,trans,by="gene_id")
unique(eqtls_wb,by="gene_id")[is.na(gene)]#0/3179


#summary :
#eqtr by gene
eqtls_wb[,n.eqtr.gene:=length(unique(na.omit(eQTR))),by="gene_id"]

summary(unique(eqtls_wb,by='gene_id')$n.eqtr.gene)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.000   2.000   6.000   9.718  12.000 191.000

eqtls_wb[n.eqtr.gene==191] # ZFP57 Transcription regulator required to maintain maternal and paternal gene imprinting

eqtls_wb[n.eqtr.gene==0] #69 eQTL so negligeable

#save useful info before reduce df to eQTR level
eqtls_wb[minLocal==T,pos.min.local:=pos]
eqtls_wb[,tss_dist:=tss_distance]
eqtls_wb[,pval_link:=pval_nominal]
eqtls_wb[minLocal==T,avg.mlog10pval:=mean(-log10(pval_link)),by=c("gene_id","eQTR")]
eqtls_wb[!is.na(gene),eqtr_id:=paste(gene,eQTR,sep="-")]
summary(unique(eqtls_wb[minLocal==T],by=c("eqtr_id"))$avg.mlog10pval)
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 3.255   6.810  11.308  18.710  21.265 277.112

eqtls_wb[,length.eqtr:=end.eQTR-start.eQTR]
summary(unique(eqtls_wb[minLocal==T],by=c("eqtr_id"))$length.eqtr)
 # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
 #      0    1000    3172    4045    6146   13955     416

fwrite(eqtls_wb,fp(out,"eqtls_wb_with_eQTR.csv.gz"))
eqtls_wb<-fread(fp(out,"eqtls_wb_with_eQTR.csv.gz"))

#and save eQTR in a bed files :
eqtrs_symbol<-unique(eqtls_wb[!is.na(eQTR)&!is.na(gene)&minLocal==T][order(gene,eQTR,abs(tss_dist))][,.(chr,start.eQTR,end.eQTR,eqtr_id,pos.min.local,pval_link,avg.mlog10pval,tss_dist,gene)])
eqtrs_symbol#183k eqtr
unique(eqtrs_symbol,by="gene") #12348 genes
fwrite(eqtrs_symbol,fp(out,"whole_blood_eQTR_symbol.bed"),sep = "\t")



#II) With meta eqtls
#see 02A1a-meta_eQTL
eqtls_meta<-fread("outputs/02A1a-meta_eQTL/tissue_wide_signif_variant_gene_pairs_hg19.csv.gz")
eqtls_meta #720574 signif_variant_gene_pairs
length(unique(eqtls_meta$gene)) #5254 genes
eqtls_meta[,n_eqtl.gene:=.N,by="gene"]

summary(unique(eqtls_meta,by="gene")$n_eqtl.gene)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  #   1.0     5.0    21.0   137.1    87.0 24115.0

#calculate dist snps asso to a same gene
eqtls_meta[,dist_to_next:=GetDistToNext(pos),by="gene"]
summary(eqtls_meta$dist_to_next)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#       0       0       0    3166     722 1868353    5254

#definition of eQTR : 1 to 10kb region around most signif eQTL locally

# 1) find min.local, i.e eqtl with pval in top25% of eQTL of the gene, and first eQTL at +/-5kb

eqtls_meta[,minLocal.cand:=PVALUE_FE<=quantile(PVALUE_FE,0.25),by=c("gene")]

eqtls_meta[minLocal.cand==T,minLocal:=is.minLocal(pos,PVALUE_FE),by=c("gene")]
eqtls_meta[minLocal==T] #98k min local

eqtls_meta[,n.local.gene:=sum(minLocal,na.rm = T),by="gene"]
summary(unique(eqtls_meta,by="gene")$n.local.gene)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 1.00    2.00    5.00   18.71   15.00 1683.00

#cluster snps by region
eqtls_meta[,eQTR:=findeQTR(minLocal,PVALUE_FE,pos,seed = 2000,log10.pval.dt.thr = 4,length.max=10000,length.min = 1000),by="gene"]
eqtls_meta[!is.na(eQTR)] #320k / 720k eQTL on eQTR

eqtls_meta[!is.na(eQTR),start.eQTR:=min(pos),by=.(gene,eQTR)]
eqtls_meta[!is.na(eQTR),end.eQTR:=max(pos),by=.(gene,eQTR)]

#n eQTL in eQTR
eqtls_meta[!is.na(eQTR),n.eQTL:=.N,by=.(gene,eQTR)]
summary(unique(eqtls_meta[!is.na(eQTR)],by=c("gene","eQTR"))$n.eQTL)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 1.00    1.00    3.00    9.02    8.00 1144.00

#for eQTR with only 1 eQTL, increase regions at +/-500pb from the eQTL
eqtls_meta[n.eQTL==1,start.eQTR:=start.eQTR-500]
eqtls_meta[n.eQTL==1,end.eQTR:=end.eQTR+500]

#summary :
#eqtr by gene
eqtls_meta[,n.eqtr.gene:=length(unique(na.omit(eQTR))),by="gene"]

summary(unique(eqtls_meta,by='gene')$n.eqtr.gene)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   # 1.00    1.00    3.00    7.08    7.00  166.00

eqtls_meta[n.eqtr.gene==166] # LY6G5B, gene du CMH III, https://www.genecards.org/cgi-bin/carddisp.pl?gene=LY6G5B


#eQTR_id for eQTR link to gene symbol
eqtls_meta[minLocal==T,pos.min.local:=pos]
eqtls_meta[,pval_link:=PVALUE_FE]

eqtls_meta[minLocal==T,avg.mlog10pval:=mean(-log10(pval_link)),by=c("gene","eQTR")]
eqtls_meta[!is.na(gene),eqtr_id:=paste(gene,eQTR,sep="-")]
summary(eqtls_meta[minLocal==T]$avg.mlog10pval)


eqtls_meta[,length.eqtr:=end.eQTR-start.eQTR]
summary(unique(eqtls_meta[minLocal==T])$length.eqtr)

fwrite(eqtls_meta,fp(out,"eqtls_meta_with_eQTR.csv.gz"))

#and save in a bed files :

eqtrs_symbol<-unique(eqtls_meta[!is.na(eQTR)&!is.na(gene)&minLocal==T][order(gene,eQTR,abs(tss_dist))][,.(chr,start.eQTR,end.eQTR,eqtr_id,pos.min.local,pval_link,avg.mlog10pval,tss_dist,gene)])
eqtrs_symbol#46k eqtr
unique(eqtrs_symbol,by="gene") #5254 genes
fwrite(eqtrs_symbol,fp(out,"tissue_wide_eQTR_symbol.bed"),sep = "\t")


#III) merge eQTRs
eqtrs<-rbind(fread(fp(out,"tissue_wide_eQTR_symbol.bed"))[,tissue:="tissue_wide"],
      fread(fp(out,"whole_blood_eQTR_symbol.bed"))[,tissue:="whole_blood"])
fwrite(eqtrs,fp(out,"all_eQTR_symbol.bed"))
