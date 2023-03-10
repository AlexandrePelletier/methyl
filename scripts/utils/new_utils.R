
library("data.table")
library("stringr")
library("ggplot2")
library("ggrepel")
library("patchwork")
library("here")
fp<-function(...)file.path(...)

ps<-function(...,sep="",collapse = NULL)paste(...,sep=sep,collapse = collapse)

bed_inter<- function(a, b, opt1="-wa", opt2="-wb",out_dir=".", select=NULL, col.names=NULL){
  require(data.table)
  l<-list(a,b)
  files_to_rm<-c(FALSE,FALSE)
  file_paths<-sapply(1:2, function(i){
    x<-l[[i]]
    if(is.data.frame(x)){
      files_to_rm[i]<-TRUE
      file_path<-fp(out_dir,paste0("temp",i,".bed"))
      fwrite(x,file_path,sep="\t",col.names = FALSE)
      
    }else{
      file_path<-x
    }
    return(file_path)
  })
  
  out_file<-fp(out_dir,"temp_inter.bed")
  cmd<-paste("bedtools intersect -a",file_paths[1],"-b",file_paths[2], opt1, opt2,">",out_file)
  message("run in shell : ",cmd)
  system(cmd)
  message("done.")
  
  if(!is.null(col.names)){
    dt<-fread(out_file,select = select,col.names = col.names)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }else{
    dt<-fread(out_file,select = select)
    file.remove(out_file)
    file.remove(file_paths[files_to_rm])
    return(dt)
  }
}




GetMartGenes<-function()biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
GetMartReg<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")
GetMartMotif<-function()biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_motif_feature")



GetBiomartAttrs<-function(mart)data.table::data.table(biomaRt::listAttributes(mart))
GetBiomartFilter<-function()data.table::data.table(biomaRt::listFilter(mart))

GetBiomartAttrs_reg<-function()data.table::data.table(listAttributes(biomaRt::useEnsembl(biomart = "regulation",dataset = "hsapiens_regulatory_feature")))
GetBiomartFilter_reg<-function()data.table::data.table(biomaRt::listFilter(biomaRt::useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature")))

TransNMtoSymbol<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
      filters = 'refseq_mrna', 
      values = refseq_ids, 
      mart = GetMartGenes())))
  }

TransSymboltoNM<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
      filters = 'hgnc_symbol', 
      values = hgnc_symbols, 
      mart = GetMartGenes())))
}

TransEnsembltoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'ensembl_gene_id', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransSymboltoEnsembl<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'hgnc_symbol', 
      values = hgnc_symbols, 
      mart = GetMartGenes())))
}

TransEnsemblVerstoSymbol<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
      filters = 'ensembl_gene_id_version', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransSymboltoEnsemblVers<-function(hgnc_symbols){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'hgnc_symbol'),
      filters = 'hgnc_symbol', 
      values = hgnc_symbols, 
      mart = GetMartGenes())))
}



TransEnsembltoNM<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna'),
      filters = 'ensembl_gene_id', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransNMtoEnsembl<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id', 'refseq_mrna'),
      filters = 'refseq_mrna', 
      values = refseq_ids, 
      mart = GetMartGenes())))
}

TransEnsemblVerstoNM<-function(ensembl_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna'),
      filters = 'ensembl_gene_id_version', 
      values = ensembl_ids, 
      mart = GetMartGenes())))
  }

TransNMtoEnsemblVers<-function(refseq_ids){
  require(biomaRt)
  require(data.table)
      return(data.table::data.table(biomaRt::getBM(attributes = c('ensembl_gene_id_version', 'refseq_mrna'),
      filters = 'refseq_mrna', 
      values = refseq_ids, 
      mart = GetMartGenes())))
  }



tr<-function(ids_sepBySlash,retourne="all",sep="/",tradEntrezInSymbol=FALSE,uniqu=TRUE){
  IDs<-as.vector(strsplit(ids_sepBySlash,sep)[[1]])
  if(retourne=="all"){
    ret<-1:length(IDs)
  }else{
    ret<-retourne
  }
  if(tradEntrezInSymbol){
    require(clusterProfiler)
    library(org.Hs.eg.db)
    if(retourne=="all"){
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL)
    }
    else{
      return(clusterProfiler::bitr(IDs, fromType = "ENTREZID",toType =  "SYMBOL",OrgDb = org.Hs.eg.db)$SYMBOL[ret])
    }
    
  }else{
    if(uniqu){
      return(unique(IDs[ret]))
    }else{
      return(IDs[ret])
    }
    
    
  }
  
}


GetVarPCs<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}





# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x,return_dt=T){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  if(return_dt)return(data.table(genesV2))

  mousex <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(mousex))
  return(mousex)
}

convertMouseGeneList <- function(x,return_dt=T){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  if(return_dt)return(data.table(genesV2))
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}


#Over repre / GSEA


OR<-function(set1,set2,size_universe){
  if(any(duplicated(set1))){
    set1<-unique(set1)
  }
  if(any(duplicated(set2))){
    set2<-unique(set2)
    }
  phyper(q=sum(set1%in%set2)-1, 
       #number of white balls drawn without replacement from an urn which contains both black and white balls. 
       #<=> number of tcf4 target in DEGs. "-1" because normally give P[X > x] but we want P[X >= x])
       m=length(set2), #the number of white balls in the urn. <=> number of DEGs 
       n=size_universe-length(set2), #the number of black balls in the urn. <=> number of genes tested - number of DEGs
       k=length(set1), #the number of balls drawn from the urn <=> numbers of tcf4 targets 
       lower.tail=FALSE)  # if TRUE (default), probabilities are P[X ??? x] (under-representation), otherwise, P[X > x] (over-representation).

}



OR2<-function(querys,terms_list,size_universe,min.term.size=0,max.term.size=Inf,overlap_column=TRUE,verbose=FALSE){
  if(is.list(querys)){
    return(Reduce(rbind,lapply(names(querys),
                             function(q)OR2(querys = querys[[q]],
                                              terms_list = terms_list,
                                              size_universe = size_universe,
                                            min.term.size = min.term.size,
                                            max.term.size = max.term.size,
                                            overlap_column = overlap_column,
                                            verbose = verbose )[,query:=q])))
  }else{
  res_or<-data.table(term=names(terms_list),term.size=sapply(terms_list,length))
  res_or<-res_or[term.size<=max.term.size]
  n_terms<-nrow(res_or)
  if(verbose)message(length(terms_list)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
  res_or<-res_or[term.size>=min.term.size]
  if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
  
  res_or[,n.query:=length(querys)]
  res_or[,n.overlap:=sum(querys%in%terms_list[[term]]),by="term"]
  if(overlap_column==TRUE){
    res_or[,genes.overlap:=paste(querys[querys%in%terms_list[[term]]],collapse="|"),by="term"]
  }
  res_or[,pct.query.overlap:=n.overlap/n.query]
  res_or[,precision:=pct.query.overlap]

  res_or[,pct.term.overlap:=n.overlap/term.size]
  
  res_or[,background_size:=size_universe]

  res_or[,pct.background:=term.size/size_universe] #TO IMPROVE (Here size universe can be the intersection between 2 "univers" while n.query is the n.query is the query universe. 
  

  res_or[,pval:=phyper(q=n.overlap-1, 
                     m=n.query, 
                     n=size_universe-n.query, 
                     k=term.size, 
                     lower.tail=FALSE),
       by="term"]
  res_or[,padj:=p.adjust(pval,method = 'BH')]
  if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
  return(res_or)
      }
  
}

OR3<-function(querys,terms_list,background,min.term.size=0,max.term.size=Inf,overlap_column=TRUE,verbose=FALSE){
  if(is.list(querys)){
    dt<-Reduce(rbind,lapply(names(querys),
                             function(q)OR3(querys = querys[[q]],
                                              terms_list = terms_list,
                                              background = background,
                                            min.term.size = min.term.size,
                                            max.term.size = max.term.size,
                                            overlap_column = overlap_column,
                                            verbose = verbose )[,query:=q]))
    
    return(dt[,query.:=query][,.SD,.SDcols=c(ncol(dt),1:(ncol(dt)-1))])
  }else{
  queryf<-intersect(querys,background)
  terms_listf<-lapply(terms_list,function(x)intersect(x,background))
  res_or<-data.table(term=names(terms_listf),term.size=sapply(terms_listf,length))
  res_or<-res_or[term.size<=max.term.size]
  n_terms<-nrow(res_or)
  if(verbose)message(length(terms_listf)-n_terms, " terms were filtered due to term.size above the limit of ",max.term.size," genes")
  res_or<-res_or[term.size>=min.term.size]
  n_terms<-nrow(res_or)
  if(verbose)message(n_terms-nrow(res_or), " terms were filtered due to term.size below the limit of ",min.term.size," genes")
  
  res_or[,n.query:=length(queryf)]
  res_or[,n.overlap:=sum(queryf%in%terms_listf[[term]]),by="term"]
  
  res_or[,pct.query.overlap:=n.overlap/n.query]
  res_or[,precision:=pct.query.overlap]

  res_or[,pct.term.overlap:=n.overlap/term.size]
  
  res_or[,background_size:=length(background)]

  res_or[,pct.term.background:=term.size/background_size] 
  

  res_or[,pval:=phyper(q=n.overlap-1, 
                     m=term.size, 
                     n=background_size-term.size, 
                     k=n.query, 
                     lower.tail=FALSE),
       by="term"]
  res_or[,padj:=p.adjust(pval,method = 'BH')]
  res_or[,fold.enrichment:=pct.query.overlap/pct.term.background]
  if(overlap_column==TRUE){
    res_or[,genes.overlap:=paste(queryf[queryf%in%terms_listf[[term]]],collapse="|"),by="term"]
  }
  if(verbose)message(nrow(res_or[padj<0.05])," terms enriched in your genes of interest with padj<0.05")
  return(res_or)
      }
  
}

#trans in hg38
hg19to38<-function(x){
  in_file<-"outputs/temp_hg19.bed"
  out_file<-"outputs/temp_hg38.bed"
  
  fwrite(x,in_file,col.names = F,sep="\t")
  
  system(paste("CrossMap.py bed ref/hg19ToHg38.over.chain.gz",in_file,out_file))
  trans<-fread(out_file,select=c(1,2,3,4),col.names = c("chr","start","end","id"))
  file.remove(c(in_file,out_file))
  return(trans)
  }


FindGOGenes<-function(terms_or_ids){
  require("biomaRt")
  require("stringr")

  if(!all(str_detect(terms_or_ids,"^GO:")))terms_or_ids=FindGO_ID(term_description=terms_or_ids)
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
  #gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
  gene.data <- getBM(attributes=c('hgnc_symbol','go_id'),
                   filters = 'go', values = terms_or_ids, mart = ensembl,uniqueRows = T)
  return(data.table(gene.data)[go_id%in%terms_or_ids])
    }



FindGO_ID<-function(term_descriptions){
  require("GO.db")
  terms<-Term(GOTERM)
  ids<-names(terms)[match(term_descriptions,unlist(terms))]
  return(ids)
}

####Signac####
GetMotifIDs<-function(object,motif.names,assay=NULL,return_dt=FALSE){
  if(is.null(assay))assay<-DefaultAssay(object)
  idx<-match(motif.names,object@assays[[assay]]@motifs@motif.names)
  if(return_dt){
    return(
      data.table(motif.name=motif.names,
                 motif.id=names(object@assays[[assay]]@motifs@motif.names[idx]))
      )
    }else{
  return(names(object@assays[[assay]]@motifs@motif.names[idx]))
    }
  
}

CheckMotif<-function(object,peaks,motif.name,assay = NULL,return.peaks=FALSE){
  require("Signac")
  if(is.null(assay))assay<-DefaultAssay(object)
  motif<-GetMotifIDs(object,motif.name,assay=assay)
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  
  motifs_peaks_tf <- motif.all[peaks,motif , drop = FALSE]
  if(return.peaks){
    motifs_peaks_tf<-rownames(motifs_peaks_tf)[as.vector(motifs_peaks_tf==1)]
    return(motifs_peaks_tf)
  }else{
    motifs_peaks_tf_vec<-as.vector(motifs_peaks_tf==1)
    names(motifs_peaks_tf_vec)<-rownames(motifs_peaks_tf)
    return(motifs_peaks_tf_vec)
  }
  
 
}
GetMotifIDs<-function(object,motif.names,assay=NULL,return_dt=FALSE){
  if(is.null(assay))assay<-DefaultAssay(object)
  idx<-match(motif.names,object@assays[[assay]]@motifs@motif.names)
  if(return_dt){
    return(
      data.table(motif.name=motif.names,
                 motif.id=names(object@assays[[assay]]@motifs@motif.names[idx]))
      )
    }else{
  return(names(object@assays[[assay]]@motifs@motif.names[idx]))
    }
  
  }
GetMotif<-function(object,peaks,motifs=NULL,assay = "peaks"){
#return tfmotif peak data.table
  require("Signac")
  require("data.table")

  if(is.null(assay))assay<-DefaultAssay(object)
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  if(is.null(motifs)) motifs<-colnames(motif.all)
  motif.filtered<-motif.all[peaks,motifs,drop=F]
  motif_dt<-melt(data.table(as.matrix(motif.filtered),keep.rownames = "peak"),
                 variable.name="motif",value.name = "presence")
  motif_dt[,presence:=as.logical(presence)]
  motif_dt[(presence)] 
  motif_dt<-motif_dt[(presence)][,-"presence"]
  
  motifsnames<-data.table(motif.name=object@assays$peaks@motifs@motif.names,motif=names(object@assays$peaks@motifs@motif.names))
  motif_dt<-merge(motif_dt,motifsnames)
  return(motif_dt)
}

  
 

start<-function(x)sapply(x,function(x)as.numeric(strsplit(x,"-")[[1]][2]))
end<-function(x)sapply(x,function(x)as.numeric(strsplit(x,"-")[[1]][3]))
seqid<-function(x)sapply(x,function(x)strsplit(x,"-")[[1]][1])
  

MethChangeReg<-function(res_meth,region){
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
  res_meth_reg[,start:=pos][,end:=pos+1]
  return(res_meth_reg)
  }
MethChangePlot<-function(res_meth,region,limits=NULL,breaks=waiver()){
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  res_meth_reg<-res_meth[chr==chromosome&pos>start.pos&pos<end.pos]
  res_meth_reg[,start:=pos][,end:=pos+1]
  p<-ggplot(data = res_meth_reg) + geom_segment(aes(x = start, y = 0, 
        xend = end, yend = logFC,col=-log10(P.Value)), size = 2, data = res_meth_reg)+
    scale_color_gradient(low = "white",high = "black",limits=limits,breaks=breaks)
  
  p<-p+ theme_classic() + ylab(label = "Methylation change") + 
      xlab(label = paste0(chromosome, " position (bp)")) + 
      xlim(c(start.pos, end.pos))

  return(p)
}

TFMotifPlot<-function(object,region,motif.name,assay=NULL){
  if(is.null(assay))assay<-DefaultAssay(object)
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)
  ranges<-object@assays[[assay]]@motifs@positions[[GetMotifIDs(object,motif.names = motif.name)]]
  dt <- data.table(as.data.frame(ranges))
  dt_reg<-dt[seqnames==chromosome&start>start.pos&end<end.pos]
  
  p<-ggplot(data = dt_reg) + geom_segment(aes(x = start, y = 0, 
        xend = end, yend = 0),col="black", size = 2, data = dt_reg)
  
  p<-p+ theme_classic() + ylab(label = motif.name) + 
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
      xlab(label = paste0(chromosome, " position (bp)")) + 
      xlim(c(start.pos, end.pos))
  return(p)
}

TFsMotifPlot<-function(object,region,motif.names,assay=NULL,size=2,alpha=1,pad=0){
  if(is.null(assay))assay<-DefaultAssay(object)
  start.pos <- start(region)
  end.pos <- end(region)
  chromosome <- seqid(region)

  dt_region <-Reduce(rbind,lapply(motif.names,function(x){
    
    ranges<-object@assays[[assay]]@motifs@positions[GetMotifIDs(object,motif.names = x)]

    dt<-data.table(as.data.frame(ranges))[seqnames==chromosome&start>start.pos&end<end.pos][,motif.name:=x]
    return(dt)
    }
    ))
  
  
  p<-ggplot(data = dt_region) + geom_segment(aes(x = start-pad, y = 0, 
        xend = end+pad, yend = 0,col=motif.name), size = size,alpha=alpha, data = dt_region)
  
  p<-p+ theme_classic() + ylab(label = "TF motif") + 
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
      xlab(label = paste0(chromosome, " position (bp)")) + 
      xlim(c(start.pos, end.pos))
  return(p)
}
 

pctPC<-function(pca,rngPCs="all"){
  if(is.character(rngPCs)){
    rngPCs<-1:length(pca$sdev)
  }
  pct.varPCs<-pca$sdev[rngPCs]^2/sum(pca$sdev^2)
  names(pct.varPCs)<-rngPCs
  return( pct.varPCs)
}
