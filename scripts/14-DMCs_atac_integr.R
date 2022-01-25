#hypermeth on OCR regulating EGR1 network ?
out<-'outputs/14-DMCs_atac_integr'
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
#renv::install("Signac")
library(Signac)

#I)DMCs overlap in OCR ?
#1) overall OCR
atacs<-readRDS("../atac/outputs/cbps_merged/cbps_atac1-4_merged_qc.rds")
ocrs<-data.table(peaks=rownames(atacs))
ocrs[,chr:=str_extract(peaks,"chr[0-9XY]+"),]
ocrs[,start:=strsplit(peaks,"-")[[1]][2],by="peaks"]
ocrs[,end:=strsplit(peaks,"-")[[1]][3],by="peaks"]

res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
cpgs<-fread("ref/all_HELP_tagging_cpgs_hg19_pos.csv")
cpgs[is.na(chr)]
#need trans in hg38 pos
cpgs38<-hg19to38(cpgs[,.(chr,pos,pos,cpg_id)][order(chr,pos)])
cpgs38[,cpg_id:=id]
cpgs38[,pos:=start]
res_cpgs<-merge(res_meth[,.(cpg_id,logFC,AveExpr,P.Value,adj.P.Val)],unique(cpgs38[,.(cpg_id,chr,pos)]))
res_cpgs[is.na(pos)] #ok
fwrite(res_cpgs,fp(out,"res_cpgs_hg38.cs.gz"))
cpgs_in_ocrs<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          ocrs[,.(chr,start,end,peaks)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peaks","start","end"))

length(unique(cpgs_in_ocrs$cpg_id))#300/750k
fwrite(cpgs_in_ocrs,fp(out,"cpgs_in_bulk_OCRs.csv.gz"))
cpgs_in_ocrs<-fread(fp(out,"cpgs_in_bulk_OCRs.csv.gz"))

res_cpgs_ocrs<-merge(res_meth[,.(cpg_id,logFC,AveExpr,P.Value,adj.P.Val)],cpgs_in_ocrs)

res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25] #3910
fwrite(res_cpgs_ocrs,fp(out,"res_meth_in_bulk_OCRs.csv.gz"))

fwrite(res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25],fp(out,"DMCs_in_bulk_OCRs.csv.gz"))

res_meth[P.Value<0.001&abs(logFC)>25] #3910 / 4815 DMCs
res_cpgs_ocrs[,peak_size:=end-start]
ggplot(res_cpgs_ocrs)+geom_density(aes(peak_size))

#TF motif enrichment in this peaks
renv::install("bioc::JASPAR2020")
renv::install("bioc::TFBSTools")

library(JASPAR2020)
library(TFBSTools)

renv::install("bioc::BSgenome.Hsapiens.UCSC.hg38")
renv::install("bioc::motifmatchr")

library(BSgenome.Hsapiens.UCSC.hg38)

# Get a list of motif position frequency matrices from the JASPAR database
#?getMatrixSet
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
#?AddMotifs
atacs <- AddMotifs(
  object = atacs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

peaks_with_dmcs<-unique(res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25]$peaks)
peaks_without<-unique(res_cpgs_ocrs[!peaks%in%peaks_with_dmcs]$peaks)
enriched.motifs <- FindMotifs(
  object = atacs,
  features = peaks_with_dmcs ,
  background=peaks_without)


enriched.motifs<-data.table(enriched.motifs)
fwrite(enriched.motifs,fp(out,"enriched_motifs_in_CpGs_overlaping_OCRs_with_DMCs_vs_without.csv"))

enriched.motifs[pvalue<0.000001]$motif.name
enriched.motifs[,motif.name:=factor(motif.name,levels=enriched.motifs[order(pvalue)]$motif.name)]

ggplot(enriched.motifs[order(pvalue)][1:24])+
  geom_point(aes(x=motif.name,size=fold.enrichment,y=observed,col=-log10(pvalue+10^-316)))

MotifPlot(
  object = atacs,
  motifs = head(enriched.motifs$motif)
)

# peaks with EGR1/KLF2/KLF4/.. motif close to which genes? genes of the EGRN ?  genes Downregulé ? Genes of the EGRN downregules ?
#EGR1
dmcs_ocrs<-fread(fp(out,"DMCs_in_bulk_OCRs.csv.gz"))
enriched.motifs<-fread(fp(out,"enriched_motifs_in_CpGs_overlaping_OCRs_with_DMCs_vs_without.csv"))


peaks_dmcs<-unique(dmcs_ocrs$peaks)

motif.all <- GetMotifData(
    object = atacs, assay = "ATAC", slot = "data"
  )
motifs_dmcs <- motif.all[peaks_dmcs, , drop = FALSE]

egr1_motif_dmcs<-motifs_dmcs[,enriched.motifs[motif.name=="EGR1"]$motif,drop=F]
peaks_dmcs_egr1<-rownames(egr1_motif_dmcs)[as.vector(egr1_motif_dmcs==1)]
length(peaks_dmcs_egr1) #2363
df_genes_dmcs_egr1<-ClosestFeature(atacs,regions = peaks_dmcs_egr1)
genes_dmcs_egr1<-unique(df_genes_dmcs_egr1$gene_name)
length(genes_dmcs_egr1) #2222

#genes of the EGRN ?  
tftargets<-fread("outputs/13-GRN_integr/tf_target_interactions.csv")
EGRn<-fread("outputs/13-GRN_integr/egr1_network_tf_target_interactions.csv")
egrn_genes<-unique(EGRn$target)
length(egrn_genes)#208
intersect(genes_dmcs_egr1,egrn_genes) 
# [1] "ID3"      "TINAGL1"  "RLF"      "ATP1B1"   "IER5"     "ATF3"     "PTGER4"  
#  [8] "TBCC"     "IFRD1"    "KLF6"     "ZEB1"     "DDIT4"    "TRIM8"    "CD151"   
# [15] "CDKN1C"   "AHNAK"    "ETS1"     "C12orf57" "UBC"      "TSC22D1"  "INTS6"   
# [22] "NFKBIA"   "ZFP36L1"  "MOAP1"    "EIF5"     "KLF13"    "HEXIM1"   "TOB1"    
# [29] "SOCS3"    "TGIF1"    "PRNP"  
#32/208
OR(genes_dmcs_egr1,egrn_genes,size_universe = length(unique(tftargets$target))) #p=1

#genes downregules ?
res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")
dn_genes<-res_degs[lineage=="HSC"&padj<0.05&log2FoldChange<(-0.5)]$gene
length(dn_genes)#285
intersect(genes_dmcs_egr1,dn_genes) 
# [1] "IFFO2"     "LUZP1"     "C1orf109"  "PLK3"      "BEND5"     "EFNA3"    
#  [7] "SEMA4A"    "ATP1B1"    "IER5"      "SLC30A1"   "TAF5L"     "FAM43A"   
# [13] "RPP38"     "BAMBI"     "UBE2D1"    "PPRC1"     "TAF5"      "DUSP5"    
# [19] "BAG3"      "PHLDA2"    "C11orf96"  "SLC3A2"    "STIP1"     "H2AFX"    
# [25] "C2CD2L"    "RAB3IP"    "IKBIP"     "CRY1"      "EFNB2"     "NFKBIA"   
# [31] "BRMS1L"    "GCH1"      "JDP2"      "CYP1A1"    "TNFRSF12A" "CFAP20"   
# [37] "ATP6V0D1"  "SPATA2L"   "DERL2"     "PSMC3IP"   "PLEKHH3"   "HEXIM1"   
# [43] "SOCS3"     "TOP1"       
#43/285 DEGs
OR(genes_dmcs_egr1,dn_genes,size_universe = length(unique(res_degs$gene))) #p=0.81

#Genes of the EGRN downregules ?
genes_dn_egr1<-intersect(egrn_genes,dn_genes)
length(genes_dn_egr1) #32
intersect(genes_dmcs_egr1,genes_dn_egr1) 
#"ATP1B1" "IER5"   "NFKBIA" "HEXIM1" "SOCS3" 5/32

#32/208
OR(genes_dmcs_egr1,genes_dn_egr1,size_universe = length(unique(intersect(tftargets$target,res_degs$gene)))) #p=1


#KLF2
enriched.motifs[order(pvalue)]$motif.name
query_motif_dmcs<-motifs_dmcs[,enriched.motifs[motif.name=="KLF2"]$motif,drop=F]
peaks_dmcs_query<-rownames(query_motif_dmcs)[as.vector(query_motif_dmcs==1)]
length(peaks_dmcs_query) #2445
df_genes_dmcs_query<-ClosestFeature(atacs,regions = peaks_dmcs_query)
genes_dmcs_klf2<-unique(df_genes_dmcs_query$gene_name)
length(genes_dmcs_klf2) #2271
#genes of the EGRN ?  
intersect(genes_dmcs_klf2,egrn_genes)
#  [1] "ID3"      "TINAGL1"  "RLF"      "ATP1B1"   "ATF3"     "PTGER4"   "TBCC"    
#  [8] "IFRD1"    "KLF6"     "ZEB1"     "DDIT4"    "TRIM8"    "CD151"    "CDKN1C"  
# [15] "AHNAK"    "ETS1"     "C12orf57" "UBC"      "TSC22D1"  "INTS6"    "NFKBIA"  
# [22] "ZFP36L1"  "MOAP1"    "EIF5"     "KLF13"    "HEXIM1"   "TOB1"     "SOCS3"   
# [29] "TGIF1"    "PRNP"    
OR(genes_dmcs_klf2,egrn_genes,size_universe = length(unique(tftargets$target))) #p=1

#genes DEGs ?
intersect(genes_dmcs_klf2,dn_genes) 
#  [1] "IFFO2"     "LUZP1"     "PLK3"      "BEND5"     "GSTM3"     "TUFT1"    
#  [7] "EFNA3"     "ATP1B1"    "SLC30A1"   "TAF5L"     "FAM43A"    "RPP38"    
# [13] "BAMBI"     "UBE2D1"    "TAF5"      "DUSP5"     "BAG3"      "PHLDA2"   
# [19] "SLC3A2"    "STIP1"     "CDK2AP2"   "C2CD2L"    "RAB3IP"    "IKBIP"    
# [25] "CRY1"      "EFNB2"     "NFKBIA"    "BRMS1L"    "GCH1"      "JDP2"     
# [31] "TRIP11"    "CYP1A1"    "KCTD5"     "TNFRSF12A" "CFAP20"    "ATP6V0D1" 
# [37] "SPATA2L"   "DERL2"     "PSMC3IP"   "PLEKHH3"   "HEXIM1"    "METTL2A"  
# [43] "SOCS3"  
OR(genes_dmcs_klf2,dn_genes,size_universe = length(unique(res_degs$gene))) #p=0.89

#KLF4
enriched.motifs[order(pvalue)]$motif.name
query_motif_dmcs<-motifs_dmcs[,enriched.motifs[motif.name=="KLF4"]$motif,drop=F]
peaks_dmcs_query<-rownames(query_motif_dmcs)[as.vector(query_motif_dmcs==1)]
length(peaks_dmcs_query) #2713
df_genes_dmcs_query<-ClosestFeature(atacs,regions = peaks_dmcs_query)
genes_dmcs_klf4<-unique(df_genes_dmcs_query$gene_name)
length(genes_dmcs_klf4) #2485
#genes of the EGRN ?  
intersect(genes_dmcs_klf4,egrn_genes) 
#  [1] "ID3"      "TINAGL1"  "RLF"      "ATP1B1"   "ATF3"     "PTGER4"   "TBCC"    
#  [8] "PHACTR2"  "IFRD1"    "KLF6"     "ZEB1"     "DDIT4"    "TRIM8"    "CD151"   
# [15] "CDKN1C"   "AHNAK"    "ETS1"     "C12orf57" "UBC"      "TSC22D1"  "INTS6"   
# [22] "NFKBIA"   "ZFP36L1"  "MOAP1"    "EIF5"     "KLF13"    "HEXIM1"   "TOB1"    
# [29] "H3F3B"    "SOCS3"    "TGIF1"    "SERTAD1"  "PRNP" 
OR(genes_dmcs_klf4,egrn_genes,size_universe = length(unique(tftargets$target))) #p=1

#genes DEGs ?
intersect(genes_dmcs_klf4,dn_genes) 
#  [1] "IFFO2"     "LUZP1"     "C1orf109"  "PLK3"      "BEND5"     "GSTM3"    
#  [7] "TUFT1"     "EFNA3"     "SEMA4A"    "ATP1B1"    "SLC30A1"   "TAF5L"    
# [13] "FAM43A"    "RPP38"     "BAMBI"     "UBE2D1"    "PPRC1"     "TAF5"     
# [19] "DUSP5"     "BAG3"      "PHLDA2"    "C11orf96"  "SLC3A2"    "STIP1"    
# [25] "CDK2AP2"   "C2CD2L"    "RAB3IP"    "IKBIP"     "CRY1"      "EFNB2"    
# [31] "NFKBIA"    "BRMS1L"    "GCH1"      "JDP2"      "CYP1A1"    "KCTD5"    
# [37] "TNFRSF12A" "CFAP20"    "ATP6V0D1"  "SPATA2L"   "DERL2"     "PLEKHH3"  
# [43] "HEXIM1"    "SOCS3" 
OR(genes_dmcs_klf4,dn_genes,size_universe = length(unique(res_degs$gene))) #p=0.96



#2) OCR by lineage

#renv::install("harmony")
library(harmony)
DefaultAssay(atacs)<-"ATAC"
atacs <- RunHarmony(
  object = atacs,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)
# re-compute the UMAP using corrected LSI embeddings

atacs <- RunUMAP(atacs, dims = 2:30, reduction = 'harmony',reduction.name = "humap")
DimPlot(atacs, group.by = 'dataset',reduction = "humap", pt.size = 0.1)

atacs <- FindNeighbors(object = atacs, reduction = 'harmony', dims = 2:30)
atacs <- FindClusters(object = atacs, verbose = FALSE, algorithm = 3 )
DimPlot(object = atacs, reduction = "humap",label = TRUE) + NoLegend()

#map hematomap
hmap<-readRDS('outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds')
DefaultAssay(hmap)<-"integrated"

# quantify gene activity

gene.activities <- GeneActivity(atacs, features = VariableFeatures(hmap))

# add gene activities as a new assay
atacs[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atacs) <- "ACTIVITY"
atacs <- SCTransform(atacs,residual.features = rownames(atacs),assay = "ACTIVITY" )
# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = hmap, query = atacs, features = VariableFeatures(object = hmap),
    reference.assay = "integrated", query.assay = "ACTIVITY", reduction = "cca",normalization.method = "SCT")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = hmap$lineage,
    weight.reduction = atacs[["lsi"]], dims = 2:30)

atacs <- AddMetaData(atacs, metadata = celltype.predictions)
head(atacs[[]])
DimPlot(atacs, group.by = "predicted.id", label = TRUE,reduction = "humap")
saveRDS(atacs,fp(out,"cbps_atacs.rds"))

VlnPlot(atacs,"prediction.score.max",group.by="predicted.id")


#OCR by lineage
#run 14A - OCRby lineaea
atacs<-readRDS(fp(out,"cbps_atacs.rds"))
#need macs2

renv::use_python()
#reticulate::install_miniconda()
reticulate::use_miniconda()
reticulate::py_install(packages ="MACS2")

peaks <- CallPeaks(
  object = atacs,
  group.by = "predicted.id",
  macs2.path = "renv/python/virtualenvs/renv-python-3.7.3/bin/macs2"
)
head(peaks)
saveRDS(peaks,fp(out,"cbps_atacs_peaks_by_lineage.rds"))
peaks_dt<-data.table(as.data.frame(peaks))
peaks_dt[,peaks:=paste(seqnames,start,end,sep="-")]
peaks_dt[,chr:=seqnames]
fwrite(peaks_dt,fp(out,"cbps_atacs_peaks_by_lineage.csv.gz"))
peaks_dt2<-Reduce(rbind,lapply(unique(atacs$predicted.id), function(x){
  ocrs<-peaks_dt[str_detect(peak_called_in,x)]
  return(ocrs[,lineage:=x][,-'peak_called_in'])
  }))
fwrite(peaks_dt2,fp(out,"cbps_atacs_peaks_by_lineage2.csv.gz"))
#find lin spe OCR overlapping DMCs

cpgs_in_ocrs_lin<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          peaks_dt[,.(chr,start,end,peaks)][order(chr,start)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peaks","start","end"))

length(unique(cpgs_in_ocrs_lin$cpg_id))#275k/750k cpgs
length(unique(cpgs_in_ocrs_lin$peaks)) #66k/215k peaks

fwrite(cpgs_in_ocrs_lin,fp(out,"cpgs_in_lin_OCRs.csv.gz"))


res_cpgs_ocrs_lin<-merge(res_meth[,.(cpg_id,logFC,AveExpr,P.Value,adj.P.Val)],cpgs_in_ocrs_lin)
res_cpgs_ocrs_lin<-merge(res_cpgs_ocrs_lin,peaks_dt[,.(peaks,peak_called_in)],by="peaks")

res_cpgs_ocrs_lin2<-Reduce(rbind,lapply(unique(atacs$predicted.id), function(x){
  ocrs<-res_cpgs_ocrs_lin[str_detect(peak_called_in,x)]
  return(ocrs[,lineage:=x][,-'peak_called_in'])
  }))

fwrite(res_cpgs_ocrs_lin2,fp(out,"res_cpgs_in_lin_OCRs.csv.gz"))
res_cpgs_ocrs_lin2[,DMC:=P.Value<0.001&abs(logFC)>25]
ggplot(res_cpgs_ocrs_lin2)+geom_bar(aes(x=lineage,fill=lineage))
ggplot(res_cpgs_ocrs_lin2[DMC==T])+geom_bar(aes(x=lineage,fill=lineage))

fwrite(res_cpgs_ocrs_lin2[DMC==T],fp(out,"DMCs_in_lin_OCRs.csv.gz"))

res_cpgs_ocrs_lin2[,peak_size:=end-start]
ggplot(res_cpgs_ocrs_lin2)+geom_density(aes(peak_size))

#enrichment for dmcs/cpgs containing OCRs in lineage OCRs
ocrs_list<-split(peaks_dt2$peaks,peaks_dt2$lineage)
meth_list<-list(dmcs=unique(res_cpgs_ocrs_lin2[DMC==T]$peaks),
                cpgs=unique(res_cpgs_ocrs_lin2$peaks))

res_or_ocr<-OR2(querys =meth_list ,
            terms_list =ocrs_list,
            size_universe = nrow(peaks_dt),
            overlap_column = FALSE )
ggplot(res_or_ocr[term!="18"])+geom_col(aes(x=term,y=precision,fill=query),position = "dodge")

#enrichment for lineage OCR overlapping dmcs in DMCs
dmcs<-res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id

ocr_cpgs<-split(res_cpgs_ocrs_lin2$cpg_id,res_cpgs_ocrs_lin2$lineage)

res_or_cpgs<-OR2(querys =dmcs ,
            terms_list =ocr_cpgs,
            size_universe = nrow(res_meth),
            overlap_column = FALSE )
res_or_cpgs#33% des cpgs sont dans des OCRs HSC, alors que 72% des DMCs sont dans des OCRs HSC
ggplot(res_or_cpgs[term!="18"])+geom_col(aes(x=term,y=n.overlap,fill=precision),position = "dodge")


peaks

FeatureMatrix(peaks)

#II) methylOCR regul EGR1 network ?


