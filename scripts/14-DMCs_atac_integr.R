#hypermeth on OCR regulating EGR1 network ?
out<-'outputs/14-DMCs_atac_integr'
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
#renv::install("Signac")
library(Signac)

#I) BULK OCRs
#1)DMCs overlap in OCR ?
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
#for EGR1
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

#for all
genes_close_tf_meth<-Reduce(rbind,lapply(c("EGR1","KLF2","KLF4"),function(x){
  tf_motif_dmcs<-motifs_dmcs[,enriched.motifs[motif.name==x]$motif,drop=F]
  peaks_dmcs_tf<-rownames(tf_motif_dmcs)[as.vector(tf_motif_dmcs==1)]
  dt_genes_dmcs_tf<-data.table(ClosestFeature(atacs,regions = peaks_dmcs_tf))
  dt_genes_dmcs_tf[,n.gene:=length(unique(gene_id))]
  dt_genes_dmcs_tf[,dmc_region:=T]
  return(dt_genes_dmcs_tf[,tf.motif:=x])
  
  }))
fwrite(genes_close_tf_meth,fp(out,"res_closest_genes_tfmotif_dmcs_ocr.csv"))

#genes of the EGRN / EGR1/KLF2/KLF4 targets, xxx DEGs ? xxx Downregules ?

tftargets<-fread("outputs/13-GRN_integr/tf_target_interactions.csv")
EGRn<-fread("outputs/13-GRN_integr/egr1_network_tf_target_interactions.csv")
res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")

trs_list<-list(EGRN_target=unique(EGRn$target),
               EGR1_target=EGRn[tf=="EGR1"]$target,
               KLF4_target=EGRn[tf=="KLF4"]$target,
               KLF2_target=EGRn[tf=="KLF2"]$target)

trs_de_list<-lapply(trs_list,function(x)intersect(x,res_degs[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5]$gene))
names(trs_de_list)<-paste0(names(trs_list),"_degs")

trs_dn_list<-lapply(trs_list,function(x)intersect(x,res_degs[lineage=="HSC"&padj<0.05&log2FoldChange<(-0.5)]$gene))
names(trs_dn_list)<-paste0(names(trs_list),"_down")

degs_list<-list(degs=res_degs[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0.5]$gene,
                down=res_degs[lineage=="HSC"&padj<0.05&log2FoldChange<(-0.5)]$gene,
                up=res_degs[lineage=="HSC"&padj<0.05&log2FoldChange>0.5]$gene)
all_trs_list<-c(trs_list,trs_de_list,trs_dn_list,degs_list)

genes_list<-split(genes_close_tf_meth$gene_name,genes_close_tf_meth$tf.motif)
genes_list<-lapply(genes_list, unique)

res_or_trs<-OR3(all_trs_list,genes_list,background = unique(res_degs$gene)) 
res_or_trs[padj<0.05]
#          query term term.size n.query n.overlap pct.query.overlap precision pct.term.overlap background_size pct.term.background
# 1: KLF2_target EGR1      1433      98        18         0.1836735 0.1836735       0.01256106           12877           0.1112837
# 2: KLF2_target KLF2      1479      98        18         0.1836735 0.1836735       0.01217039           12877           0.1148559
# 3: KLF2_target KLF4      1598      98        19         0.1938776 0.1938776       0.01188986           12877           0.1240972
# 4:        down EGR1      1433     285        44         0.1543860 0.1543860       0.03070482           12877           0.1112837
#          pval       padj
# 1: 0.02198625 0.03104840
# 2: 0.02915675 0.03104840
# 3: 0.03104840 0.03104840
# 4: 0.01523962 0.04571885
#                                                                                                                                                                                                                                                                                        genes.overlap
# 1:                                                                                                                                                                                      AHNAK|ATF3|CD151|DDIT4|ETS1|HEXIM1|ID3|INTS6|KLF13|PTGER4|SOCS3|TINAGL1|TOB1|TRIM8|TSC22D1|UBC|ZFP36L1|IFRD1
# 2:                                                                                                                                                                                      AHNAK|ATF3|CD151|DDIT4|ETS1|HEXIM1|ID3|INTS6|KLF13|PTGER4|SOCS3|TINAGL1|TOB1|TRIM8|TSC22D1|UBC|ZFP36L1|IFRD1
# 3:                                                                                                                                                                              AHNAK|ATF3|CD151|DDIT4|ETS1|HEXIM1|ID3|INTS6|KLF13|PTGER4|SERTAD1|SOCS3|TINAGL1|TOB1|TRIM8|TSC22D1|UBC|ZFP36L1|IFRD1
# 4: C1orf109|TAF5L|ATP1B1|EFNA3|BEND5|IER5|LUZP1|IFFO2|SLC30A1|PLK3|SEMA4A|FAM43A|BAMBI|RPP38|SLC3A2|STIP1|C2CD2L|PHLDA2|C11orf96|CRY1|UBE2D1|NFKBIA|BRMS1L|EFNB2|RAB3IP|GCH1|DUSP5|JDP2|TAF5|PPRC1|BAG3|IKBIP|TNFRSF12A|PLEKHH3|CFAP20|DERL2|PSMC3IP|CYP1A1|SPATA2L|ATP6V0D1|SOCS3|HEXIM1|TOP1|H2AFX
#          query
# 1: KLF2_target
# 2: KLF2_target
# 3: KLF2_target
# 4:        down

fwrite(res_or_trs,fp(out,"res_trs_network_enrichment_for_closest_genes_of_dmcs_containing_tf_bulk_ocrs.csv"))


#II) OCR by lineage
#0) compute lin OCRs
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


#1) DMCs enrichment in linOCRs
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


#2) lin spe methylOCR regul EGR1 network ?
#2a) For HSC peaks
#feature matrix #[to do]

#assay HSC_peaks

#TF motif enrichment in dmcs peaks #[to update]
peaks_with_dmcs<-unique(res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25]$peaks)
peaks_without<-unique(res_cpgs_ocrs[!peaks%in%peaks_with_dmcs]$peaks)
enriched.motifs <- FindMotifs(
  object = atacs,
  features = peaks_with_dmcs ,
  background=peaks_without)


enriched.motifs<-data.table(enriched.motifs)
fwrite(enriched.motifs,fp(out,"enriched_motifs_in_CpGs_overlaping_OCRs_with_DMCs_vs_without.csv"))



#genes of the EGRN / EGR1/KLF2/KLF4 targets, xxx DEGs ? xxx Downregules ?

genes_close_tf_meth<-Reduce(rbind,lapply(c("EGR1","KLF2","KLF4"),function(x){
  tf_motif_dmcs<-motifs_dmcs[,enriched.motifs[motif.name==x]$motif,drop=F]
  peaks_dmcs_tf<-rownames(tf_motif_dmcs)[as.vector(tf_motif_dmcs==1)]
  dt_genes_dmcs_tf<-data.table(ClosestFeature(atacs,regions = peaks_dmcs_tf))
  dt_genes_dmcs_tf[,n.gene:=length(unique(gene_id))]
  dt_genes_dmcs_tf[,dmc_region:=T]
  return(dt_genes_dmcs_tf[,tf.motif:=x])
  
  }))
fwrite(genes_close_tf_meth,fp(out,"res_closest_genes_tfmotif_dmcs_ocr.csv"))

genes_list<-split(genes_close_tf_meth$gene_name,genes_close_tf_meth$tf.motif)
genes_list<-lapply(genes_list, unique)

res_or_trs<-OR3(all_trs_list,genes_list,background = unique(res_degs$gene)) 
res_or_trs[padj<0.05]
#          query term term.size n.query n.overlap pct.query.overlap precision pct.term.overlap background_size pct.term.background
# 1: KLF2_target EGR1      1433      98        18         0.1836735 0.1836735       0.01256106           12877           0.1112837
# 2: KLF2_target KLF2      1479      98        18         0.1836735 0.1836735       0.01217039           12877           0.1148559
# 3: KLF2_target KLF4      1598      98        19         0.1938776 0.1938776       0.01188986           12877           0.1240972
# 4:        down EGR1      1433     285        44         0.1543860 0.1543860       0.03070482           12877           0.1112837
#          pval       padj
# 1: 0.02198625 0.03104840
# 2: 0.02915675 0.03104840
# 3: 0.03104840 0.03104840
# 4: 0.01523962 0.04571885
#                                                                                                                                                                                                                                                                                        genes.overlap
# 1:                                                                                                                                                                                      AHNAK|ATF3|CD151|DDIT4|ETS1|HEXIM1|ID3|INTS6|KLF13|PTGER4|SOCS3|TINAGL1|TOB1|TRIM8|TSC22D1|UBC|ZFP36L1|IFRD1
# 2:                                                                                                                                                                                      AHNAK|ATF3|CD151|DDIT4|ETS1|HEXIM1|ID3|INTS6|KLF13|PTGER4|SOCS3|TINAGL1|TOB1|TRIM8|TSC22D1|UBC|ZFP36L1|IFRD1
# 3:                                                                                                                                                                              AHNAK|ATF3|CD151|DDIT4|ETS1|HEXIM1|ID3|INTS6|KLF13|PTGER4|SERTAD1|SOCS3|TINAGL1|TOB1|TRIM8|TSC22D1|UBC|ZFP36L1|IFRD1
# 4: C1orf109|TAF5L|ATP1B1|EFNA3|BEND5|IER5|LUZP1|IFFO2|SLC30A1|PLK3|SEMA4A|FAM43A|BAMBI|RPP38|SLC3A2|STIP1|C2CD2L|PHLDA2|C11orf96|CRY1|UBE2D1|NFKBIA|BRMS1L|EFNB2|RAB3IP|GCH1|DUSP5|JDP2|TAF5|PPRC1|BAG3|IKBIP|TNFRSF12A|PLEKHH3|CFAP20|DERL2|PSMC3IP|CYP1A1|SPATA2L|ATP6V0D1|SOCS3|HEXIM1|TOP1|H2AFX
#          query
# 1: KLF2_target
# 2: KLF2_target
# 3: KLF2_target
# 4:        down

fwrite(res_or_trs,fp(out,"res_trs_network_enrichment_for_closest_genes_of_dmcs_containing_tf_bulk_ocrs.csv"))

#enrichment for OCRs/ DMCs OCRs / EGR1 OCRSs / DMCs EGR1 OCRs for HSC opening ?

