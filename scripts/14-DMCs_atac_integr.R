#hypermeth on OCR regulating EGR1 network ?
out<-'outputs/14-DMCs_atac_integr'
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
#renv::install("Signac")
library(Signac)

#I) BULK OCRs####
#1)DMCs overlap in OCR ?
atacs<-readRDS("../atac/outputs/cbps_merged/cbps_atac1-4_merged_qc.rds")
nrow(atacs) #126925
ocrs<-data.table(peaks=rownames(atacs))
ocrs[,chr:=str_extract(peaks,"chr[0-9XY]+"),]
ocrs[,start:=strsplit(peaks,"-")[[1]][2],by="peaks"]
ocrs[,end:=strsplit(peaks,"-")[[1]][3],by="peaks"]

#QC peaks :
#ChromHMM anno dans peaks ?
chrine_feat<-fread("ref/Chromatin_Annot/CBP/CD34_all_chromatin_feature.csv")
#trans in hg38
chrine38<-hg19to38(chrine_feat[,.(chr,start,end,chrReg)][order(chr,start)])
chrine38[,chrReg:=id]
chrine38<-merge(chrine38[,-"id"],chrine_feat[,.(chrReg,type)],by="chrReg")
chrine38[is.na(start)] #ok
fwrite(chrine38,fp(out,"CD34_all_chromatin_feature_hg38.cs.gz"))


peaks_chrine<-bed_inter(a=ocrs[,.(chr,start,end,peaks)][order(chr,start)],
          b=chrine38[,.(chr,start,end,type)][order(chr,start)],
          select = c(1,2,3,4,6,7,8),col.names = c("chr","start.ocrs","end.ocrs","peaks","start.chrine","end.chrine","chrine_feat"))
peaks_chrine[,size_overlap:=min(c(end.ocrs,end.chrine))-max(c(start.ocrs,start.chrine)),by=.(peaks,chrine_feat)]

peaks_chrine[,bigger.overlap:=size_overlap==max(size_overlap),by="peaks"]
fwrite(peaks_chrine,fp(out,"peaks_chrine_feature_anno.csv.gz"))
peaks_chrine[,chrine_feat:=as.factor(chrine_feat)]
p1<-ggplot(peaks_chrine[bigger.overlap==T])+geom_bar(aes(x=chrine_feat,fill=chrine_feat))+ggtitle("Peaks-overlapping chromHMM regions")
chrine_feat[,type:=as.factor(type)]
chrine_feat[,chrine_feat:=type]
p2<-ggplot(chrine_feat)+geom_bar(aes(x=chrine_feat,fill=chrine_feat))+ggtitle("All chromHMM Regions")
p2+p1+plot_layout(guides = "collect")


#peaks overlap with cpgs
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
res_cpgs<-fread(fp(out,"res_cpgs_hg38.cs.gz"))

cpgs_in_ocrs<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          ocrs[,.(chr,start,end,peaks)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peaks","start","end"))

length(unique(cpgs_in_ocrs$cpg_id))#300k/750k
fwrite(cpgs_in_ocrs,fp(out,"cpgs_in_bulk_OCRs.csv.gz"))
cpgs_in_ocrs<-fread(fp(out,"cpgs_in_bulk_OCRs.csv.gz"))

res_cpgs_ocrs<-merge(res_meth[,.(cpg_id,logFC,AveExpr,P.Value,adj.P.Val)],cpgs_in_ocrs)

res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25] #3910
fwrite(res_cpgs_ocrs,fp(out,"res_meth_in_bulk_OCRs.csv.gz"))
res_cpgs_ocrs<-fread(fp(out,"res_meth_in_bulk_OCRs.csv.gz"))
fwrite(res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25],fp(out,"DMCs_in_bulk_OCRs.csv.gz"))
res_meth[P.Value<0.001&abs(logFC)>25] #3910 / 4815 DMCs

#QC 
#n DMCs / peaks
res_cpgs_ocrs[,dmc:=P.Value<0.001&abs(logFC)>25]
res_cpgs_ocrs[,count:=.N,by=.(dmc,peaks)]
ggplot(res_cpgs_ocrs)+geom_density(aes(x=count,fill=dmc))

ggplot(res_cpgs_ocrs[dmc==T])+geom_histogram(aes(x=count),bins=7,col="white")

res_cpgs_ocrs[,peak_size:=end-start]
ggplot(res_cpgs_ocrs)+geom_density(aes(peak_size))

#DMCs peaks in 4/6 regions ?
res_cpgs_ocrs_anno<-merge(res_cpgs_ocrs,peaks_chrine[,.(peaks,chrine_feat)],allow.cartesian=T)
ggplot(res_cpgs_ocrs_anno)+geom_bar(aes(x=dmc,fill=chrine_feat),position = "fill")+scale_y_continuous(labels = scales::percent)

#DMCs close TSS ?
cpgs_in_peaks<-unique(res_cpgs_ocrs$cpg_id)
dmcs_in_peaks<-unique(res_cpgs_ocrs[dmc==T]$cpg_id)

cpgs_anno<-fread("outputs/02A-CpGs_annotations/cpgs_annot.csv.gz")
cpgs_anno_tss<-cpgs_anno[in_eQTR==F]
cpgs_anno_tss[,in_peaks:=cpg_id%in%cpgs_in_peaks]
cpgs_anno_tss[,dmc:=cpg_id%in%res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id]
cpgs_anno_tss[,dmc_in_peaks:=dmc&in_peaks]
p1<-ggplot(cpgs_anno_tss)+geom_density(aes(x=tss_dist))+ggtitle("all cpgs")
p2<-ggplot(cpgs_anno_tss[in_peaks==T])+geom_density(aes(x=tss_dist))+ggtitle("cpgs in peaks")
p3<-ggplot(cpgs_anno_tss[dmc==T])+geom_density(aes(x=tss_dist))+ggtitle("DMCs")
p4<-ggplot(cpgs_anno_tss[dmc_in_peaks==T])+geom_density(aes(x=tss_dist))+ggtitle("DMCs in peaks")
(p1+p2)/(p3+p4)
cpgs_anno_tss[,cpg_type:=ifelse(dmc_in_peaks==T,"DMC in peak",ifelse(in_peaks==T,"CpG in peak",ifelse(dmc==T,"DMC not in peak","CpG not in peak")))]
ggplot(cpgs_anno_tss)+geom_density(aes(x=log10(abs(tss_dist)+1),fill=cpg_type))+facet_wrap("cpg_type")
ggplot(cpgs_anno_tss)+geom_density(aes(x=abs(tss_dist)+1,fill=cpg_type))+facet_wrap("cpg_type")+scale_x_log10()

#2) TF motif enrichment in this peaks
#renv::install("bioc::JASPAR2020")
#renv::install("bioc::TFBSTools")

library(JASPAR2020)
library(TFBSTools)

#renv::install("bioc::BSgenome.Hsapiens.UCSC.hg38")
#renv::install("bioc::motifmatchr")

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
enriched.motifs<-fread(fp(out,"enriched_motifs_in_CpGs_overlaping_OCRs_with_DMCs_vs_without.csv"))
head(enriched.motifs,25)
enriched.motifs[pvalue<0.000001]$motif.name
enriched.motifs[,motif.name:=factor(motif.name,levels=enriched.motifs[order(pvalue)]$motif.name)]

ggplot(enriched.motifs[order(pvalue)][1:25])+
  geom_point(aes(x=motif.name,size=fold.enrichment,y=observed,col=-log10(pvalue+10^-316)))

MotifPlot(
  object = atacs,
  motifs = head(enriched.motifs$motif)
)


#3) Downregulated genes enriched for EGR1 + methylated HSC peaks ? compared to EGR1 peaks, methylated peaks
#get EGR1 peaks and closest genes 
motif.all <- GetMotifData(
    object = atacs, assay = "ATAC", slot = "data"
  )
enriched.motifs<-fread(fp(out,"enriched_motifs_in_CpGs_overlaping_OCRs_with_DMCs_vs_without.csv"))

egr1_peaks<-motif.all[,enriched.motifs[motif.name=="EGR1"]$motif,drop=F]
#sum(egr1_peaks>1)
egr1_peaks<-rownames(egr1_peaks)[as.vector(egr1_peaks==1)]
length(egr1_peaks) #22k
egr1_genes_dt<-data.table(ClosestFeature(atacs,regions = egr1_peaks))
fwrite(egr1_genes_dt,fp(out,"genes_of_egr1_containing_bulk_peaks.csv.gz"))
egr1_genes<-unique(egr1_genes_dt$gene_name)
length(egr1_genes) #13k

#get methylated peaks
dmcs_ocrs<-fread(fp(out,"DMCs_in_bulk_OCRs.csv.gz"))
dmcs_peaks<-unique(dmcs_ocrs$peaks)
length(dmcs_peaks) #3170
dmcs_genes_dt<-data.table(ClosestFeature(atacs,regions = dmcs_peaks))
fwrite(dmcs_genes_dt,fp(out,"genes_of_dmcs_containing_bulk_peaks.csv.gz"))
dmcs_genes<-unique(dmcs_genes_dt$gene_name)
length(dmcs_genes) #2854


#EGR1 + methylated peaks/genes
egr1_dmcs_peaks<-intersect(egr1_peaks,dmcs_peaks)
length(egr1_dmcs_peaks) #2363

egr1_dmcs_genes<-intersect(egr1_genes,dmcs_genes)
length(egr1_dmcs_genes) #2494

#downregulated genes enriched for genes close to egr1 dmcs peaks ?

res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")
degs_list<-list(degs=res_degs[lineage=="HSC"&padj<0.05&abs(log2FoldChange)>0]$gene,
                down=res_degs[lineage=="HSC"&padj<0.05&log2FoldChange<(-0.5)]$gene,
                up=res_degs[lineage=="HSC"&padj<0.05&log2FoldChange>0.5]$gene)

egr1_met_list<-list(egr1=egr1_genes,
                meth=dmcs_genes,
                egr1_meth=egr1_dmcs_genes)

res_or<-OR3(degs_list,egr1_met_list,background = unique(res_degs$gene)) 
fwrite(res_or,fp(out,"res_degs_enrichment_for_egr1_and_or_dmcs_peaks.csv"))

res_or<-fread(fp(out,"res_degs_enrichment_for_egr1_and_or_dmcs_peaks.csv"))
ggplot(res_or[,.SD,.SDcols=-ncol(res_or)])+geom_point(aes(x=term,y=-log10(pval),size=n.overlap,col=pct.term.overlap))+facet_wrap("query")


length(intersect(egr1_genes,res_degs$gene))
length(unique(res_degs$gene))#8k/13k ac EGR1 motif

#down of the EGRN enriched for egr1 dmcs peaks ?
egrn<-fread("outputs/13-GRN_integr/egr1_network_tf_target_interactions.csv")
egrn_genes<-unique(egrn$target)
length(egrn_genes) #208
degs_egrn_list<-lapply(degs_list, function(x)intersect(x,egrn_genes))
# $degs
#  [1] "ATP1B1"  "DUSP2"   "IER5"    "ABHD5"   "MAP9"    "TNFAIP3" "EGR1"    "CDKN1A"  "DNAJB9" 
# [10] "IER3"    "BRD2"    "HSPA1B"  "HSPA1A"  "DNAJA1"  "GADD45G" "NFKBIA"  "UBL3"    "PNP"    
# [19] "ARIH1"   "SOCS3"   "HEXIM1"  "GADD45B" "ID1"     "KLF2"    "ZFP36"   "NFATC1"  "DNAJB1" 
# [28] "IER2"    "SERTAD3" "JUNB"    "MYADM"   "TOB2"   
# 
# $down
#  [1] "ATP1B1"  "DUSP2"   "IER5"    "ABHD5"   "MAP9"    "TNFAIP3" "EGR1"    "CDKN1A"  "DNAJB9" 
# [10] "IER3"    "BRD2"    "HSPA1B"  "HSPA1A"  "DNAJA1"  "GADD45G" "NFKBIA"  "UBL3"    "PNP"    
# [19] "ARIH1"   "SOCS3"   "HEXIM1"  "GADD45B" "ID1"     "KLF2"    "ZFP36"   "NFATC1"  "DNAJB1" 
# [28] "IER2"    "SERTAD3" "JUNB"    "MYADM"   "TOB2"   
# 
# $up
# character(0)


res_or2<-OR3(degs_egrn_list,egr1_met_list,background = unique(res_degs$gene)) 
fwrite(res_or2,fp(out,"res_egrn_degs_enrichment_for_egr1_and_or_dmcs_peaks.csv"))
res_or2<-fread(fp(out,"res_egrn_degs_enrichment_for_egr1_and_or_dmcs_peaks.csv"))

ggplot(res_or2[,.SD,.SDcols=-ncol(res_or2)])+geom_point(aes(x=term,y=-log10(pval),size=n.overlap,col=pct.term.overlap))+facet_wrap("query")
res_orm<-rbind(res_or[,.SD,.SDcols=-ncol(res_or)][,genes:="All DEGs"],res_or2[,.SD,.SDcols=-ncol(res_or2)][,genes:="DEGs of the EGR1 network"])
ggplot(res_orm)+geom_point(aes(x=term,y=-log10(pval),size=n.overlap,col=pct.term.overlap))+facet_grid(genes~query)


# peaks methylé with EGR1/KLF2/KLF4/.. motif close to which genes? 
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
genes_close_tf_meth<-fread(fp(out,"res_closest_genes_tfmotif_dmcs_ocr.csv"))


# genes of the EGRN ?  genes Downregulé ? Genes of the EGRN downregules ?
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

ggplot(res_or_trs[,.SD,.SDcols=-ncol(res_or_trs)])+geom_point(aes(x=term,y=-log10(pval),size=n.overlap,col=pct.term.overlap))+facet_wrap("query")

fwrite(res_or_trs,fp(out,"res_trs_network_enrichment_for_closest_genes_of_dmcs_containing_tf_bulk_ocrs.csv"))
res_or_trs<-fread(fp(out,"res_trs_network_enrichment_for_closest_genes_of_dmcs_containing_tf_bulk_ocrs.csv"))
res_or_trs[padj<0.1]


#II) ANalysis by lineage####

#1) need first annotate cells ~ lineage hmap
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

#2) compute lin OCRs
atacs<-readRDS(fp(out,"cbps_atacs.rds"))
#need macs2

renv::use_python()
#reticulate::install_miniconda()
#reticulate::use_miniconda()
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
peaks_dt<-fread(fp(out,"cbps_atacs_peaks_by_lineage.csv.gz"))


peaks_dt2<-Reduce(rbind,lapply(str_replace(unique(atacs$predicted.id),"/","_"), function(x){
  ocrs<-peaks_dt[str_detect(peak_called_in,x)]
  return(ocrs[,lineage:=x][,-'peak_called_in'])
  }))
peaks_dt2[,n.lin:=.N,by="peaks"]
table(peaks_dt2[n.lin==1]$lineage)
    #  18          DC Erythro-Mas         HSC    Lymphoid       Mk_Er    MPP_LMPP 
    #    4608        2138       12377        1470        3975           2       41604 
    # Myeloid 
    #    1567 
fwrite(peaks_dt2,fp(out,"cbps_atacs_peaks_by_lineage2.csv.gz"))
peaks_dt2<-fread(fp(out,"cbps_atacs_peaks_by_lineage2.csv.gz"))

#3) cpgs overlapping peaks
cpgs_in_ocrs_lin<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          peaks_dt[,.(chr,start,end,peaks)][order(chr,start)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peaks","start","end"))

length(unique(cpgs_in_ocrs_lin$cpg_id))#275k/750k cpgs
length(unique(cpgs_in_ocrs_lin$peaks)) #66k/215k peaks

fwrite(cpgs_in_ocrs_lin,fp(out,"cpgs_in_lin_OCRs.csv.gz"))

#% DMCs overlapping linOCR
cpgs_in_ocrs_lin<-fread(fp(out,"cpgs_in_lin_OCRs.csv.gz"))
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")

peaks_meth<-merge(cpgs_in_ocrs_lin,res_meth,by="cpg_id")
nrow(peaks_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]) #74%

#vs % CpGs overlapping 
nrow(peaks_meth)/nrow(res_meth) #36%

(nrow(peaks_meth[P.Value<0.001&abs(logFC)>25])/nrow(res_meth[P.Value<0.001&abs(logFC)>25]))/(nrow(peaks_meth)/nrow(res_meth)) #2 foldenrichment


OR(set1 = res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id,
   set2= peaks_meth$cpg_id,size_universe = nrow(res_meth))

#4) linspe OCRs (findmarkers)
#   a) new assays
peaks<-readRDS(fp(out,"cbps_atacs_peaks_by_lineage.rds"))
peaks_mat<-FeatureMatrix(Fragments(atacs),
                             features =peaks,
                         process_n = 2000,
                             cells=colnames(atacs))

saveRDS(peaks_mat,fp(out,"lin_spe_peaks_counts.rds"))
#Create Assay
atacs[["lin_peaks"]]<-CreateChromatinAssay(peaks_mat)
saveRDS(atacs@assays$lin_peaks,fp(out,"cbps_lin_spe_peaks_assay.rds"))


#   b) findallmarkers
DefaultAssay(atacs)<-"lin_peaks"
Idents(atacs)<-"predicted.id"
peaks_lin<-FindAllMarkers(atacs,
                          only.pos=T,
                          min.pct = 0.05,
                          test.use = "LR",
                          latent.vars = "peak_region_fragments")
peaks_lin<-data.table(peaks_lin)
table(peaks_lin[p_val_adj<0.001]$cluster)
  # MPP/LMPP    Lymphoid         HSC      B cell     Myeloid Erythro-Mas          DC          18      LT-HSC      T cell   Mk/Er 
  #         1         622         465        1574         359        1932        2030        2173           2         656   0

fwrite(peaks_lin,fp(out,"peaks_markers_lineage.csv.gz"))
peaks_lin<-fread(fp(out,"peaks_markers_lineage.csv.gz"))

#  c) DMCs enrichemnt for HSC spec peaks ?
cpgs_in_ocrs_lin<-fread(fp(out,"cpgs_in_lin_OCRs.csv.gz"))
length(unique(cpgs_in_ocrs_lin$peaks))
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
dmcs<-res_meth[P.Value<0.001&abs(logFC)>25]$cpg_id
lin_ofint<-c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")
lin_spe_cpgs<-lapply(lin_ofint,
                              function(lin){
                                ocrs_lin_spe<-peaks_lin[cluster==lin&p_val_adj<0.001]$gene
                                cpgs_lin_spe<-unique(cpgs_in_ocrs_lin[peaks%in%ocrs_lin_spe]$cpg_id)
                                return(cpgs_lin_spe)
                                })
names(lin_spe_cpgs)<-lin_ofint
lapply(lin_spe_cpgs, length)
# $HSC
# [1] 4786
# 
# $`MPP/LMPP`
# [1] 3
# 
# $`Erythro-Mas`
# [1] 2124
# 
# $Lymphoid
# [1] 1198
# 
# $Myeloid
# [1] 337
res_dmcs_enr<-OR3(querys =lin_spe_cpgs,
            terms_list =list(dmcs=dmcs),
            background = res_meth$cpg_id,
            overlap_column = FALSE )
#         query. term term.size n.query n.overlap pct.query.overlap   precision pct.term.overlap
# 1:         HSC dmcs      4815    4786        63       0.013163393 0.013163393     0.0130841121
# 2:    MPP/LMPP dmcs      4815       3         0       0.000000000 0.000000000     0.0000000000
# 3: Erythro-Mas dmcs      4815    2124        14       0.006591337 0.006591337     0.0029075805
# 4:    Lymphoid dmcs      4815    1198        13       0.010851419 0.010851419     0.0026998962
# 5:     Myeloid dmcs      4815     337         2       0.005934718 0.005934718     0.0004153686
#    background_size pct.term.background         pval         padj       query
# 1:          754931         0.006378066 1.444455e-07 1.444455e-07         HSC
# 2:          754931         0.006378066 1.000000e+00 1.000000e+00    MPP/LMPP
# 3:          754931         0.006378066 4.871420e-01 4.871420e-01 Erythro-Mas
# 4:          754931         0.006378066 4.735086e-02 4.735086e-02    Lymphoid
# 5:          754931         0.006378066 6.339120e-01 6.339120e-01     Myeloid
#yes !
ggplot(res_dmcs_enr)+geom_point(aes(x=query,size=n.overlap,y=-log10(pval),col=precision))

fwrite(res_dmcs_enr,fp(out,"res_dmcs_enrichment_for_lineage_specific_cpgs.csv"))

#same but at the peak level
#rq : distribution of n dmcs / peaks containing dmcs
table(as.vector(table(cpgs_in_ocrs_lin[cpg_id%in%dmcs]$peaks)))
#    1    2    3    4    5 
# 2584  392   63    9    3 
sum(c(2584,392,63,9,3)) #3051
#end rq
peaks_dmcs<-unique(cpgs_in_ocrs_lin[cpg_id%in%dmcs]$peaks)
lin_spe_peaks<-lapply(lin_ofint,
                              function(lin){
                                ocrs_lin_spe<-peaks_lin[cluster==lin&p_val_adj<0.001]$gene
                                return(ocrs_lin_spe)
                                })
names(lin_spe_peaks)<-lin_ofint
lapply(lin_spe_peaks, length)
# $HSC
# [1] 465
# 
# $`MPP/LMPP`
# [1] 1
# 
# $`Erythro-Mas`
# [1] 1932
# 
# $Lymphoid
# [1] 622
# 
# $Myeloid
# [1] 359

res_dmcs_peaks_enr<-OR3(querys = lin_spe_peaks,
            terms_list =list(peaks_dmcs=peaks_dmcs),
            background = unique(cpgs_in_ocrs_lin$peaks),
            overlap_column = FALSE )
res_dmcs_peaks_enr
#         query.       term term.size n.query n.overlap pct.query.overlap  precision
# 1:         HSC peaks_dmcs      3051     438        48        0.10958904 0.10958904
# 2:    MPP/LMPP peaks_dmcs      3051       1         0        0.00000000 0.00000000
# 3: Erythro-Mas peaks_dmcs      3051     722        13        0.01800554 0.01800554
# 4:    Lymphoid peaks_dmcs      3051     427        13        0.03044496 0.03044496
# 5:     Myeloid peaks_dmcs      3051     158         2        0.01265823 0.01265823
#    pct.term.overlap background_size pct.term.background         pval         padj
# 1:     0.0157325467           66196           0.0460904 3.519755e-08 3.519755e-08
# 2:     0.0000000000           66196           0.0460904 1.000000e+00 1.000000e+00
# 3:     0.0042608981           66196           0.0460904 9.999863e-01 9.999863e-01
# 4:     0.0042608981           66196           0.0460904 9.589900e-01 9.589900e-01
# 5:     0.0006555228           66196           0.0460904 9.950412e-01 9.950412e-01
#          query
# 1:         HSC
# 2:    MPP/LMPP
# 3: Erythro-Mas
# 4:    Lymphoid
# 5:     Myeloid

ggplot(res_dmcs_peaks_enr)+geom_point(aes(x=query,size=n.overlap,col=-log10(padj),y=pct.term.overlap))
fwrite(res_dmcs_peaks_enr,fp(out,"res_dmcs_peaks_enrichment_for_lineage_specific_peaks.csv"))

#annot 48 peaks HSC with DMCs
annotations<-readRDS("../atac/ref/gene_annotations_hg38_GRanges.rds")

# add the gene information to the object
Annotation(atacs) <- annotations

peaks_hsc_dmcs<-intersect(peaks_dmcs,lin_spe_peaks$HSC)
peaks_hsc_dmcs_anno<-ClosestFeature(atacs,peaks_hsc_dmcs)
peaks_hsc_dmcs_anno$gene_name
fwrite(peaks_hsc_dmcs_anno,fp(out,"anno_hsc_dmcs_peaks.csv"))



#  d) EGR1 enrichment in HSC spec peaks ?
# Get a list of motif position frequency matrices from the JASPAR database

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
atacs <- AddMotifs(
  object = atacs,assay = "lin_peaks",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

saveRDS(atacs@assays$lin_peaks@motifs,fp(out,"atacs_cbps_lin_peaks_motif_object.rds"))
peaks_lin_list<-split(peaks_lin[p_val_adj<0.001]$gene,peaks_lin[p_val_adj<0.001]$cluster)

res_motif_lin<-Reduce(rbind,lapply(lin_ofint, function(lin){
  enriched.motifs <- FindMotifs(object = atacs,assay = "lin_peaks",
                                features = peaks_lin_list[[lin]] )
  return(data.table(enriched.motifs)[,lineage:=lin])
  }))


fwrite(res_motif_lin,fp(out,"enriched_motifs_by_lineage_specific_peaks.csv"))

lapply(lin_ofint, function(lin)res_motif_lin[lineage==lin][1:20])

p_list<-lapply(lin_ofint, function(lin){
  res<-res_motif_lin[lineage==lin][1:6]
  p<-ggplot(res)+geom_point(aes(x=motif.name,size=fold.enrichment,y=observed,col=-log10(pvalue+10^-316)))+ggtitle(lin)
  return(p)
  })

patchwork::wrap_plots(p_list,ncol = 2)

res_motif_lin[,top6:=rank(pvalue)<=6,by="lineage"]
res_motif_lin[,pvalue:=pvalue+10^-316]
ggplot(res_motif_lin[top6==T])+
  geom_point(aes(x=motif.name,size=observed,col=-log10(pvalue),y=percent.observed))+
facet_wrap("lineage",scales = "free_x",ncol = 3)+
  scale_color_gradient2(low = "grey",mid = "red",high = "darkred",limits=c(1,317),midpoint = 150)+
  scale_y_continuous(limits=c(0,100))

res_motif_lin[,top6.hsc:=motif.name%in%motif.name[(top6)&lineage=="HSC"]]
ggplot(res_motif_lin[top6.hsc==T])+
  geom_point(aes(x=lineage,size=observed,y=-log10(pvalue),col=percent.observed))+
  scale_x_discrete(limits=lin_ofint)+
facet_wrap("motif.name")

res_motif_lin[lineage=='HSC'][1:20]

# DMCs enirchment in TF containing peaks 
res_dmcs_tf_enr_lin<-Reduce(rbind,lapply(lin_ofint, function(lin){
  motifs_lin <- motif.all[peaks_lin_list[[lin]], , drop = FALSE]
  top10tf<-res_motif_lin[lineage==lin][order(pvalue)][1:10]$motif.name
  
  tf_peaks_list<-lapply(top10tf,function(tf){
    tf_peaks<-motifs_hsc[,unique(res_motif_lin[motif.name==tf]$motif),drop=F]
    tf_peaks<-rownames(tf_peaks)[as.vector(tf_peaks>0)]
    return(tf_peaks)
    })
  names(tf_peaks_list)<-top10tf
  dmcs_peaks<-intersect(peaks_dmcs,peaks_lin_list[[lin]])
  
  res_dmcs_tf_enr<-OR3(querys =tf_peaks_list ,
            terms_list =list(dmcs=dmcs_peaks),
            background = peaks_lin_list[[lin]],
            overlap_column = FALSE )
  return(res_dmcs_tf_enr[,lineage:=lin])
  
  }))

res_dmcs_tf_enr_lin[padj<0.05]
#     query term term.size n.query n.overlap pct.query.overlap precision pct.term.overlap background_size pct.term.background         pval         padj
# 1:  KLF15 dmcs        48     343        45         0.1311953 0.1311953       0.93750000             465          0.10322581 0.0002445452 0.0002445452
# 2:    SP9 dmcs        48     333        42         0.1261261 0.1261261       0.87500000             465          0.10322581 0.0054758517 0.0054758517
# 3:  KLF16 dmcs        48     355        43         0.1211268 0.1211268       0.89583333             465          0.10322581 0.0131948449 0.0131948449
# 4: ZNF148 dmcs        48     395        45         0.1139241 0.1139241       0.93750000             465          0.10322581 0.0471975676 0.0471975676
# 5:    SP3 dmcs        48     324        41         0.1265432 0.1265432       0.85416667             465          0.10322581 0.0071404726 0.0071404726
# 6:   EGR3 dmcs        48     254        34         0.1338583 0.1338583       0.70833333             465          0.10322581 0.0119781305 0.0119781305
# 7:   EBF1 dmcs        13       2         1         0.5000000 0.5000000       0.07692308             622          0.02090032 0.0413967721 0.0413967721
# 8:  ASCL1 dmcs        13       2         1         0.5000000 0.5000000       0.07692308             622          0.02090032 0.0413967721 0.0413967721
#     query  lineage
# 1:  KLF15      HSC
# 2:    SP9      HSC
# 3:  KLF16      HSC
# 4: ZNF148      HSC
# 5:    SP3      HSC
# 6:   EGR3      HSC
# 7:   EBF1 Lymphoid
# 8:  ASCL1 Lymphoid

ggplot(res_dmcs_tf_enr_lin[,.SD,.SDcols=-(ncol(res_dmcs_tf_enr_lin)-1)])+
  geom_point(aes(x=query,size=n.overlap,y=-log10(pval),col=pct.term.overlap))+
facet_wrap("lineage",scales = "free_x")

#Leukemic marrow infiltration reveals a novel role for 
#Egr3 as a potent inhibitor of normal hematopoietic stem cell proliferation [https://ashpublications.org/blood/article/126/11/1302/34392/Leukemic-marrow-infiltration-reveals-a-novel-role]

res_dmcs_tf_enr_lin[query=="EGR1"]
fwrite(res_dmcs_tf_enr_lin,fp(out,"res_dmcs_enrichment_in_top10_lin_spe_tf_peaks.csv"))


#TF motif enrichment in HSC DMCs peaks vs CpGs peaks
#peaks_dt<-fread("outputs/14-DMCs_atac_integr/cbps_atacs_peaks_by_lineage2.csv.gz")
hsc_peaks<-peaks_dt[lineage=="HSC"]$peaks
peaks_dmcs<-unique(cpgs_in_ocrs_lin[cpg_id%in%dmcs]$peaks)
peaks_hsc_with_dmcs<-intersect(hsc_peaks,peaks_dmcs)
length(peaks_hsc_with_dmcs)#2940
peaks_with_cpgs<-unique(cpgs_in_ocrs_lin[!peaks%in%peaks_dmcs]$peaks)
peaks_hsc_with_cpgs<-intersect(peaks_with_cpgs,hsc_peaks)
length(peaks_hsc_with_cpgs)#49113
res_motif_hsc_dmcs<-FindMotifs(object = atacs,assay = "lin_peaks",
                                features = peaks_hsc_with_dmcs,
                               background = peaks_hsc_with_cpgs)

res_motif_hsc_dmcs<-data.table(res_motif_hsc_dmcs)
res_motif_hsc_dmcs[pvalue<10^-50]
fwrite(res_motif_hsc_dmcs,fp(out,"res_tf_motif_enrichment_in_DMCs_containing_vs_CpGs_containing_hsc_peaks.csv"))


#overlap DMCs containing and EGR1 containing HSC peaks ?
motif.all <- GetMotifData(
    object = atacs, assay = "lin_peaks", slot = "data"
  )
motifs_hsc <- motif.all[hsc_peaks$peaks, , drop = FALSE]
egr1_peaks_hsc<-motifs_hsc[,unique(res_motif_lin[motif.name=="EGR1"]$motif),drop=F]
#sum(egr1_peaks_hsc>0)
egr1_peaks_hsc<-rownames(egr1_peaks_hsc)[as.vector(egr1_peaks_hsc==1)]
length(egr1_peaks_hsc) #22815

dmcs_egr1_hsc_peaks<-intersect(peaks_dmcs,egr1_peaks_hsc)
length(dmcs_egr1_hsc_peaks) #2106



# DEGs integr####
# peaks containing DMCs enriched for DEGs ? peaks with EGR1/KLF2/KLF4 ? peaks with EGR1/KLF2/KLF4+DMCs ? 
#peaks_genes_anno
atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
DefaultAssay(atacs)<-"lin_peaks"
#annot peaks
annotations<-readRDS("../atac/ref/gene_annotations_hg38_GRanges.rds")
Annotation(atacs) <- annotations

peaks_dt<-fread(fp(out,"cbps_atacs_peaks_by_lineage2.csv.gz"))

peaks_hsc_genes<-data.table(ClosestFeature(atacs,regions =peaks_dt[lineage=="HSC"]$peaks ))
fwrite(peaks_hsc_genes,fp(out,"peaks_hsc_genes_anno.csv.gz"))

#merge with res_degs
res_degs<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")
peaks_hsc_genes[,gene:=gene_name]
peaks_expr<-merge(peaks_hsc_genes,res_degs[lineage=="HSC"],by = "gene")
peaks_expr#63959 peaks asso to express genes
length(unique(peaks_expr$gene)) #10157 genes

#merge with res_cpgs
peaks_cpgs<-fread("outputs/14-DMCs_atac_integr/cpgs_in_lin_OCRs.csv.gz")
peaks_meth<-merge(peaks_cpgs,fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz"))
length(unique(peaks_meth$peaks)) #66196 peaks
length(unique(peaks_meth$cpg_id)) #274516 CpGs

#merge with tf motif peaks
res_motif_lin<-fread('outputs/14-DMCs_atac_integr/res_tf_motif_enrichment_in_DMCs_containing_vs_CpGs_containing_hsc_peaks.csv')
atacs@assays$lin_peaks@motifs<-readRDS("outputs/14-DMCs_atac_integr/atacs_cbps_lin_peaks_motif_object.rds")
motif.all <- GetMotifData(
    object = atacs, assay = "lin_peaks", slot = "data"
  )

tfs<-c("EGR1","KLF2","KLF4")
peaks_tf<-Reduce(rbind,lapply(tfs, function(tf){
  tf_peaks<-motif.all[,unique(res_motif_lin[motif.name==tf]$motif),drop=F]
  tf_peaks<-rownames(tf_peaks_hsc)[as.vector(tf_peaks_hsc==1)]
  return(data.table(peaks=tf_peaks,motif=tf))
  }))
fwrite(peaks_tf,fp(out,"peaks_containing_egr1_klf2_klf4.csv.gz"))
peaks_tf_hsc<-peaks_tf[peaks%in%peaks_dt[lineage=="HSC"]$peaks]

#term list : peaks DMCs/Tf
meth_tf_peaks_list<-list(CpGs=unique(peaks_meth$peaks),
                         DMCs=unique(peaks_meth[P.Value<0.001&abs(logFC)>25]$peaks),
                         EGR1=peaks_tf_hsc[motif=="EGR1"]$peaks,
                         KLF2=peaks_tf_hsc[motif=="KLF2"]$peaks,
                         KLF4=peaks_tf_hsc[motif=="KLF4"]$peaks,
                         EGR1_DMCs=intersect(peaks_tf_hsc[motif=="EGR1"]$peaks,unique(peaks_meth[P.Value<0.001&abs(logFC)>25]$peaks)),
                         KLF2_DMCs=intersect(peaks_tf_hsc[motif=="KLF2"]$peaks,unique(peaks_meth[P.Value<0.001&abs(logFC)>25]$peaks)),
                         KLF4_DMCs=intersect(peaks_tf_hsc[motif=="KLF2"]$peaks,unique(peaks_meth[P.Value<0.001&abs(logFC)>25]$peaks)))
lapply(meth_tf_peaks_list, length)
lapply(meth_tf_peaks_list, head)

#query list : degs up dn
degs_list<-list(degs=unique(peaks_expr[padj<0.05&abs(log2FoldChange)>0.5]$query_region),
                down=unique(peaks_expr[padj<0.05&log2FoldChange<(-0.5)]$query_region),
                up=unique(peaks_expr[padj<0.05&log2FoldChange>0.5]$query_region))
lapply(degs_list, length)
lapply(degs_list, head)

#background : peaks with expr genes 
background<-peaks_expr$query_region
length(background) #63959/122506 peaks (52%)
head(background)

res_or_degs_dmcs_tf<-OR3(degs_list,
                    meth_tf_peaks_list,
                    background = background,
                    overlap_column = F) 

fwrite(res_or_degs_dmcs_tf,fp(out,"res_degs_enrichment_in_dmcs_tf_hsc_peaks.csv"))
res_or_degs_dmcs_tf[padj<0.05] #

ggplot(res_or_degs_dmcs_tf)+geom_col(aes(y=fold.enrichment,x=term,fill=query),position = "dodge")


#DEGs enrichment gene lvl
#genes for peaks DMCs
peaks_hsc_genes[,peaks:=query_region]
peaks_meth_genes<-merge(peaks_hsc_genes,peaks_meth,by = "peaks")
length(unique(peaks_meth_genes$peaks))#52k peaks asso to peaks DMCs
length(unique(peaks_meth_genes$gene)) #19k genes

#genes for peaks tf
peaks_tf_genes<-merge(peaks_hsc_genes,peaks_tf_hsc,by = "peaks")
length(unique(peaks_tf_genes$peaks))#52k peaks asso to peaks DMCs
length(unique(peaks_tf_genes$gene)) #19k genes

#querylist : genes DMCs/Tf
meth_tf_genes_list<-list(CpGs=unique(peaks_meth_genes$gene),
                         DMCs=unique(peaks_meth_genes[P.Value<0.001&abs(logFC)>25]$gene),
                         EGR1=unique(peaks_tf_genes[motif=="EGR1"]$gene),
                         KLF2=unique(peaks_tf_genes[motif=="KLF2"]$gene),
                         KLF4=unique(peaks_tf_genes[motif=="KLF4"]$gene),
                         EGR1_DMCs=intersect(unique(peaks_tf_genes[motif=="EGR1"]$gene),unique(peaks_meth_genes[P.Value<0.001&abs(logFC)>25]$gene)),
                         KLF2_DMCs=intersect(unique(peaks_tf_genes[motif=="KLF2"]$gene),unique(peaks_meth_genes[P.Value<0.001&abs(logFC)>25]$gene)),
                         KLF4_DMCs=intersect(unique(peaks_tf_genes[motif=="KLF2"]$gene),unique(peaks_meth_genes[P.Value<0.001&abs(logFC)>25]$gene)))
lapply(meth_tf_genes_list, length)
lapply(meth_tf_genes_list, head)

#term list : degs up dn
degs_list<-list(degs=unique(peaks_expr[padj<0.05&abs(log2FoldChange)>0.5]$gene),
                down=unique(peaks_expr[padj<0.05&log2FoldChange<(-0.5)]$gene),
                up=unique(peaks_expr[padj<0.05&log2FoldChange>0.5]$gene))
lapply(degs_list, length)
lapply(degs_list, head)

#background :  expr genes 
background<-unique(peaks_expr$gene)
length(background) #10k
head(background)

res_or_degs_dmcs_tf_genes<-OR3(meth_tf_genes_list,
                         degs_list,
                    background = background,
                    overlap_column = F) 

fwrite(res_or_degs_dmcs_tf_genes,fp(out,"res_degs_enrichment_in_dmcs_tf_hsc_genes.csv"))
res_or_degs_dmcs_tf_genes[padj<0.05] #not sig


#GO analysis of HSC peaks containing DMCs
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
peaks_hsc_genes<-fread(fp(out,"peaks_hsc_genes_anno.csv.gz"))

cpgs_in_ocrs_lin<-fread(fp(out,"cpgs_in_lin_OCRs.csv.gz"))
res_meth<-fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz")
res_meth[P.Value<0.001&abs(logFC)>25]
peaks_meth<-merge(cpgs_in_ocrs_lin,res_meth,by="cpg_id")
peaks_hsc_genes[,peaks:=query_region]
peaks_hsc_genes_meth<-merge(peaks_hsc_genes,peaks_meth,by="peaks")
unique(peaks_hsc_genes_meth[P.Value<0.001&abs(logFC)>25]$gene_name)

res_go_bp<-enrichGO(bitr(unique(peaks_hsc_genes_meth[P.Value<0.001&abs(logFC)>25]$gene_name),
                         fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID,ont = "BP",
                 OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1,maxGSSize = 800,
                 universe =bitr(unique(peaks_hsc_genes_meth$gene_name),fromType = "SYMBOL",
                             toType = "ENTREZID",OrgDb = org.Hs.eg.db)$ENTREZID)

res_go_dt<-data.table(as.data.frame(res_go_bp))
res_go_dt[p.adjust<0.05]#36
res_go_dt[p.adjust<0.1]$Description
emapplot(pairwise_termsim(res_go_bp),showCategory = 36,cex_label_category=0.66)
res_go_dt[Description=="negative regulation of growth"]
res_go_dt[Description=="regulation of growth"] #nop
saveRDS(res_go_bp,fp(out,"res_go_bp_dmcs_hsc_peaks.rds"))
fwrite(res_go_dt,fp(out,"res_go_bp_dmcs_hsc_peaks.csv.gz"))

#gsego

peaks_hsc_genes_meth[P.Value<0.001&abs(logFC)>25,dmcscore.peak:=max(-log10(P.Value)*abs(logFC)),by=.(peaks)]
peaks_hsc_genes_meth[P.Value<0.001&abs(logFC)>25,dmcscore.gene:=mean(dmcscore.peak),by=.(gene_name)]
res_dt1<-merge(peaks_hsc_genes_meth,
              data.table(bitr(unique(peaks_hsc_genes_meth$gene_name),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db,drop=F))[,gene_name:=SYMBOL],by="gene_name")[!is.na(ENTREZID)]
genelist<-unique(res_dt1[!is.na(dmcscore.gene)][order(-dmcscore.gene)],by="ENTREZID")$dmcscore.gene
names(genelist)<-unique(res_dt1[!is.na(dmcscore.gene)][order(-dmcscore.gene)],by="ENTREZID")$ENTREZID

res_gsea_go<- gseGO(geneList     =genelist , 
                    ont="BP",
                    exponent = 1,
                        minGSSize    = 10,maxGSSize = 500,
                        pvalueCutoff = 1,
                        eps = 0,scoreType="pos",
                        OrgDb = org.Hs.eg.db)
res_go_dt<-data.table(as.data.frame(res_gsea_go))
res_go_dt[p.adjust<0.05]#0
