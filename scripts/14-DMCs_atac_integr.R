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
res_cpgs<-merge(res_meth[,.(cpg_id,logFC,AveExpr,P.Value,adj.P.Val)],unique(cpgs[,.(cpg_id,chr,pos)]))
res_cpgs[is.na(pos)] #ok
cpgs_in_ocrs<-bed_inter(res_cpgs[,start:=pos][,end:=pos+1][,.(chr,start,end,cpg_id)][order(chr,start)],
          ocrs[,.(chr,start,end,peaks)],
          select = c(4,1,2,8,6,7),col.names = c("cpg_id","chr","pos","peaks","start","end"))

fwrite(cpgs_in_ocrs,fp(out,"cpgs_in_bulk_OCRs.csv.gz"))
res_cpgs_ocrs<-merge(res_meth[,.(cpg_id,logFC,AveExpr,P.Value,adj.P.Val)],cpgs_in_ocrs)

res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25] #369 
fwrite(res_cpgs_ocrs[P.Value<0.001&abs(logFC)>25],fp(out,"DMCs_in_bulk_OCRs.csv.gz"))

res_meth[P.Value<0.001&abs(logFC)>25] #369 / 4815 DMCs
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
)

enriched.motifs<-data.table(enriched.motifs)
fwrite(enriched.motifs,fp(out,"enriched_motifs_CpGs_overlaping_OCRs_with_DMCs_vs_without.csv"))

enriched.motifs[pvalue<0.0001]$motif.name
enriched.motifs[,motif.name:=factor(motif.name,levels=enriched.motifs[order(pvalue)]$motif.name)]

ggplot(enriched.motifs[order(pvalue)][1:20])+
  geom_point(aes(x=motif.name,size=fold.enrichment,y=observed,col=-log10(pvalue)))

MotifPlot(
  object = cbps,
  motifs = head(rownames(enriched.motifs))
)



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
#need macs2

renv::use_python()
#reticulate::install_miniconda()

reticulate::use_miniconda()

reticulate::py_install(packages ="MACS2")

peaks <- CallPeaks(
  object = atacs,
  group.by = "predicted.id",
  macs2.path = ""
)


#II) methylOCR regul EGR1 network ?

