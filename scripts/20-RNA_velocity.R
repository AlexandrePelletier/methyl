out<-"outputs/20-RNA_velocity"
dir.create(out)

source("scripts/utils/new_utils.R")

#infer RNA velocity using scVelo to 1) find if LGA HSC have a differentiation bias ?
#2) study which genes (IEGs ?) are freshly transcribed (++ premRNA/ mRNA) 
#and 3) if there is transcription differences between LGA and control

####need first install python####
renv::snapshot()
#rm python files

renv::use_python()

library(reticulate)

reticulate::py_install(packages = c("numpy" ,"scipy", "cython" ,"numba" ,"matplotlib", "scikit-learn", "h5py", "click"))
reticulate::py_install(packages ="umap-learn")
reticulate::py_install(packages ="pysam",pip=T)
reticulate::py_install(packages = c("velocyto"),pip = T)
reticulate::py_install(packages ="scvelo",pip=T)

#need downgrade numba to 0.52
reticulate::py_install(packages ="numba==0.52")
reticulate::py_run_string("
import scvelo as scv
scv.logging.print_version()
                          ") ##Running scvelo 0.2.4 (python 3.7.3) on 2022-04-07 13:23

renv::install("bioc::LoomExperiment")
renv::snapshot()

#need also samtools
system("../singlecell/tools/samtools-1.12/samtools")


#run velocyto counting pipeline
  #first, need dl genome repeat sequence to mask here : https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf
#run 20A- [to get from singlecell project scripts]


#then use scVelo (https://scvelo.readthedocs.io/getting_started.html) in python to generate the RNA velocity matrix
# calculate and formulate the anndata

#test for 1 : CBP2
reticulate::py_run_string("
import scvelo as scv
filename='../singlecell/outputs/09-Velocity/velocyto_counts/cbp2/possorted_genome_bam_97NV3.loom'
out='outputs/20-RNA_velocity'
#load loom file
adata = scv.read(filename, cache=True)

#adata.var_names_make_unique()
#scv.pl.proportions(adata)

# #basic preprocessing
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#recover dynmics
scv.tl.recover_dynamics(adata)

scv.tl.velocity(adata, mode='dynamical')
adata.write(out+'/cbp2.h5ad', compression='gzip')

scv.tl.velocity_graph(adata)

adata.write(out+'/cbp2.h5ad', compression='gzip')

                          ")

#get velocity matrix 
#run 20B
velo_mat<-fread("outputs/20-RNA_velocity/cbp2_velocity_vec_matrix.csv")
sum(!is.na(as.matrix(velo_mat)))

#for all HTO data
#need first merge all velocito object #see https://github.com/basilkhuder/Seurat-to-RNA-Velocity#integrating-loom-file-and-meta-data
#with cells kept in seurat analysis
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")

mtd[,cell_id:=paste(str_extract(bc,"[ATCG]+"),orig.ident,sep="_")]
fwrite(mtd,"outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
#update reticulat because of bug
renv::install("rstudio/reticulate")

#get umap coord
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
umap_coord<-data.table(cbps@reductions$ref.umap@cell.embeddings,keep.rownames = "bc")
umap_coord<-merge(umap_coord,mtd[,.(bc,cell_id)])
fwrite(umap_coord,"outputs/06-integr_singlecell_cbps/umap_cbps.csv")

table(mtd[hto==T]$batch)

#run 20B-scvelo
#get velocity matrix
library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad", dest = "h5seurat", overwrite = TRUE)
cbps_velo <- LoadH5Seurat("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5seurat") #Error: Missing required datasets 'levels' and 'values'

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
head(cbps[[]])

velo<-fread("outputs/20-RNA_velocity/cbps_hto_velocity_matrix.csv")
velo[1:10,.(cell_id)]
velo<-as.matrix(t(data.frame(velo,row.names = "cell_id")))

mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdf<-mtd[cell_id%in%colnames(velo)]
cbps_h<-cbps[,mtdf$bc]
rm(cbps)
cbps_h<-RenameCells(cbps_h,new.names = mtdf$cell_id )

cbps_h[["velocity"]]<-CreateAssayObject(data =velo[2:nrow(velo),] )

mtdvelo<-fread("outputs/20-RNA_velocity/cbps_hto_velocity_metadata.csv")
mtdvelo[,cell_id:=V1]
cbps_h<-AddMetaData(cbps_h,metadata = data.frame(mtdvelo[,-c("V1","batch")],row.names = "cell_id"))
saveRDS(cbps_h,fp(out,"cbps_hto_with_velocity_assay.rds"))

#next : analyze diff bias, diff gene regulation by lineage (IEGs upreg in HSC ?), reduce upregulation IEGs in LGA HSC ?
cbps_h<-readRDS(fp(out,"cbps_hto_with_velocity_assay.rds"))

#diff bias
#diff "intensity" by lineage
mtd<-data.table(cbps_h@meta.data,keep.rownames = "cell_id")
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")
mtdl<-mtd[lineage_hmap%in%lins]
mtdl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdl[,velo_len_avg:=mean(velocity_length),by=.(lineage_hmap,sample)]

ggplot(unique(mtdl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=velocity_length,fill=group))

#Differentiation Analysis : ponderate average of the pseudotime of the predicted cells transition
#need=> transition graph : correl change of expression between 2 cells and predict change based on the RNA velocity
# need get velocity_graph
trans<-fread("outputs/20-RNA_velocity/cbps_hto_transition_matrix.csv")
trans<-as.matrix(t(data.frame(trans,row.names = "cell_id")))

trans[1:10,1:10]
trans<-trans[2:nrow(trans),]

mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdf<-mtd[cell_id%in%colnames(trans)]
#add pseudotime
pseudo<-fread("outputs/12-Pseudotime/metadata_pseudotime.csv")
pseudo[,cell_id:=ps(str_extract(bc,"[ATCG]+"),orig.ident,sep='_')]
mtdf<-merge(mtdf,pseudo[,.(cell_id,pseudotime)])
mtdff<-mtdf[pseudotime!=Inf]
summary(mtdff$pseudotime)
mtdff[,pseudo_bias:=sapply(1:.N,function(i){
  cell<-cell_id[i]
  return(sum(pseudotime*trans[cell,cell_id]))
  })]

ggplot(mtdff)+geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

mtdffl<-mtdff[lineage_hmap%in%lins]
mtdffl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdffl[,avg_pseudobias:=mean(pseudo_bias),by=.(lineage_hmap,sample)]

ggplot(unique(mtdffl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=avg_pseudobias,fill=group))
unique(mtdffl[lineage_hmap=="MPP/LMPP"],by=c("sample","lineage_hmap"))

#[to continue]

