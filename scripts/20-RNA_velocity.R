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
out<-"outputs/09-Velocity/sc_velo_outputs"
dir.create(out)

