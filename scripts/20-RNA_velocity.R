out<-"outputs/20-RNA_velocity"
dir.create(out)

source("scripts/utils/new_utils.R")

#infer RNA velocity using scVelo to 1) find if LGA HSC have a differentiation bias ?
#2) study which genes (IEGs ?) are freshly transcribed (++ premRNA/ mRNA) 
#and 3) if there is transcription differences between LGA and control

#need first install python

renv::use_python()

library(reticulate)
# Warning message:
# In use_condaenv(condaenv = condaenv, conda = miniconda_conda(),  :
#   multiple Conda environments found; the first-listed will be chosen.
#           name                                                                          python
# 1 r-reticulate /disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin/python
# 2 r-reticulate          /home/apelletier/.local/share/r-miniconda/envs/r-reticulate/bin/python

reticulate::py_install(packages = c("numpy" ,"scipy", "cython" ,"numba" ,"matplotlib", "scikit-learn", "h5py", "click"))
reticulate::py_install(packages ="umap-learn")
reticulate::py_install(packages ="pysam",pip=T)
reticulate::py_install(packages = c("velocyto"),pip = T)
renv::snapshot()

#need also samtools
system("../singlecell/tools/samtools-1.12/samtools")


#run velocyto counting pipeline
  #first, need dl genome repeat sequence to mask here : https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf
#run 20A- [to get from singlecell project scripts]

renv::install("bioc::LoomExperiment")

#then use scVelo (https://scvelo.readthedocs.io/getting_started.html) in python to generate the RNA velocity matrix
# calculate and formulate the anndata
#need install scVelo
reticulate::py_install(packages ="scvelo",pip=T)
#need downgrade numba to 0.52
reticulate::py_install(packages ="numba==0.52")
reticulate::py_run_string("
import scvelo as scv
scv.logging.print_version()
                          ") #Running scvelo 0.2.3 (python 3.6.13) on 2021-04-26 12:29.


out<-"outputs/09-Velocity/sc_velo_outputs"
dir.create(out)

