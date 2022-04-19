import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
#%load_ext rpy2.ipython

cbps=anndata.read("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad")

scv.pp.moments(cbps,n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(cbps)


scv.tl.velocity(cbps, mode = "dynamical")
scv.tl.velocity_graph(cbps)

cbps.write("outputs/20-RNA_velocity/cbps_hto_dynamical_velo.h5ad")

scv.pl.velocity_embedding(cbps, basis = 'umap',save="outputs/20-RNA_velocity/umap_dynamical_velocity.pdf")

scv.pl.velocity_embedding_grid(cbps, basis='umap',save="outputs/20-RNA_velocity/umap_dynamical_grid_velocity.pdf")

scv.pl.velocity_embedding_stream(cbps, basis='umap',save="outputs/20-RNA_velocity/umap_dynamical_stream_velocity.pdf")

#save velocity matrix :
#velocities are obtained by modeling transcriptional dynamics of splicing kinetics.
#For each gene, a steady-state-ratio of pre-mature (unspliced) and mature (spliced) mRNA counts is fitted,
#which constitutes a constant transcriptional state. Velocities are then obtained as residuals from this ratio. 
#Positive velocity indicates that a gene is up-regulated, which occurs for cells that show higher abundance of 
#unspliced mRNA for that gene than expected in steady state. 
#Conversely, negative velocity indicates that a gene is down-regulated.

#The computed velocities are stored in adata.layers just like the count matrices.
velo=pd.DataFrame(cbps.layers['velocity'],columns=cbps.var.index) #add genes name 
velo.shape
velo["cell_id"]=cbps.obs.index.values

velo.to_csv("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_matrix.csv")

#speed/ rate of differentiation
scv.tl.velocity_confidence(cbps)

#latent time (~ pseudotime)
scv.tl.latent_time(cbps)
scv.pl.scatter(cbps, color='latent_time', color_map='gnuplot', size=80,save="outputs/20-RNA_velocity/umap_latent_time_dynamical_velocity.pdf")
cbps.write("outputs/20-RNA_velocity/cbps_hto_dynamical_velo.h5ad")


#save with metadata
cbps.obs.to_csv("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_metadata.csv")

#save velocity_graph (correl cell trans and velocity vec)

velog=pd.DataFrame.sparse.from_spmatrix(cbps.uns['velocity_graph'],columns=cbps.obs.index) 
velog.shape

velog["cell_id"]=cbps.obs.index.values

velog.to_csv("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_graph_matrix.csv")

#transition proba matrix 
trans=scv.utils.get_transition_matrix(cbps)

trans=pd.DataFrame.sparse.from_spmatrix(trans,columns=cbps.obs.index)  
trans.head()
trans.shape

trans["cell_id"]=cbps.obs.index.values
trans.to_csv("outputs/20-RNA_velocity/cbps_hto_dynamical_transition_matrix.csv")


#save velocity umap coord 
velo_u=pd.DataFrame(cbps.obsm['velocity_umap'],columns=['humap_1','humap_2']) #add genes name 
velo_u.shape
velo_u["cell_id"]=cbps.obs.index.values

velo_u.to_csv("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_umap_coord.csv")


#interpret important genes
scv.pl.velocity(cbps, ['EGR1',  'KLF2', 'SOCS3', 'JUNB'], ncols=2,save="outputs/20-RNA_velocity/dynamical_velocity_important_genes.pdf")

