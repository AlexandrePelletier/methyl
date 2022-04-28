import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
#%load_ext rpy2.ipython


cbp0 = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp0_ctrl/CBP547_possorted_genome_bam_G59L8.loom")
cbp6a = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp6a/possorted_genome_bam_8UU1X.loom")
cbp6b = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp6b/possorted_genome_bam_VUJYA.loom")
cbp6c = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp6c/possorted_genome_bam_KIPKM.loom")
cbp7a = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp7a/possorted_genome_bam_LBOVY.loom")
cbp7b = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp7b/possorted_genome_bam_1M4YO.loom")
cbp7c = anndata.read_loom("../singlecell/outputs/09-Velocity/velocyto_counts/cbp7c/possorted_genome_bam_PYBII.loom")

#filter for cells that pass Seurat QC
#need homogene cell_id first
cbp0.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp0") for item in cbp0.obs.index])
cbp6a.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp6a") for item in cbp6a.obs.index])
cbp6b.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp6b") for item in cbp6b.obs.index])
cbp6c.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp6c") for item in cbp6c.obs.index])

cbp7a.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp7a") for item in cbp7a.obs.index])
cbp7b.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp7b") for item in cbp7b.obs.index])
cbp7c.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp7c") for item in cbp7c.obs.index])

mtd=pd.read_csv("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtd=mtd.loc[mtd['group'] == "ctrl"]
cbp0 = cbp0[np.isin(cbp0.obs.index,mtd["cell_id"])]
cbp6a = cbp6a[np.isin(cbp6a.obs.index,mtd["cell_id"])]
cbp6b = cbp6b[np.isin(cbp6b.obs.index,mtd["cell_id"])]
cbp6c = cbp6c[np.isin(cbp6c.obs.index,mtd["cell_id"])]

cbp7a = cbp7a[np.isin(cbp7a.obs.index,mtd["cell_id"])]
cbp7b = cbp7b[np.isin(cbp7b.obs.index,mtd["cell_id"])]
cbp7c = cbp7c[np.isin(cbp7c.obs.index,mtd["cell_id"])]

cbp0.var_names_make_unique()
cbp6a.var_names_make_unique()
cbp6b.var_names_make_unique()
cbp6c.var_names_make_unique()

cbp7a.var_names_make_unique()
cbp7b.var_names_make_unique()
cbp7c.var_names_make_unique()

cbps = cbp0.concatenate(cbp6a, cbp6b, cbp6c,cbp7a,cbp7b,cbp7c,index_unique=None)

cbps.obs.index

umap_cord = pd.read_csv("outputs/06-integr_singlecell_cbps/umap_cbps.csv")

cbps_index=pd.DataFrame(cbps.obs.index).rename(columns={0:'cell_id'})
umap_ordered=cbps_index.merge(umap_cord, on = "cell_id")
umap_ordered=umap_ordered.iloc[:,2:]
cbps.obsm['X_umap'] = umap_ordered.values

scv.pp.filter_and_normalize(cbps)
scv.pp.moments(cbps)



scv.pp.moments(cbps,n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(cbps)


scv.tl.velocity(cbps, mode = "dynamical")
scv.tl.velocity_graph(cbps)

cbps.write("outputs/20-RNA_velocity/cbps_ctrl_non_hto/cbps_ctrl_non_hto_dynamical_velo.h5ad")

scv.pl.velocity_embedding(cbps, basis = 'umap',save="outputs/20-RNA_velocity/cbps_ctrl_non_hto/umap_dynamical_velocity.pdf")

scv.pl.velocity_embedding_grid(cbps, basis='umap',save="outputs/20-RNA_velocity/cbps_ctrl_non_hto/umap_dynamical_grid_velocity.pdf")

scv.pl.velocity_embedding_stream(cbps, basis='umap',save="outputs/20-RNA_velocity/cbps_ctrl_non_hto/umap_dynamical_stream_velocity.pdf")

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

velo.to_csv("outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_velocity_matrix.csv")

#speed/ rate of differentiation
scv.tl.velocity_confidence(cbps)

#latent time (~ pseudotime)
scv.tl.latent_time(cbps)
scv.pl.scatter(cbps, color='latent_time', color_map='gnuplot', size=80,save="outputs/20-RNA_velocity/cbps_ctrl_non_hto/umap_latent_time_dynamical_velocity.pdf")
cbps.write("outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_velo.h5ad")


#save with metadata
cbps.obs.to_csv("outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_velocity_metadata.csv")

#save velocity_graph (correl cell trans and velocity vec)

velog=pd.DataFrame.sparse.from_spmatrix(cbps.uns['velocity_graph'],columns=cbps.obs.index) 
velog.shape

velog["cell_id"]=cbps.obs.index.values

velog.to_csv("outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_velocity_graph_matrix.csv")

#transition proba matrix 
trans=scv.utils.get_transition_matrix(cbps)

trans=pd.DataFrame.sparse.from_spmatrix(trans,columns=cbps.obs.index)  
trans.head()
trans.shape

trans["cell_id"]=cbps.obs.index.values
trans.to_csv("outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_transition_matrix.csv")


#save velocity umap coord 
velo_u=pd.DataFrame(cbps.obsm['velocity_umap'],columns=['humap_1','humap_2']) #add genes name 
velo_u.shape
velo_u["cell_id"]=cbps.obs.index.values

velo_u.to_csv("outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_velocity_umap_coord.csv")


#interpret important genes
scv.pl.velocity(cbps, ['EGR1',  'KLF2', 'SOCS3', 'JUNB'], ncols=2,save="outputs/20-RNA_velocity/cbps_ctrl_non_hto/dynamical_velocity_important_genes.pdf")

