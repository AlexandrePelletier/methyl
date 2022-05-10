import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
#%load_ext rpy2.ipython

cbp2 = anndata.read_loom("outputs/20-RNA_velocity/velocyto_counts/cbp2/possorted_genome_bam_97NV3.loom")
cbp3 = anndata.read_loom("outputs/20-RNA_velocity/velocyto_counts/cbp3/possorted_genome_bam_CF6IK.loom")
cbp4 = anndata.read_loom("outputs/20-RNA_velocity/velocyto_counts/cbp4/possorted_genome_bam_QIW9J.loom")
cbp8 = anndata.read_loom("outputs/20-RNA_velocity/velocyto_counts/cbp8/possorted_genome_bam_B7ECE.loom")

#filter for cells that pass Seurat QC
#need homogene cell_id first
cbp2.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp2") for item in cbp2.obs.index])
cbp3.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp3") for item in cbp3.obs.index])
cbp4.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp4") for item in cbp4.obs.index])
cbp8.obs.index=pd.Index([item.split(":")[1].replace("x","_cbp8") for item in cbp8.obs.index])

mtd=pd.read_csv("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
cbp2 = cbp2[pd.Index(set(cbp2.obs.index)&set(mtd.cell_id.values))]

cbp3 = cbp3[np.isin(cbp3.obs.index,mtd["cell_id"])]
cbp4 = cbp4[np.isin(cbp4.obs.index,mtd["cell_id"])]
cbp8 = cbp8[np.isin(cbp8.obs.index,mtd["cell_id"])]

cbp2.var_names_make_unique()
cbp3.var_names_make_unique()
cbp4.var_names_make_unique()
cbp8.var_names_make_unique()

cbps = cbp2.concatenate(cbp3, cbp4, cbp8,index_unique=None)

cbps.obs.index

umap_cord = pd.read_csv("outputs/06-integr_singlecell_cbps/umap_cbps.csv")

cbps_index=pd.DataFrame(cbps.obs.index).rename(columns={0:'cell_id'})
umap_ordered=cbps_index.merge(umap_cord, on = "cell_id")
umap_ordered=umap_ordered.iloc[:,2:]
cbps.obsm['X_umap'] = umap_ordered.values

cbps.write("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad")

scv.pp.filter_and_normalize(cbps)
scv.pp.moments(cbps)
scv.tl.velocity(cbps, mode = "stochastic")
scv.tl.velocity_graph(cbps)
cbps.write("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad")

scv.pl.velocity_embedding(cbps, basis = 'umap',save="outputs/20-RNA_velocity/umap_velocity.pdf")

scv.pl.velocity_embedding_grid(cbps, basis='umap',save="outputs/20-RNA_velocity/umap_grid_velocity.pdf")

scv.pl.velocity_embedding_stream(cbps, basis='umap',save="outputs/20-RNA_velocity/umap_stream_velocity.pdf")

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

velo.to_csv("outputs/20-RNA_velocity/cbps_hto_velocity_matrix.csv")

#speed/ rate of differentiation
scv.tl.velocity_confidence(cbps)
cbps.write("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad")

#save with metadata
cbps.obs.to_csv("outputs/20-RNA_velocity/cbps_hto_velocity_metadata.csv")

#save velocity_graph (correl cell trans and velocity vec)
cbps=anndata.read("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad")

velog=pd.DataFrame.sparse.from_spmatrix(cbps.uns['velocity_graph'],columns=cbps.obs.index) 
velog.shape

velog["cell_id"]=cbps.obs.index.values

velog.to_csv("outputs/20-RNA_velocity/cbps_hto_velocity_graph_matrix.csv")

#transition proba matrix 
trans=scv.utils.get_transition_matrix(cbps)

trans=pd.DataFrame.sparse.from_spmatrix(trans,columns=cbps.obs.index)  
trans.head()
trans.shape

trans["cell_id"]=cbps.obs.index.values
trans.to_csv("outputs/20-RNA_velocity/cbps_hto_transition_matrix.csv")


#save velocity umap coord 
velo_u=pd.DataFrame(cbps.obsm['velocity_umap'],columns=['humap_1','humap_2']) #add genes name 
velo_u.shape
velo_u["cell_id"]=cbps.obs.index.values

velo_u.to_csv("outputs/20-RNA_velocity/cbps_hto_velocity_umap_coord.csv")


#interpret important genes
scv.pl.velocity(cbps, ['EGR1',  'KLF2', 'SOCS3', 'JUNB'], ncols=2,save="outputs/20-RNA_velocity/velocity_important_genes.pdf")

