
import scvelo as scv
import scanpy as sc
import pandas as pd
import math
import anndata
adata = sc.read_h5ad("outputs/20-RNA_velocity/cbp2.h5ad")
adata.var.index.values
velo=pd.DataFrame(adata.layers['velocity'],columns=adata.var.index) #add genes name 
velo.shape
#get cell id
velo["cell_id"]=adata.obs.index.values


velo.to_csv("outputs/20-RNA_velocity/cbp2_velocity_vec_matrix.csv")
