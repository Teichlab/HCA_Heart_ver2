#!/usr/bin/env python
# coding: utf-8

# ref: https://docs.scvi-tools.org/en/stable/user_guide/notebooks/MultiVI_tutorial.html

# ## Import modules

# In[ ]:


import scvi
import anndata
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

scvi.settings.seed = 420
scvi.settings.num_threads = 24


# In[ ]:


import session_info
session_info.show()


# ## Read in data

# In[ ]:


# adata_mvi=sc.read_h5ad('/nfs/team205/kk18/data/6region_v2/MultiVI/adata_mvi_downsized.h5ad')

path_adata = '/nfs/team205/heart/anndata_objects/'
directory = path_adata + 'MultiVI'
adata_mvi=sc.read_h5ad(directory + '/adata_mvi.h5ad')

print(adata_mvi.X.data[:10])
adata_mvi


# In[ ]:


adata_mvi.obs['modality'].value_counts()


# In[ ]:


adata_mvi.obs.donor_cellnuc.isna().sum()


# In[ ]:


adata_mvi.var.modality.isna().sum()


# In[ ]:


adata_mvi.var['modality'].value_counts()


# ## Setup and Training MultiVI

# In[ ]:


scvi.data.setup_anndata(adata_mvi, batch_key='modality', categorical_covariate_keys=['donor_cellnuc'])
# scvi.data.setup_anndata(adata_mvi, batch_key='modality', categorical_covariate_keys=['Donor'])\


# In[ ]:


scvi.data.view_anndata_setup(adata_mvi)


# In[ ]:


# When creating the object, we need to specify how many of the features are genes, and how many are genomic regions. 
# This is so MultiVI can determine the exact architecture for each modality.
mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes=(adata_mvi.var['modality']=='Gene Expression').sum(),
    n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
)


# In[ ]:


mvi


# In[ ]:


mvi.train()


# In[ ]:


# mvi.save("/nfs/team205/kk18/data/6region_v2/MultiVI/model/trained_multivi_downsized_20210922")
mvi.save(directory + "/trained_multivi_full_2021_10_06")


# In[ ]:


adata_mvi.obsm["MultiVI_latent"] = mvi.get_latent_representation()


# In[ ]:


# adata_mvi.write('/nfs/team205/kk18/data/6region_v2/MultiVI/adata_post-multivi_downsized.h5ad')
adata_mvi.write(directory + '/adata_post-multivi_full_2021_10_06.h5ad')


# In[ ]:


sc.pp.neighbors(adata_mvi, use_rep="MultiVI_latent")
sc.tl.umap(adata_mvi, min_dist=0.2)
# sc.pl.umap(adata_mvi, color=['modality','cell_states'])
sc.pl.umap(adata_mvi, color=['modality','region'])

