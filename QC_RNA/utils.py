from typing import Dict
import numpy as np
import tables
import scipy.sparse as sp

from statsmodels.stats.multitest import multipletests
import scrublet as scr
import scanpy as sc
import pandas as pd
import numpy as np
import scipy


# From Kazumasa
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.

    """

    try:
        import anndata
    except ImportError:
        raise ImportError('The anndata package must be installed to use the '
                          'function anndata_from_h5()')

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the count matrix.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)})
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata


# From Krzysztof:
# some functions that Ni uses in scanpy scripts to run scrublet
# which in turn are inspired by my original notebook on the matter
# (extracted from scanpy_scripts 0.2.10 to get around scanpy version incompatibility)
def test_outlier(x, upper_mad_only=True):
    med = np.median(x)
    if upper_mad_only:
        mad = np.median(x[x>med] - med) * 1.4826
    else:
        mad = np.median(np.abs(x - med)) * 1.4826
    pvals = 1 - scipy.stats.norm.cdf(x, loc=med, scale=mad)
    bh_pvals = multipletests(pvals, method='fdr_bh')[1]
    return pvals, bh_pvals

def run_scrublet(adata, resolution_function=None):
    
    old_verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1
    if resolution_function is None:
        resolution_function = lambda x: np.maximum(np.maximum(np.log10(x)-1, 0)**2, 0.1)
    scrub = scr.Scrublet(adata.X)
    #this has the potential to brick for poor quality data
    #if so, abort it and everything downstream
    try:
        ds, pd = scrub.scrub_doublets(verbose=False)
    except:
        return
    adata.obs['scrublet_score'] = ds

    adata_copy = adata.copy()
    sc.pp.filter_genes(adata_copy, min_cells=3)
    sc.pp.normalize_total(adata_copy, target_sum=1e4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(adata_copy, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True)
    sc.pp.scale(adata_copy, zero_center=False)
    sc.pp.pca(adata_copy, svd_solver='arpack', zero_center=False)
    sc.pp.neighbors(adata_copy, n_pcs=30)
    sc.tl.umap(adata_copy)
    sc.tl.leiden(adata_copy, resolution=1)
    for clst in np.unique(adata_copy.obs['leiden']):
        clst_size = sum(adata_copy.obs['leiden'] == clst)
        sc.tl.leiden(adata_copy, restrict_to=('leiden', [clst]), 
                     resolution=resolution_function(clst_size), key_added='leiden_R')
        adata_copy.obs['leiden'] = adata_copy.obs['leiden_R']
    clst_meds = []
    for clst in np.unique(adata_copy.obs['leiden']):
        k = adata_copy.obs['leiden'] == clst
        clst_med = np.median(adata_copy.obs.loc[k, 'scrublet_score'])
        adata_copy.obs.loc[k, 'cluster_scrublet_score'] = clst_med
        clst_meds.append(clst_med)
    clst_meds = np.array(clst_meds)
    pvals, bh_pvals = test_outlier(clst_meds)
    for i, clst in enumerate(np.unique(adata_copy.obs['leiden'])):
        k = adata_copy.obs['leiden'] == clst
        adata_copy.obs.loc[k, 'pval'] = pvals[i]
        adata_copy.obs.loc[k, 'bh_pval'] = bh_pvals[i]
    sc.settings.verbosity = old_verbosity
    #need to also export the clustering, for soupx purposes
    adata.obs['scrublet_leiden'] = adata_copy.obs['leiden']
    adata.obs['scrublet_score'] = adata_copy.obs['scrublet_score']
    adata.obs['cluster_scrublet_score'] = adata_copy.obs['cluster_scrublet_score']
    adata.obs['doublet_pval'] = adata_copy.obs['pval']
    adata.obs['doublet_bh_pval'] = adata_copy.obs['bh_pval']
    del adata_copy
    
