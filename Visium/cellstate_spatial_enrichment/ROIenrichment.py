import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from scipy.stats import chi2_contingency
import statsmodels.stats.multitest as smt


# Function to calculate odds ratio from 2-by-2 contingency table: (ROI, non-ROI)-by-(target cell, other cells)
# Also computes the chi-square statistic and p-value with `scipy.stats.chi2_contingency` 

# @param adata: cell2location output anndata with slides of interest
# @param roi: Region of interest listed in annotation csv files
# @param target: Target cell type
# @param ref: Reference cell types including target cell. 'all' or list of reference cell types you want to test (eg. immune cells)
# @param category_name: The column name in adata.obs which includes ROIs
# @param prop_header: The header of the culumn in adata.obs which has cell proportion values

def ROI_enrichment(adata, roi, target, ref='all', category_name='Annotation', prop_header='prop_'):
    
    # cells other than target cell type
    if ref=='all':
        allcells = [x.replace(prop_header,'') for x in adata.obs.columns if prop_header in x]
        others = list(set(allcells) - set(target))
    else:
        others = list(set(ref) - set(target))
    
    # extract proportion data to a dataframe
    df=adata.obs[[x for x in adata.obs.columns if prop_header in x]]
    df.columns=[x.replace(prop_header,'') for x in df.columns]
    
    # separate data into roi and non-roi spots
    df_roi = df.loc[adata.obs[adata.obs[category_name]==roi].index]
    df_nonroi = df.loc[adata.obs[adata.obs[category_name]!=roi].index]
    
    # calculate total proportion value per condition
    # contingency table: 2 annotations (ROI, non-ROI) x 2 cell categories (target cell, other cells) 
    roi_target = np.sum(df_roi[[target]])[0]
    nonroi_target = np.sum(df_nonroi[[target]])[0]
    roi_others = np.sum(np.sum(df_roi[others], axis=1))
    nonroi_others = np.sum(np.sum(df_nonroi[others], axis=1))
    
    # calculate logOR, SE, and error bar for plotting
    logOR = np.log((roi_target/nonroi_target)/(roi_others/nonroi_others))
    SE = np.sqrt((1/roi_target)+(1/nonroi_target)+(1/roi_others)+(1/nonroi_others))
    ebar = 1.96*SE
    
    # calculate chi-squared value and p-value
    table = np.array([[roi_target,nonroi_target],
                      [roi_others,nonroi_others]])
    x2, p, dof, expected = chi2_contingency(table)
    
    return dict({'logOR':logOR, 'ebar':ebar, 'x2':x2, 'p':p})


# Calculate enrichment for all cell types in the reference with adjusted p-value
# @param method: Method used for adjustment of pvalues.
## using `statsmodels.stats.multitest.multipletests`(https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html)

def ROI_enrichment_allcomparison(adata, roi, ref='all', category_name='Annotation', prop_header='prop_', method='fdr_bh'):
    if ref=='all':
        allcells = [x.replace(prop_header,'') for x in adata.obs.columns if prop_header in x]
        res = [ROI_enrichment(adata, roi=roi, target=x, ref=allcells, category_name=category_name, prop_header=prop_header) for x in allcells]
        res = pd.DataFrame(res)
        res['cell'] = allcells
    else:
        res = [ROI_enrichment(adata, roi=roi, target=x, ref=ref, category_name=category_name, prop_header=prop_header) for x in ref]
        res = pd.DataFrame(res)
        res['cell'] = ref
    
    res['p_adj']=smt.multipletests(res['p'], method=method, is_sorted=False, returnsorted=False)[1]
    
    return res


# Plot ROI_enrichment (regions x cell states) 
def plot_ROIenrichment(adata,
                       roi,
                       category_name='Annotation',
                       prop_header='prop_',
                       ref='all',
                       method='fdr_bh',
                       figure_dpi=200,
                       dendrogram_ratio=(.3),
                       cbar_pos=(0.983, .45, .005, .2),
                       figsize=(30,3.5),
                       title_fontsize=24,
                       title_x=0.6,
                       title_y=1.05,
                       select_significant=False,
                       p_adj_thresh=0.05,
                       dendrogram=True,
                       regions_order=['SAN','RA','LA','RV','LV','SP','AX']
                      ):
    
    # select regions which 
    ct=pd.crosstab(adata.obs['region'],adata.obs[category_name])
    regions=list(ct.index[ct[roi]>0])
    
    for reg in regions:
        # subset region anndata
        adata_reg=adata[adata.obs['region']==reg].copy()
        # remove cell states not exist in the region
        ## drop nan
        adata_reg.obs.dropna(axis='columns',how='all',inplace=True)
        ## drop zero proportion
        rm_list = [x for x in np.sum(adata_reg.obs)[np.sum(adata_reg.obs)==0].index if prop_header in x]
        adata_reg.obs.drop(rm_list, axis=1, inplace=True)
        # ROI enrichment
        res = ROI_enrichment_allcomparison(adata_reg,roi=roi,ref=ref,method=method,category_name=category_name)
        res['region']=reg
        if reg==regions[0]:
            res_all=res.copy()
        else:
            res_all=pd.concat([res_all,res])

    if select_significant:
        # select cell states which is significantly enriched (postively) in at least one region
        res_all['significance_positive']=(res_all['p_adj']<p_adj_thresh)&(res_all['logOR']>0)
        res_grouped=res_all.groupby('cell').max()
        cell_selected=res_grouped[res_grouped['significance_positive']>0].index
        if len(cell_selected)==0:
            raise Exception("no positively enriched cell states")
        
        # filter dataframe
        res_all=res_all[res_all['cell'].isin(cell_selected)]
    else:
        pass
        
    # pivot and fill nan with 0
    df=pd.pivot(res_all,index='region',columns='cell',values='logOR',)
    df.fillna(0,inplace=True)
    
    # masking color
    cmap = mpl.cm.get_cmap("bwr")
    cmap.set_bad(color='darkgrey')
    
    # plotting
    plt.rcParams["figure.dpi"] = figure_dpi
    
    if dendrogram:
        ax = sns.clustermap(df, cmap=cmap, mask=(df==0),
                            linewidths=.75,dendrogram_ratio=dendrogram_ratio,
                           cbar_pos=cbar_pos,
                            figsize=figsize,center=0,
                        cbar_kws={'label': 'logOR'}
                           )
        ax.ax_row_dendrogram.remove()  
    else:
        cell_order=res_all.groupby('cell').sum().sort_values('logOR',ascending=False).index
        df=df.loc[[x for x in regions_order if x in df.index],cell_order]
        ax = sns.clustermap(df, cmap=cmap, mask=(df==0),
                            linewidths=.75,row_cluster=False, col_cluster=False,
                           cbar_pos=cbar_pos,
                            figsize=figsize,center=0,
                        cbar_kws={'label': 'logOR'}
                           )
        
    ax.fig.suptitle(roi,fontsize=title_fontsize, y=title_y, x=title_x) 
    ax.ax_heatmap.yaxis.set_ticks_position("left")
    
    ## remove axis label
    axh=ax.ax_heatmap
    axh.set_ylabel('')
    axh.set_xlabel('')
    
    plt.setp(ax.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.show()
    