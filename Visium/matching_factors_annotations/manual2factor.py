import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc
import os
import pickle
import matplotlib as mpl
import warnings
warnings.filterwarnings("ignore",category=mpl.cbook.mplDeprecation)

from fastprogress.fastprogress import progress_bar

from scipy import stats
import statsmodels.stats.multitest as smt
from math import sqrt

# save pickle file
def save_pkl(obj, file):
    with open(file, "wb") as tf:
        pickle.dump(obj, tf)

# read pickle file
def read_pkl(file):
    with open(file, "rb") as tf:
        obj = pickle.load(tf)

    return obj

# for calculate effect size, cohen's d    
def cohend(d1, d2):
    # calculate the size of samples
    n1, n2 = len(d1), len(d2)
    # calculate the variance of the samples
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    # calculate the pooled standard deviation
    s = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1, u2 = np.mean(d1), np.mean(d2)
    # calculate the effect size
    return (u1 - u2) / s


# for matching manual structural annotations and unbiased factors from cell2loc-NMF
def manual2factor(
    adata,
    path_to_cell2locNMF,
    manual_annotation_column,
    mannual_rois='all',
    n_permute=1000,
    adjp_thresh = 0.01,
    n_factor_min = 5,
    n_factor_max = 15,
    plot_heatmap=True,
):

    if mannual_rois=='all':
        mannual_ROIs = set(adata.obs[manual_annotation_column])
    else:
        if len(set(mannual_rois)-set(adata.obs[manual_annotation_column])) > 0: # checking whether mannual_rois are all in the column manual_annotation_column
            raise KeyError(f'{mannual_rois}: invalid <mannual_rois>')
        else:
            mannual_ROIs = mannual_rois
    
    ad_obs = adata.obs.copy()

    # add annotation as boolean data type 
    df = ad_obs[[manual_annotation_column]].reset_index()
    df['values']=int(1)
    df = df.pivot(index='spot_id', columns=manual_annotation_column, values='values')
    df.fillna(int(0),inplace=True)
    ad_obs = pd.concat([ad_obs,df.reindex(ad_obs.index)],axis=1)

    # splot number for each ROI
    spotn_dict={}
    for ROI in mannual_ROIs:
        spotn_dict[ROI]=sum(ad_obs[manual_annotation_column]==ROI)
    

    for ii, nfact in enumerate(np.arange(n_factor_min,n_factor_max)):
        print(f'matching n_fact{nfact}')
        
        fname = os.listdir(path_to_cell2locNMF)[0]
        mod_sk = read_pkl(f'{path_to_cell2locNMF}/{fname}/models/model_n_fact{nfact}.p')['mod']

        # rename factor names
        mod_sk.cell_type_fractions.columns = [x.replace('mean_cell_type_factors','') for x in mod_sk.cell_type_fractions.columns]
        mod_sk.location_factors_df.columns = [x.replace('mean_nUMI_factors','') for x in mod_sk.location_factors_df.columns]

        # add factor values in adata
        ad_obs[f'n_fact{nfact}|'+mod_sk.location_factors_df.columns]=mod_sk.location_factors_df.reindex(ad_obs.index)

        for i, ROI in enumerate(mannual_ROIs):

            rois_obs = ad_obs[ad_obs[manual_annotation_column]==ROI]
            nonrois_obs = ad_obs[ad_obs[manual_annotation_column]!=ROI]

            cd ={}
            for fact in f'n_fact{nfact}|'+mod_sk.location_factors_df.columns:

                # effect size of 'mean_nUMI_factors' between zero and nonzero
                cd[fact] = cohend(rois_obs[fact],nonrois_obs[fact])

            # convert to dataframe
            cd_df_each = pd.DataFrame.from_dict(cd,orient='index').rename(columns={0:ROI})

            # concatate
            if i ==0:
                cd_df = cd_df_each.copy()
            else:
                cd_df = pd.concat([cd_df,cd_df_each], axis=1)

            del cd_df_each

        # permutation test
        fact_list = f'n_fact{nfact}|'+mod_sk.location_factors_df.columns
        n_ROI = len(mannual_ROIs)
        
        for p in progress_bar(range(1000)):

            for i, ROI in enumerate(mannual_ROIs):

                spotn_roi=spotn_dict[ROI]

                # cd={}
                for j, factor_name in enumerate(fact_list):

                    vals_roi,vals_nonroi=np.split(np.random.permutation(ad_obs[factor_name]),[spotn_roi])

                    if j==0:
                        cd_1d=cohend(vals_roi,vals_nonroi)
                    else:
                        cd_1d=np.append(cd_1d,cohend(vals_roi,vals_nonroi))

                    del vals_roi,vals_nonroi

                # append
                if i==0:
                    cd_2d=cd_1d.reshape(nfact,1).copy()
                else:
                    cd_2d=np.append(cd_2d,cd_1d.reshape(nfact,1),axis=1)

                del cd_1d,spotn_roi 

            if p==0:
                cd_3d=cd_2d.reshape(nfact,n_ROI,1).copy()
            else:
                cd_3d=np.append(cd_3d,cd_2d.reshape(nfact,n_ROI,1).copy(),axis=2)

            del cd_2d


        cd_df_p=pd.DataFrame(data=np.nan,index=cd_df.index,columns=cd_df.columns)
        for i in range(cd_df.shape[0]):
            for j in range(cd_df.shape[1]):
                if cd_df.iloc[i,j]>0:
                    cd_df_p.iloc[i,j] = np.sum(cd_3d[i,j,:]>cd_df.iloc[i,j])/n_permute
                else:
                    cd_df_p.iloc[i,j] = np.sum(cd_3d[i,j,:]<cd_df.iloc[i,j])/n_permute

        cd_df_adjp=cd_df_p*cd_df.shape[0]*cd_df.shape[1]
        cd_df_adjp[cd_df_adjp>1]=1

        cd_df = cd_df[cd_df_adjp<adjp_thresh]


        # concatenate all cd_df
        if ii==0:
            cd_df_all=cd_df.copy()
        else:
            cd_df_all=pd.concat([cd_df_all,cd_df])

        # plot
        if plot_heatmap:
            with mpl.rc_context({'font.size': 10,
                                'figure.figsize': [6,4]}):

                # plt.figure(figsize = (3,5))
                cmap = mpl.cm.get_cmap("seismic")
                cmap.set_bad(color='darkgrey')
                mask = cd_df.isnull()
                sns.heatmap(cd_df, cmap=cmap, mask=mask, center=0, cbar_kws={'label': "cohens'd"})
                plt.title(f'n_fact: {nfact}')
                plt.show()

        del cd_df, fact_list, mod_sk
        
    # add to adata.obs 
    adata.obs[[x for x in ad_obs.columns if 'n_fact'==x[0:6]]]=ad_obs[[x for x in ad_obs.columns if 'n_fact'==x[0:6]]].reindex(adata.obs_names)
        
    return adata, cd_df_all


# to find good factors from factor-by-manualannotation dataframe which has cohens d value as above function

def find_good_factors(
    factor_manual_df,
    cohensd_thresh=0.8,
    ratio_to_top=0.3,
    n_factor_min = 5,
    n_factor_max = 15,
    n_dividing=2
):
    
    if n_dividing < 2:
        raise KeyError(f'{n_dividing}: invalid <n_dividing>, required to be 2 or larger')
        
    # select structures which has a given value of max(cohens d)
    DF = factor_manual_df[factor_manual_df.columns[np.max(factor_manual_df)>cohensd_thresh]]

    # identify factors which has a given ration of cohens d value to the value of the top factor (best matched)
    DF = DF.sub(np.max(DF)*ratio_to_top, axis='columns')>0

    for i, nfact in enumerate(np.arange(n_factor_min,n_factor_max)):
        df=pd.DataFrame(data=np.sum(DF.loc[[x for x in DF.index if f'n_fact{nfact}' in x]]),
                        index=DF.columns,
                         columns=[f'n_fact{nfact}']
                        )
        df=df.T

        if i==0:
            n_goodfactor_df=df.copy()
        else:
            n_goodfactor_df=pd.concat([n_goodfactor_df,df])
        del df

    # good factors for each ROI
    goodfactors_dict={}
    for ROI in DF.columns:
        if np.min(n_goodfactor_df[ROI])>n_dividing:
            f = n_goodfactor_df.index[n_goodfactor_df[ROI]==np.min(n_goodfactor_df[ROI])][0]
        elif np.max(n_goodfactor_df[ROI])>(n_dividing-1):
            f = n_goodfactor_df.index[n_goodfactor_df[ROI]==n_dividing][0]
        elif np.max(n_goodfactor_df[ROI])==1:
            f = factor_manual_df.index[factor_manual_df[DF][ROI]==np.max(factor_manual_df[DF][ROI])][0]
        else: # np.max(n_goodfactor_df[ROI])<n_dividing but not 1
            f = n_goodfactor_df.index[n_goodfactor_df[ROI]==np.max(n_goodfactor_df[ROI])][0]

        s = DF.loc[[x for x in DF.index if f in x],ROI]
        goodfactors_dict[ROI] = list(s.index[s])
        del f,s

    return n_goodfactor_df, goodfactors_dict