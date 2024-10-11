# -*- coding: utf-8 -*-
"""
@File    :   _rank_psis_groups.py
@Time    :   2024/09/01
@Author  :   Dawn
@Version :   1.0
@Desc    :   DPSI for scCyclone
"""


import numpy as np
import pandas as pd
import anndata as ad
from joblib import Parallel, delayed
from typing import Union


from . import _utils
     

def _compute_dpsi(e, data_a, data_b):
    """
    Compute dpsi values for a given event.

    Parameters
    ----------
    e : str
        Event identifier.
    data_a : pd.DataFrame
        Data frame for group A.
    data_b : pd.DataFrame
        Data frame for group B.

    Returns
    -------
    dict
        Dictionary containing the event, shuffled dpsi values, and observed dpsi value.
    """
    sub_a = data_a[[e]].dropna()
    sub_b = data_b[[e]].dropna()
    
    min_shape = min(sub_a.shape[0], sub_b.shape[0])
    dpsi_observed = np.round((sub_b.median() - sub_a.median()).values[0], 3)
    
    if min_shape >= 50:
        sub_data = pd.concat([sub_a, sub_b])
        label_list = np.array(["ref"] * sub_a.shape[0] + ["target"] * sub_b.shape[0])
        
        dpsi_list = []
        for _ in range(100):
            np.random.shuffle(label_list)
            shuffle_a = sub_data[label_list == "ref"]
            shuffle_b = sub_data[label_list == "target"]
            
            dpsi = np.round((shuffle_b.mean() - shuffle_a.mean()).values[0], 3)
            dpsi_list.append(dpsi)
    else:
        dpsi_list = [None]
    
    return {
        "event": e,
        "dpsi_shuffle": dpsi_list,
        "dpsi_observed": dpsi_observed
    }
    
    

def _compute_pvalue(row):
    """
    Compute p-value based on observed and shuffled dpsi values.

    Parameters
    ----------
    row : pandas.Series
        A row containing 'dpsi_shuffle' and 'dpsi_observed' values.

    Returns
    -------
    float or None
        Computed p-value.
    """
    dpsi_shuffle = row['dpsi_shuffle']
    dpsi_observed = row['dpsi_observed']
    
    if len(set(dpsi_shuffle)) > 1:
        count_greater = np.sum(np.where(np.array(dpsi_shuffle) > dpsi_observed, 1, 0))
        pvalue = count_greater / len(dpsi_shuffle)
    else:
        pvalue = 1

    return pvalue



def rank_psis_groups(
    adata: ad.AnnData,
    groupby: str,
    groups: Union[str, list]="all",
    reference: str = "rest",
    key_added: Union[str, None] = None):
    """
    Rank psi for characterizing groups.

    Expects logarithmized data.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The key of the observations grouping to consider.
    groups
        Subset of groups, e.g. ['g1', 'g2', 'g3'], to which comparison shall be restricted, or 'all' (default), for all groups.
    reference
        If 'rest', compare each group to the union of the rest of the group. If a group identifier, compare with respect to this group.
    key_added
        The key in adata.uns information is saved to.

    Returns
    -------
    ad.AnnData
        Annotated data matrix with rank information stored in adata.uns[key_added].
    """

    groups_order = _utils.check_groups(adata,groupby,groups,reference)
    print(groups_order)


    # Initialize adata.uns[key_added]
    
    if key_added is None:
        key_added = "rank_psis_groups"
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = {"groupby": groupby, "reference": reference}

    # Initialize empty dictionaries for storing results
    data_event_dict = {}
    data_dpsi_observed_dict = {}
    data_pval_dict = {}
    data_pval_adj_dict = {}
    

    # Iterate over groups
    for i in groups_order:
        if i != reference:
            print("Group {} start!".format(i))
            event_data_list = []
            adata_target = adata[adata.obs[groupby] == i]
            adata_ref = adata[adata.obs[groupby] != i] if reference == "rest" else adata[adata.obs[groupby].isin([reference])]
            data_target = adata_target.to_df()
            data_ref = adata_ref.to_df()

            # Parallel computation of dpsi
            print("Compute dpsi...")
            rs = Parallel(n_jobs=-1)(delayed(_compute_dpsi)(e, data_ref, data_target) for e in adata.var.index)
            event_data_list.extend(rs)
            result = pd.DataFrame(event_data_list)

            # Compute p-values and sort results
            print("Compute pvalue...")
            result['pvals'] = result.apply(_compute_pvalue, axis=1)
            result = result.sort_values(by=['dpsi_observed', 'pvals'], ascending=[False, True])

            # Store results in dictionaries
            data_event_dict[i] = result['event'].to_list()
            data_dpsi_observed_dict[i] = result['dpsi_observed'].to_list()
            data_pval_dict[i] = result['pvals'].to_list()
            data_pval_adj_dict[i] = _utils.compute_pvalue_bonferroni(result['pvals'].to_list())

            print("Group {} complete!".format(i))
            print("-----------------------------------------")
        
    # Convert dictionaries to structured arrays
    name_data = pd.DataFrame(data_event_dict).to_records(index=False)
    dpsi_data = pd.DataFrame(data_dpsi_observed_dict).to_records(index=False)
    pval_data = pd.DataFrame(data_pval_dict).to_records(index=False)
    pval_adj_data = pd.DataFrame(data_pval_adj_dict).to_records(index=False)

    # Store results in adata.uns
    adata.uns[key_added]['names'] = name_data
    adata.uns[key_added]['dpsi'] = dpsi_data
    adata.uns[key_added]['pvals'] = pval_data
    adata.uns[key_added]['pvals_adj'] = pval_adj_data

    return adata



        
    


