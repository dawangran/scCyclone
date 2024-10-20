# -*- coding: utf-8 -*-
"""
@File    :   utils.py
@Time    :   2024/09/01
@Author  :   Dawn
@Version :   1.0
@Desc    :   Little tools for scCyclone
"""


import anndata as ad
from typing import Union
from statsmodels.stats.multitest import multipletests



def check_groups(
    adata: ad.AnnData,
    groupby: str,
    groups: Union[str, list] = "all",
    reference: str = "rest",
):
    """
    Check and validate groups based on a specific column in AnnData object.

    Parameters:
    adata (ad.AnnData): AnnData object containing the data.
    groupby (str): Column in adata.obs to group by.
    groups (Union[str, list]): Specific groups to check.
    reference (str): Reference group.

    Returns:
    list: List of validated groups.
    """
    
    # Check if groupby column exists in adata.obs
    if groupby not in adata.obs:
        raise ValueError(f"Column '{groupby}' not found in adata.obs.")
    
    cats = adata.obs[groupby].cat.categories.tolist()
    
    if groups == "all":
        groups_order = adata.obs[groupby].cat.categories
    elif isinstance(groups, (str, int)):
        raise ValueError("Specify a sequence of groups")
    else:
        groups_order = list(groups)
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
            if not set(groups_order) <= set(cats):
                raise ValueError(
                    f"groups = {groups} needs to be one of groupby = {cats}.")
        if reference != "rest" and reference not in set(groups_order):
            groups_order.append(reference)
    
    if reference != "rest" and reference not in adata.obs[groupby].cat.categories:
        raise ValueError(
            f"reference = {reference} needs to be one of groupby = {cats}."
        )
    
    return list(groups_order)




def compute_pvalue_bonferroni(
    pvalues : list
    ):
    """
    Compute Bonferroni-corrected p-values.

    Parameters
    ----------
    pvalues : list
        List of p-values to be corrected.

    Returns
    -------
    list
        Bonferroni-corrected p-values.
    """
    
    _, pvals_corrected, _, _ = multipletests(pvalues, method='bonferroni')
    return pvals_corrected
