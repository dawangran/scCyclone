# -*- coding: utf-8 -*-
"""
@File    :   _apa_score.py
@Time    :   2025/6/17
@Author  :   Dawn
@Version :   1.0
@Desc    :   APA score for scCyclone
"""



import numpy as np
import pandas as pd
import anndata as ad
from typing import Union
import scanpy as sc

import warnings
warnings.filterwarnings('ignore')



def apa_score(
    adata: ad.AnnData,
    groupby: str,
    gene_list: list,
    groups: Union[str, list] = "all",
    var_name: str = "gene_name"
    ):
    """
    Calculate the apa score for a given list of genes.

    Parameters:
    adata (AnnData): The AnnData object containing gene expression data.
    groupby (str): The column name in `adata.obs` to group the data by.
    groups (list): The list of groups to consider from the `groupby` column.
    gene_list (list): The list of gene names for which to calculate the isoform entropy score.
    var_name (str): The column name in `adata.var` that contains the gene symbols.

    Returns:
    dict: A dictionary with gene names as keys and their corresponding isoform entropy scores as values.
    """

    if groupby not in adata.obs:
        raise ValueError(f"Column '{groupby}' not found in adata.obs.")
    
    cats = adata.obs[groupby].cat.categories.tolist()
    
    if groups == "all":
        groups_order = adata.obs[groupby].cat.categories
    else:
        if len(set(groups)&set(cats))==len(groups):
            groups_order=groups
        else:
            raise ValueError(f"groups not found in groupby.")

    
    # Subset the data to only include the specified groups
    adata_sub = adata[adata.obs[groupby].isin(groups_order)]
    adata_sub = adata_sub[:,adata_sub.var[var_name].isin(gene_list)]
    
    
    apa_id_translate_dict={}
    for g in gene_list:
        apa_data=adata_sub.var[adata_sub.var["gene_name"] == g]
        mean_len=apa_data['length'].mean()
        for l in apa_data[apa_data['length']>=mean_len].index:
            apa_id_translate_dict[l]="{}_long".format(g)
        for s in apa_data[apa_data['length']<mean_len].index:
            apa_id_translate_dict[s]="{}_short".format(g)
            
    adata_sub.var['apa_translate'] = list(adata_sub.var['apa'])
    adata_sub.var['apa_translate'] = adata_sub.var['apa_translate'].replace(apa_id_translate_dict)
    data=adata_sub.to_df().T
    data['apa_translate']=list(adata_sub.var['apa_translate'])
    data=data.groupby("apa_translate").agg(sum).T

    score_dict={}

    # Iterate over each gene in the gene list
    for g in gene_list:
        score=np.log2(data["{}_long".format(g)],data["{}_short".format(g)])
        score_dict[g]=score
        
    score_data=pd.DataFrame(score_dict)
    
    score_data = score_data.replace(-np.inf, -10)
    score_data = score_data.replace(np.inf, 10)
    adata_score=ad.AnnData(score_data)
    adata_score.obs=adata_sub.obs
    
    return adata_score