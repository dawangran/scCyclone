# -*- coding: utf-8 -*-
"""
@File    :   _iso_entropy.py
@Time    :   2024/11/08
@Author  :   Dawn
@Version :   1.0
@Desc    :   Isoform entropy for scCyclone
"""



import numpy as np
import pandas as pd
import anndata as ad
from typing import Union
import scanpy as sc

import warnings
warnings.filterwarnings('ignore')


def _compute_score(G, T):
    """
    Calculate the score of a transcript based on the Kullback-Leibler divergence.

    Parameters:
    G (numpy.ndarray): The ground truth distribution.
    T (numpy.ndarray): The predicted distribution.

    Returns:
    float: The score of the transcript, normalized by the maximum possible score.
    """
    epsilon = 1e-10  # A small value to avoid division by zero in the log function
    DKL = 0  # Initialize the Kullback-Leibler divergence sum

    # Iterate over each sample in the predicted distribution
    for i in range(T.shape[0]):
        Ti = T[i]  # Get the i-th sample of the predicted distribution
        # Calculate the Kullback-Leibler divergence for the i-th sample
        DKL_Ti_G = np.sum(Ti * np.log((Ti + epsilon) / (G + epsilon)))
        DKL += DKL_Ti_G  # Add the divergence to the total sum

    # Calculate the mean Kullback-Leibler divergence
    DKL_mean = DKL / T.shape[0]
    
    # Calculate the maximum possible Kullback-Leibler divergence for normalization
    DKL_max = -np.log(T.shape[0]) / T.shape[0]
    
    # Calculate the score by normalizing the mean divergence by the maximum divergence
    score = DKL_mean / DKL_max
    return score



def isoform_entropy_score(
    adata: ad.AnnData,
    groupby: str,
    gene_list: list,
    groups: Union[str, list] = "all",
    var_name: str = "gene_name"
    ):
    """
    Calculate the isoform entropy score for a given list of genes.

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
    
    # Initialize a dictionary to store the isoform entropy scores
    isoform_score = {}
    percent_cell = {}
    
    # Iterate over each gene in the gene list
    for g in gene_list:
        # Get the indices of isoforms for the current gene
        isoform_list = adata_sub.var[adata_sub.var[var_name] == g].index.to_list()
        if len(isoform_list)==0:
            score=None
            percent=0  
         
        else:
            # Extract the data for the current gene's isoforms and convert to a DataFrame
            data = adata_sub[:, isoform_list].to_df()

            if data.sum(0).sum()>0:
            
                # Calculate the sum of isoform expressions for each cell
                data['gene'] = data.sum(axis=1)

                if len(isoform_list)==1:
                    percent = data['gene'].ne(0).sum()/len(data['gene'])
                    score=0
                else:
                    # Calculate the normalized gene score for each cell
                    data["gene_score"] = data['gene'] / data['gene'].sum()
                    
                    # Calculate the normalized isoform score for each isoform
                    for i in isoform_list:
                        data["{}_score".format(i)] = data[i] / data['gene'].sum()
                    
                    # Prepare the predicted distribution (T) and ground truth distribution (G)
                    T_list = []
                    for i in isoform_list:
                        T_list.append(data["{}_score".format(i)].to_list())
                    G = np.array(data['gene_score'].to_list())
                    
                    # Calculate the transcript score using the provided transcript_score function
                    score = _compute_score(G, np.array(T_list))
                    percent = data['gene'].ne(0).sum()/len(data['gene'])
            elif data.sum(0).sum()==0:
                score=None
                percent=0
            else:
                raise ValueError(f"adata.X cannot contain negative number")
            
        # Store the isoform entropy score for the current gene
        isoform_score[g] = score
        percent_cell[g] = percent
    # Return the dictionary of isoform entropy scores
    return isoform_score,percent_cell