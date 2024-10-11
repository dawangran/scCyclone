# -*- coding: utf-8 -*-
"""
@File    :   _rank_ifs_groups.py
@Time    :   2024/09/04
@Author  :   Dawn
@Version :   1.0
@Desc    :   DTU for scCyclone
"""


import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import logging
import anndata as ad
from typing import Union

from . import _utils
from ..get import get




def _filter_iso(
    adata: ad.AnnData,
    groupby: str,
    groups: Union[str, list] = 'all',
    reference: str = 'rest',
    percent: float = 0.1,
    var_name: str = "gene_name"):
    
    """
    Filter genes based on the provided criteria and return the updated AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    groupby (str): Grouping variable.
    groups (Union[str, list]): Specific groups to consider.
    reference (str): Reference group.
    var_name (str): Name of the variable (e.g., gene_name).
    percent (float): Percentage threshold for filtering genes.

    Returns:
    ad.AnnData: Updated annotated data object after gene filtering.
    """
    

    # Check if groupby column exists in adata.obs
    if groupby not in adata.obs.columns:
        raise ValueError(f"cell_label '{groupby}' not found in adata.obs")
    
    # Check if percent is between 0 and 1
    if not (0 <= percent <= 1):
        raise ValueError("percent should be between 0 and 1")
    
    # Check if var_name exists in adata.var
    if var_name not in adata.var:
        raise ValueError(f"Column '{var_name}' not found in adata.var.")

    # Get the order of groups
    groups_order = _utils.check_groups(adata, groupby, groups, reference)
    
    # Filter genes based on counts
    gene_list_2 = list(pd.DataFrame(adata.var[var_name].value_counts())[pd.DataFrame(adata.var[var_name].value_counts())[var_name] > 1].index)
    adata = adata[:, adata.var[var_name].isin(gene_list_2)]

    iso_list = []
    for i in groups_order:
        sub_adata = adata[adata.obs[groupby] == i]
        sc.pp.filter_genes(sub_adata, min_cells=int(sub_adata.shape[0] * percent))
        iso_list += sub_adata.var_names.to_list()
        
    iso_list = list(set(iso_list)) 
    gene_list = list(set(adata.var[adata.var.index.isin(iso_list)][var_name]))
    
    adata = adata[:, adata.var[var_name].isin(gene_list)]   
    
        
    return adata

        

def _generate_gene_iso(
    adata: ad.AnnData,
    var_name: str = "gene_name"):
    """
    Generate a dictionary mapping gene names to their corresponding indices.

    Parameters:
    adata (AnnData): Annotated data object.
    var_name (str): Name of the variable containing gene names.

    Returns:
    dict: A dictionary mapping gene names to their indices.
    """

    # Check if var_name exists in adata.var columns
    if var_name not in adata.var.columns:
        raise ValueError(f"var_name '{var_name}' not found in adata.var")

    # Extract gene data based on var_name
    gene_data = adata.var[[var_name]]
    
    gene_iso = {}
    for k, v in zip(gene_data[var_name], gene_data.index):
        if k not in gene_iso:
            gene_iso[k] = [v]
        else:
            gene_iso[k].append(v)

    return gene_iso
    




def _generate_IF_bulk_matrix(
    adata: ad.AnnData,
    gene_iso: dict):
    """
    Generate an isoform fraction bulk matrix based on the provided gene isoform dictionary.

    Parameters:
    adata (AnnData): Annotated data object.
    gene_iso (dict): Dictionary mapping gene names to their corresponding indices.

    Returns:
    pd.DataFrame: DataFrame representing the isoform fraction bulk matrix.
    """

    # Calculate the sum of expression values for each gene across all cells
    data = pd.DataFrame(adata.to_df().sum(0)).T
    
    data_IF_list = []
    for _, isoforms in gene_iso.items():
        data_sub = data[isoforms]
        row_sums = data_sub.sum(axis=1)
        data_IF_sub = data_sub.div(row_sums, axis=0)
        data_IF_list.append(data_IF_sub)

    # Concatenate the isoform fraction data and transpose the matrix
    data_IF = pd.concat(data_IF_list, axis=1).T.fillna(0)
            
    return data_IF
    



def _generate_IF_adata(
    adata: ad.AnnData,
    var_name: str='gene_name'):
    """
    Generate an AnnData object with isoform fraction data based on the provided gene_name.

    Parameters:
    adata (AnnData): Annotated data object.
    var_name (str): Name of the variable containing gene names.

    Returns:
    AnnData: AnnData object with isoform fraction data.
    """

    # Convert AnnData to DataFrame
    data = adata.to_df().T
    data[var_name] = adata.var[var_name].to_list()

    # Group data by gene and aggregate by sum
    data_gene = data.groupby(var_name).agg(sum)

    # Initialize DataFrame for isoform fraction data
    data_IF = pd.DataFrame()

    # Calculate isoform fraction for each gene
    for i, j in enumerate(data.index):
        sub_data = data.iloc[i, :-1].values
        gene = data.iloc[i, -1]
        sub_gene = data_gene.loc[gene].values
        data_IF[j] = list(sub_data / sub_gene)

        if i % 1000 == 0:
            print("Process successful for {}".format(i))

    # Transpose DataFrame and create AnnData object
    data_IF = data_IF.T
    data_IF.columns = data_gene.columns
    adata_IF = ad.AnnData(data_IF.T)
    adata_IF.obs = adata.obs

    # Replace NaN values with 0
    adata_IF.X = np.nan_to_num(adata_IF.X)

    # Copy variable annotations from input AnnData
    adata_IF.var = adata.var

    return adata_IF
    




def _generate_iso_rank(data_IF: pd.DataFrame, gene_iso: dict):
    """
    Generate a DataFrame with ranks for isoforms based on the provided gene_iso dictionary.

    Parameters:
    data_IF (pd.DataFrame): DataFrame representing the isoform fraction data.
    gene_iso (dict): Dictionary mapping gene names to their corresponding isoforms.

    Returns:
    pd.DataFrame: DataFrame with ranks for isoforms.
    """

    data_rank_list = []

    for gene, isoforms in gene_iso.items():
        if not all(isoform in data_IF.index for isoform in isoforms):
            logging.warning(f"Not all isoforms of gene '{gene}' are present in data_IF")
            continue
        
        data_sub = data_IF.loc[isoforms]
        data_rank_sub = data_sub.apply(lambda column: column.rank(ascending=False, method='min'), axis=0)
        data_rank_list.append(data_rank_sub)
    
    data_rank = pd.concat(data_rank_list, 0)      
            
    data_rank = data_rank.astype(int)
    
    return data_rank
    

        
        

def _compute_dif(data_IF_ref: pd.DataFrame, data_IF_target: pd.DataFrame):
    """
    Compute the difference between two DataFrame objects containing isoform fraction data.

    Parameters:
    data_IF_ref (pd.DataFrame): DataFrame containing isoform fraction data for the reference.
    data_IF_target (pd.DataFrame): DataFrame containing isoform fraction data for the target.

    Returns:
    list: List of differences between the two DataFrames.
    """
    
    # Compute the difference between the two DataFrames
    dIF_list = (data_IF_target - data_IF_ref)[0].to_list()
    
    return dIF_list

        
        


def _compute_rank(data_rank_ref: pd.DataFrame, data_rank_target: pd.DataFrame):
    """
    Compute the rank comparison between two DataFrame objects containing ranks.

    Parameters:
    data_rank_ref (pd.DataFrame): DataFrame containing rank data for the reference.
    data_rank_target (pd.DataFrame): DataFrame containing rank data for the target.

    Returns:
    tuple: Tuple containing lists of rank differences, state of changes, and first ranks.
    """
    
    # Flatten the rank DataFrames and create a list of tuples
    dr_list = [(target_rank, ref_rank) for target_rank, ref_rank in zip(data_rank_target.values.flatten().tolist(), data_rank_ref.values.flatten().tolist())]

    # Determine the state of change for each rank comparison
    dr_state_list = []
    for target_rank, ref_rank in dr_list:
        if target_rank == ref_rank:
            state = "normal"
        elif target_rank > ref_rank:
            state = "up"
        else:
            state = "down"
        dr_state_list.append(state)

    # Identify if the target rank is the first rank
    dr_first_list = [target_rank == 1 for target_rank, _ in dr_list]

    return dr_list, dr_state_list, dr_first_list
    

    



def _compute_pvalue(adata_ref: ad.AnnData, adata_target: ad.AnnData):
    """
    Compute p-values using the Mann-Whitney U test between corresponding columns of two DataFrames.

    Parameters:
    data_ref_IF (pd.DataFrame): DataFrame containing reference isoform fraction data.
    data_target_IF (pd.DataFrame): DataFrame containing target isoform fraction data.

    Returns:
    list: List of p-values for each pair of corresponding columns.
    """

    pval_list = []
    data_ref=adata_ref.to_df()
    data_target=adata_target.to_df()
    

    for col in adata_ref.var.index:
        group1 = data_ref[col].to_list()
        group2 = data_target[col].to_list()
        _, p_value = mannwhitneyu(group1, group2, alternative='two-sided')
        pval_list.append(p_value)
    
    return pval_list
    


def _compute_proportion(adata_ref: ad.AnnData, adata_target: ad.AnnData):
    """
    Calculate the difference in the proportion of non-zero data for each column between two AnnData objects.

    Parameters:
    adata_ref (ad.AnnData): The reference AnnData object.
    adata_target (ad.AnnData): The target AnnData object.

    Returns:
    list: A list of differences in the proportion of non-zero data for each column.
    """
    
    # Convert AnnData objects to DataFrames
    data_ref = adata_ref.to_df()
    data_target = adata_target.to_df()
    
    # Calculate the number of non-zero data points for each column in the reference dataset
    non_zero_counts_ref = (data_ref != 0).sum()
    
    # Calculate the total number of data points for each column in the reference dataset
    total_counts_ref = data_ref.shape[0]
    
    # Calculate the proportion of non-zero data for each column in the reference dataset
    non_zero_ratios_ref = non_zero_counts_ref / total_counts_ref
    
    # Calculate the number of non-zero data points for each column in the target dataset
    non_zero_counts_target = (data_target != 0).sum()
    
    # Calculate the total number of data points for each column in the target dataset
    total_counts_target = data_target.shape[0]
    
    # Calculate the proportion of non-zero data for each column in the target dataset
    non_zero_ratios_target = non_zero_counts_target / total_counts_target
    
    # Calculate the difference in the proportion of non-zero data between the two datasets for each column
    dpr_list = (non_zero_ratios_target - non_zero_ratios_ref).tolist()
    
    return dpr_list

    

def rank_ifs_groups(
    adata: ad.AnnData,
    groupby: str,
    groups: Union[str, list] = "all",
    reference: Union[str, list] = "rest",
    key_added:  Union[None,str] = None,
    percent: float = 0.1,
    var_name: str = "gene_name"):
    """
    Rank if for characterizing groups.

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
    percent (float): Percentage threshold for filtering isoform.
         var_name (str): Name of the variable (e.g., gene_name).

    Returns
    -------
    ad.AnnData
        Annotated data matrix with rank information stored in adata.uns[key_added].
    """

    # Check and get the order of groups
    groups_order = _utils.check_groups(adata, groupby, groups, reference)

    
    # Set a default key if not provided
    if key_added is None:
        key_added = "rank_ifs_groups"
    
    # Initialize parameters in adata.uns
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = {"groupby": groupby, "reference": reference}

    # Filter data based on groups and reference
    adata_filter = _filter_iso(adata, groupby=groupby, groups=groups_order, reference=reference, var_name=var_name, percent=percent)
    gene_iso = _generate_gene_iso(adata_filter, var_name=var_name)
    iso_list = [val for sublist in gene_iso.values() for val in sublist]
    adata_filter = adata_filter[:, iso_list]
    
    # Initialize dictionaries to store results
    data_iso_dict = {}
    data_dif_dict = {}
    data_dr_dict = {}
    data_dr_state_dict = {}
    data_dr_first_dict = {}
    data_pval_dict = {}
    data_pval_adj_dict = {}
    data_dpr_dict = {}
    var_name_dict = {}
    
    
    # Iterate over groups
    for i in groups_order:
        print("Group {} start!".format(i))
        target_groups = [i]
        ref_groups = [x for x in groups_order if x != i] if reference == "rest" else [reference]

        adata_target = adata_filter[adata_filter.obs[groupby].isin(target_groups)]
        adata_ref = adata_filter[adata_filter.obs[groupby].isin(ref_groups)]
        print("Generate IF matrix...")
        data_target_IF_matrix = _generate_IF_bulk_matrix(adata_target, gene_iso=gene_iso)
        data_ref_IF_matrix = _generate_IF_bulk_matrix(adata_ref, gene_iso=gene_iso)

        print("Generate rank matrix...")
        data_target_rank = _generate_iso_rank(data_target_IF_matrix, gene_iso=gene_iso)
        data_ref_rank = _generate_iso_rank(data_ref_IF_matrix, gene_iso=gene_iso)

        print("Generate IF adata...")
        adata_target_IF = _generate_IF_adata(adata_target, var_name=var_name)
        adata_ref_IF = _generate_IF_adata(adata_ref, var_name=var_name)
        
        print("Compute dIF...")
        dif_list = _compute_dif(data_ref_IF_matrix, data_target_IF_matrix)
        dr_list, dr_state_list, dr_first_list = _compute_rank(data_target_rank, data_ref_rank)
        
        print("Compute pvalue...")
        pval_list = _compute_pvalue(adata_ref_IF, adata_target_IF)
        pval_adj_list = _utils.compute_pvalue_bonferroni(pval_list)
        
        print("Compute proportion...")
        dpr_list = _compute_proportion(adata_ref_IF, adata_target_IF)
        
        data_iso_dict[i] = iso_list
        data_dif_dict[i] = dif_list
        data_dr_dict[i] = dr_list
        data_dr_state_dict[i] = dr_state_list
        data_dr_first_dict[i] = dr_first_list
        data_pval_dict[i] = pval_list
        data_pval_adj_dict[i] = pval_adj_list
        data_dpr_dict[i] = dpr_list 
        var_name_dict[i] = adata_target.var[var_name].to_list()
        
        print("Group {} complete!".format(i))
        print("-----------------------------------------")
    
    # Convert dictionaries to structured arrays
    name_data = pd.DataFrame(data_iso_dict).to_records(index=False)
    dif_data = pd.DataFrame(data_dif_dict).to_records(index=False)
    dr_data = pd.DataFrame(data_dr_dict).to_records(index=False)
    dr_state_data = pd.DataFrame(data_dr_state_dict).to_records(index=False)
    dr_first_data = pd.DataFrame(data_dr_first_dict).to_records(index=False)
    pval_data = pd.DataFrame(data_pval_dict).to_records(index=False)
    pval_adj_data = pd.DataFrame(data_pval_adj_dict).to_records(index=False)
    dpr_data = pd.DataFrame(data_dpr_dict).to_records(index=False)
    var_name_data = pd.DataFrame(var_name_dict).to_records(index=False)
    
    
    # Store results in adata.uns
    adata.uns[key_added]['names'] = name_data
    adata.uns[key_added]['dif'] = dif_data
    adata.uns[key_added]['dr'] = dr_data
    adata.uns[key_added]['dr_state'] = dr_state_data
    adata.uns[key_added]['dr_first'] = dr_first_data
    adata.uns[key_added]['pvals'] = pval_data
    adata.uns[key_added]['pvals_adj'] = pval_adj_data
    adata.uns[key_added]['dpr'] = dpr_data
    adata.uns[key_added]['gene_name'] = var_name_data
    
    return adata

        
   