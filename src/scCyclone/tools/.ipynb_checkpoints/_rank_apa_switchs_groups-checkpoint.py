# -*- coding: utf-8 -*-
"""
@File    :   _rank_apa_swiths_groups.py
@Time    :   2025/06/17 
@Author  :   Dawn
@Version :   1.0
@Desc    :   APA Switch for scCyclone
"""


import pandas as pd
import anndata as ad
from typing import Union
import numpy as np


from ..get import get

        

def _convert_table(data):
    """
    Convert the data to a specific table format.

    Parameters:
    ----------
    
    data (DataFrame): Input data to be converted.

    Returns:
    ----------
    
    DataFrame: Converted data in a specific format.
    """
    # Initialize an empty list to store the converted data
    converted_data = []

    # Iterate over the original data
    for group, gene_name, sig, name in zip(data['group'], data['gene_name'], data['sig'], data['names']):
    
        # Extract all 'u' positions and values
        u_indices = [i for i, s in enumerate(sig) if s == 'up']
        u_values = [name[i] for i in u_indices]
        
        # Extract all 'd' positions and values
        d_indices = [i for i, s in enumerate(sig) if s == 'down']
        d_values = [name[i] for i in d_indices]
        
        # Pair each 'u' value with all 'd' values
        for u_value in u_values:
            for d_value in d_values:
                converted_data.append({'group': group, 'gene_name': gene_name, 'apaUpregulated': u_value, 'apaDownregulated': d_value})

    # Create a DataFrame
    df = pd.DataFrame(converted_data)
    return df



def rank_apa_switchs_groups(
    adata: ad.AnnData,
    group: Union[None, str, list] = None,
    key: str = "rank_apas_groups",
    key_added: Union[None, str] = None,
    pval_cutoff: Union[None, float] = 0.05,
    dpr_cutoff: Union[None, float] = None,
    tpr_cutoff: Union[None, float] = None,
    rpr_cutoff: Union[None, float] = None,
    abs_min_dpui: float = 0,
    abs_max_dpui: float = 1,
    ):
    """
    Perform ranking and switching analysis on groups of interest in AnnData.

    Parameters:
    ----------
    
    adata (ad.AnnData): Anndata object containing the data.
    group (Union[None, str, list]): Groups of interest for analysis.
    key (str): Key for the data.
    key_added (Union[None, str]): Additional key to be added.
    pval_cutoff (float): P-value cutoff for analysis.
    dpr_cutoff (float): dpr cutoff for analysis.
    tpr_cutoff (float): tpr cutoff for analysis.
    rpr_cutoff (float): rpr cutoff for analysis.
    abs_min_dpui (float): Minimum difference for analysis.
    abs_max_dpui (float): Maximum difference for analysis.


    Returns:
    ----------
    
    switch_data: Modified AnnData object after switch analysis.
    """
    
    # Set a default key if not provided
    if key_added is None:
        key_added = "rank_apa_switchs_groups"
    
    # Initialize parameters in adata.uns
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = adata.uns[key]['params']

    
    # Validate min_dif and max_dif values
    if not (0 <= abs_min_dpui <= 1) or not (0 <= abs_max_dpui <= 1):
        raise ValueError("abs_min_dpui and abs_max_dpui must be between 0 and 1.")
    
    # Ensure group is a list
    group = [group] if isinstance(group, str) else group or list(adata.uns[key]["names"].dtype.names)
    
    # Get DataFrame with rank data
    dPUI_data = get.rank_apas_groups_df(
        adata=adata, group=group, key=key, pval_cutoff=pval_cutoff,
        min_dpui=abs_min_dpui, max_dpui=abs_max_dpui,
        dpr_cutoff=dpr_cutoff,
        tpr_cutoff=tpr_cutoff,
        rpr_cutoff=rpr_cutoff,
        compare_abs=True)
    
    # Assign 'up' or 'down' based on 'dif' values
    dPUI_data['sig'] = ['up' if x > 0 else 'down' if x < 0 else '' for x in dPUI_data['dpui']]
    
    # Group by 'group' and 'gene_name' and aggregate 'sig' and 'names'
    dPUI_data = dPUI_data.groupby(["group", "gene_name"]).agg({"sig": list, "names": list}).reset_index()
    dPUI_data=dPUI_data.dropna()
    
    # Filter out rows with only one unique value in 'sig'
    dPUI_data = dPUI_data[dPUI_data['sig'].apply(lambda x: len(set(x)) > 1)]
    
    # Convert data to switch format
    switch_data = _convert_table(dPUI_data)
    
    adata.uns[key_added]['value'] = switch_data.to_records(index=False)
    
    
    # Return the modified AnnData object
    return adata


        
def rank_apa_switch_consequences_groups(
    adata: ad.AnnData, 
    var_name_list: list,
    key: str = "rank_apa_switchs_groups",
    key_added: Union[None, str] = None, 
    ):
    """
    This function ranks switch consequences groups based on specified variable names.

    Parameters:
    ----------
    
    adata (ad.AnnData): Anndata object containing variable information.
    var_name_list (list): List of variable names to be considered.
    key (str): Key for the data.
    key_added (Union[None, str]): Additional key to be added.

    Returns:
    ----------

    pd.DataFrame: Updated switch_data DataFrame.
    pd.DataFrame: Summary of switch data based on variable groups.
    """

    # Check if all variable names in var_name_list are present in adata.var.columns
    if not set(var_name_list).issubset(adata.var.columns):
        raise ValueError("var_name_list error")
    
    # Set a default key if not provided
    if key_added is None:
        key_added = "rank_apa_switchs_consequences_groups"
    
    # Initialize parameters in adata.uns
    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = adata.uns[key]['params']
    
    # Split variable names into numeric and string categories
    var_name_number_list = [i for i in var_name_list if isinstance(adata.var[i][0], (np.float64, np.int32, np.int64))]
    var_name_str_list = [i for i in var_name_list if isinstance(adata.var[i][0], np.bool_)]  
    
    switch_data_summary_list = []
    
    switch_data=get.rank_switchs_groups_df(adata,key)
    
    # Process numeric variables
    if var_name_number_list:
        for i in var_name_number_list:
            tmp_dict = {k: v for k, v in zip(adata.var['apa'], adata.var[i])}
            _switch_data = switch_data.replace(tmp_dict)

            # Determine "more" or "less" based on comparisons
            _switch_data[i] = ["more" if up > down else "less" if up < down else None for up, down in zip(_switch_data['apaUpregulated'], _switch_data['apaDownregulated'])]

            _switch_data_summary = _switch_data.groupby('group')[i].value_counts().unstack(fill_value=0)
            _switch_data_summary=_switch_data_summary.rename(columns={"more": "nUP", "less": "nDown"})
            _switch_data_summary['total'] = _switch_data_summary['nUP'] + _switch_data_summary['nDown']
            _max_total = _switch_data_summary['total'].max()
            _switch_data_summary['percentage'] = _switch_data_summary['total'] / _max_total
            _switch_data_summary['log2fc'] = np.log2(_switch_data_summary["nUP"] / _switch_data_summary["nDown"])
            _switch_data_summary['feature'] = i
            switch_data_summary_list.append(_switch_data_summary)
            

    # Process string variables
    if var_name_str_list:
        for i in var_name_str_list:
            tmp_dict = {k: v for k, v in zip(adata.var['apa'], adata.var[i])}
            _switch_data = switch_data.replace(tmp_dict)

            # Assign values based on conditions
            _switch_data[i] = [None if up == down else "more" if up == True else "less" for up, down in zip(_switch_data['apaUpregulated'], _switch_data['apaDownregulated'])]

            _switch_data_summary = _switch_data.groupby('group')[i].value_counts().unstack(fill_value=0)
            _switch_data_summary=_switch_data_summary.rename(columns={"more": "nUP", "less": "nDown"})
            _switch_data_summary['total'] = _switch_data_summary['nUP'] + _switch_data_summary['nDown']
            _max_total = _switch_data_summary['total'].max()
            _switch_data_summary['percentage'] = _switch_data_summary['total'] / _max_total 
            _switch_data_summary['log2fc'] = np.log2(_switch_data_summary["nUP"] / _switch_data_summary["nDown"])
            _switch_data_summary['feature'] = i
            switch_data_summary_list.append(_switch_data_summary)
        
    switch_data_summary = pd.concat(switch_data_summary_list)
    switch_data_summary = switch_data_summary.reset_index()
    
    
    adata.uns[key_added]['value'] = switch_data_summary.to_records(index=False)
    

    return adata
        
        
    
    
            
