# -*- coding: utf-8 -*-
"""
@File    :   utils.py
@Time    :   2024/05/29 13:39:55
@Author  :   Dawn
@Version :   1.0
@Desc    :   Get query for scCyclone
"""

import anndata as ad
import pandas as pd
from typing import Union
from packaging.version import Version
import joblib
import os
import logging
import numpy as np



__root_dir__ = os.path.abspath(os.path.dirname(__file__))


def rank_ifs_groups_df(
    adata: ad.AnnData,
    group: Union[None, str, list] = None,
    key: str = "rank_ifs_groups",
    pval_cutoff: float = 0.05,
    min_dif: float = 0,
    max_dif: float = 1,
    rank_state: Union[None, str] = None,
    first: bool = False,
    dpr_cutoff: float = 0.3,
    compare_abs: bool = False):
    """
    Rank and filter differential splicing events by group in an AnnData object.

    Parameters:
    adata (ad.AnnData): Anndata object containing the data.
    group (Union[None, str, list]): Groups of interest for analysis.
    key (str): Key for the data.
    pval_cutoff (float): P-value cutoff for analysis.
    min_dif (float): Minimum difference value.
    max_dif (float): Maximum difference value.
    rank_state (Union[None, str]): State for ranking ('up', 'down', 'normal', or None).
    first (bool): Flag to filter for the first occurrence.
    dpr_cutoff (float): dpr cutoff for analysis.
    compare_abs (bool): Flag to compare absolute values.

    Returns:
    pd.DataFrame: DataFrame containing ranked and filtered differential splicing events.
    """

    if not (0 <= min_dif <= 1):
        raise ValueError("min_dif must be between 0 and 1.")
    if not (0 <= max_dif <= 1):
        raise ValueError("max_dif must be between 0 and 1.")
    if not (0 <= dpr_cutoff <= 1):
        raise ValueError("max_dif must be between 0 and 1.")
    


    if rank_state not in ["up", "down", "normal", None]:
        raise ValueError("Invalid rank_state value. Please provide 'up', 'down', 'normal', or None.")

    if isinstance(group, str):
        group = [group]
    if group is None:
        group = list(adata.uns[key]["names"].dtype.names)

    colnames = ["names", "dif", "dr", "dr_state", "dr_first", "pvals", "pvals_adj", "dpr","gene_name"]

    d = [pd.DataFrame(adata.uns[key][c])[group] for c in colnames]
    d = pd.concat(d, axis=1, names=[None, "group"], keys=colnames)

    if Version(pd.__version__) >= Version("2.1"):
        d = d.stack(level=1, future_stack=True).reset_index()
    else:
        d = d.stack(level=1).reset_index()
    d["group"] = pd.Categorical(d["group"], categories=group)
    d = d.sort_values(["group", "level_0"]).drop(columns="level_0")

    if pval_cutoff is not None:
        d = d[d["pvals_adj"] < pval_cutoff]
        
    if dpr_cutoff is not None:
        d = d[d["dpr"] < dpr_cutoff]

    d = d[(abs(d["dif"]) > min_dif if compare_abs else d["dif"] > min_dif) & 
          (abs(d["dif"]) < max_dif if compare_abs else d["dif"] < max_dif)]

    if rank_state is not None:
        d = d[d["dr_state"] == rank_state]
    if first:
        d = d[d["dr_first"] == first]

    return d.reset_index(drop=True)



def rank_switchs_groups_df(
    adata: ad.AnnData,
    key: str = "rank_switchs_groups",
    ):
    """
    Rank and filter differential splicing events by group in an AnnData object.

    Parameters:
    adata (ad.AnnData): Anndata object containing the data.
    key (str): Key for the data.
   
    Returns:
    pd.DataFrame: DataFrame containing ranked and filtered differential splicing events.
    """
            
    d = pd.DataFrame(adata.uns[key]['value'])
    
    return d


def rank_switch_consequences_groups_df(
    adata: ad.AnnData,
    key: str = "rank_switch_consequences_groups",
    ):
    """
    Rank and filter differential splicing events by group in an AnnData object.

    Parameters:
    adata (ad.AnnData): Anndata object containing the data.
    key (str): Key for the data.
   
    Returns:
    pd.DataFrame: DataFrame containing ranked and filtered differential splicing events.
    """
            
    d = pd.DataFrame(adata.uns[key]['value'])
    
    return d

  
        
def rank_psis_groups_df(
    adata: ad.AnnData,
    group: Union[None, str, list] = None,
    key: str = "rank_psis_groups",
    pval_cutoff: float = 0.05,
    min_dpsi: float = 0,
    max_dpsi: float = 1,
    gene_symbols: Union[None, str] = None,
    compare_abs: bool = False
    ):
    
    """
    Rank and filter differential splicing events by group in an AnnData object.

    Parameters:
    adata (ad.AnnData): Anndata object containing the data.
    group (Union[None, str, list]): Groups of interest for analysis.
    key (str): Key for the data.
    pval_cutoff (float): P-value cutoff for analysis.
    min_dpsi (float): Minimum delta PSI value.
    max_dpsi (float): Maximum delta PSI value.
    gene_symbols (Union[None, str]): Gene symbol information.
    compare_abs (bool): Flag to compare absolute values.

    Returns:
    pd.DataFrame: DataFrame containing ranked and filtered differential splicing events.
    """

    if not (0 <= min_dpsi <= 1):
        raise ValueError("min_dpsi must be between 0 and 1.")
    if not (0 <= max_dpsi <= 1):
        raise ValueError("max_dpsi must be between 0 and 1.")

    if isinstance(group, str):
        group = [group]
    if group is None:
        group = list(adata.uns[key]["names"].dtype.names)

    colnames = ["names", "dpsi", "pvals", "pvals_adj"]

    d = [pd.DataFrame(adata.uns[key][c])[group] for c in colnames]
    d = pd.concat(d, axis=1, names=[None, "group"], keys=colnames)

    if Version(pd.__version__) >= Version("2.1"):
        d = d.stack(level=1, future_stack=True).reset_index()
    else:
        d = d.stack(level=1).reset_index()
    d["group"] = pd.Categorical(d["group"], categories=group)
    d = d.sort_values(["group", "level_0"]).drop(columns="level_0")

    if pval_cutoff is not None:
        d = d[d["pvals_adj"] < pval_cutoff]

    d = d[(abs(d["dpsi"]) > min_dpsi if compare_abs else d["dpsi"] > min_dpsi) & 
          (abs(d["dpsi"]) < max_dpsi if compare_abs else d["dpsi"] < max_dpsi)]

    if gene_symbols is not None:
        d = d.join(adata.var[gene_symbols], on="names")

    return d.reset_index(drop=True)



def psis_rmaps_df(
    event_list: list,
    type: str,
    gtf_file: str):
    
    if not os.path.exists(gtf_file):
        logging.error(f"Data file {gtf_file} does not exist.")

    
    gtf=pd.read_csv(gtf_file,sep="\t",header=None)
    gtf['isoform']=gtf[8].str.split('"',expand=True)[1].to_list()
    gtf['gene_id']=gtf[8].str.split('"',expand=True)[3].to_list()
    gtf=gtf[['gene_id','isoform',0,2,3,4,6]]
    gtf=gtf.rename(columns={1:"chr",2:"type",3:"start",4:"end",6:"stand"})
    gtf=gtf[gtf['type']=="exon"]
    psi_data=pd.DataFrame()
    chr_list=[]
    strand_list=[]
    
    if type=="SE":
        firstExonEnd_list=[]
        exonStart_list=[]
        exonEnd_list=[]
        secondExonStart_list=[]
        firstExonstart_list=[]
        secondExonEnd_list=[]
        for i in event_list:
            gene_id, chr, strand = i.split(";")[0], i.split(";")[1].split(":")[1], i.split(";")[1].split(":")[4]
            firstExonEnd, exonStart=map(int, i.split(";")[1].split(":")[2].split("-"))
            exonEnd, secondExonStart=map(int, i.split(";")[1].split(":")[3].split("-"))
            firstExonStart = gtf[(gtf['gene_id'] == gene_id) & (gtf['end'] == firstExonEnd)]['start'].iloc[0]
            secondExonEnd = gtf[(gtf['gene_id'] == gene_id) & (gtf['start'] == secondExonStart)]['end'].iloc[0]
            
            chr_list.append(chr)
            strand_list.append(strand)
            firstExonEnd_list.append(firstExonEnd)
            exonStart_list.append(exonStart)
            exonEnd_list.append(exonEnd)
            secondExonStart_list.append(secondExonStart)
            firstExonstart_list.append(firstExonStart)
            secondExonEnd_list.append(secondExonEnd)
            
        psi_data['chr']=chr_list
        psi_data['strand']=strand_list
        psi_data['firstExonStart']=firstExonstart_list
        psi_data['firstExonEnd']=firstExonEnd_list
        psi_data['exonStart']=exonStart_list
        psi_data['exonEnd']=exonEnd_list
        psi_data['secondExonStart']=secondExonStart_list
        psi_data['secondExonEnd']=secondExonEnd_list
        
    if type=="RI":
        upstreamExonStart_list=[]
        upstreamExonEnd_list=[]
        downstreamExonStart_list=[]
        downstreamExonEnd_list=[]
        for i in event_list:
            chr=i.split(";")[1].split(":")[1]
            strand=i.split(";")[1].split(":")[5]
            upstreamExonStart=i.split(";")[1].split(":")[2]
            upstreamExonEnd=i.split(";")[1].split(":")[3].split("-")[0]
            downstreamExonStart=i.split(";")[1].split(":")[3].split("-")[1]
            downstreamExonEnd=i.split(";")[1].split(":")[4]
            
            chr_list.append(chr)
            strand_list.append(strand)
            upstreamExonStart_list.append(upstreamExonStart)
            upstreamExonEnd_list.append(upstreamExonEnd)
            downstreamExonStart_list.append(downstreamExonStart)
            downstreamExonEnd_list.append(downstreamExonEnd)
            
        psi_data['chr']=chr_list
        psi_data['strand']=strand_list
        psi_data['riExonStart']=upstreamExonStart_list
        psi_data['riExonEnd']=downstreamExonEnd_list
        psi_data['upstreamExonStart']=upstreamExonStart_list
        psi_data['upstreamExonEnd']=upstreamExonEnd_list
        psi_data['downstreamExonStart']=downstreamExonStart_list
        psi_data['downstreamExonEnd']=downstreamExonEnd_list
        
    if type=="A3":
        flankingExonStart_list=[]
        flankingExonEnd_list=[]
        longExonStart_list=[]
        shortExonStart_list=[]
        shortExonEnd_list=[]
        for i in event_list:
            gene_id=i.split(";")[0]
            chr=i.split(";")[1].split(":")[1]
            strand=i.split(";")[1].split(":")[4]
            flankingExonEnd=i.split(";")[1].split(":")[2].split("-")[0]
            flankingExonStart=list(gtf[(gtf['gene_id']==gene_id) & (gtf['end']==int(flankingExonEnd))]['start'])[0]
            longExonStart=i.split(";")[1].split(":")[2].split("-")[1]
            shortExonStart=i.split(";")[1].split(":")[3].split("-")[1]
            shortExonEnd=list(gtf[(gtf['gene_id']==gene_id) & (gtf['start']==int(shortExonStart))]['end'])[0]
            
            chr_list.append(chr)            
            strand_list.append(strand)
            flankingExonStart_list.append(flankingExonStart)
            flankingExonEnd_list.append(flankingExonEnd)
            longExonStart_list.append(longExonStart)
            shortExonStart_list.append(shortExonStart)
            shortExonEnd_list.append(shortExonEnd)
            
        psi_data['chr']=chr_list
        psi_data['strand']=strand_list
        psi_data['longExonStart']=longExonStart_list
        psi_data['longExonEnd']=shortExonEnd_list
        psi_data['shortExonStart']=shortExonStart_list
        psi_data['shortExonEnd']=shortExonEnd_list
        psi_data['flankingExonStart']=flankingExonStart_list
        psi_data['flankingExonEnd']=flankingExonEnd_list
        

    if type=="A5":
        flankingExonStart_list=[]
        flankingExonEnd_list=[]
        longExonEnd_list=[]
        shortExonStart_list=[]
        shortExonEnd_list=[]
        for i in event_list:
            gene_id=i.split(";")[0]
            chr=i.split(";")[1].split(":")[1]
            strand=i.split(";")[1].split(":")[4]
            shortExonEnd=i.split(";")[1].split(":")[3].split("-")[0]
            shortExonStart=list(gtf[(gtf['gene_id']==gene_id) & (gtf['end']==int(shortExonEnd))]['start'])[0]
            longExonEnd=i.split(";")[1].split(":")[2].split("-")[0]
            flankingExonStart=i.split(";")[1].split(":")[2].split("-")[1]
            flankingExonEnd=list(gtf[(gtf['gene_id']==gene_id) & (gtf['start']==int(flankingExonStart))]['end'])[0]
            
            chr_list.append(chr)            
            strand_list.append(strand)
            flankingExonStart_list.append(flankingExonStart)
            flankingExonEnd_list.append(flankingExonEnd)
            longExonEnd_list.append(longExonEnd)
            shortExonStart_list.append(shortExonStart)
            shortExonEnd_list.append(shortExonEnd)
                                        
        psi_data['chr']=chr_list
        psi_data['strand']=strand_list
        psi_data['longExonStart']=shortExonStart_list
        psi_data['longExonEnd']=longExonEnd_list
        psi_data['shortExonStart']=shortExonStart_list
        psi_data['shortExonEnd']=shortExonEnd_list
        psi_data['flankingExonStart']=flankingExonStart_list
        psi_data['flankingExonEnd']=flankingExonEnd_list
        
    if type=="MX":
        firstExonStart_list=[]
        firstExonEnd_list=[]
        secondExonStart_list=[]
        secondExonEnd_list=[]
        upstreamExonStart_list=[]
        upstreamExonEnd_list=[]
        downstreamExonStart_list=[]
        downstreamExonEnd_list=[]
        
        for i in event_list:
            gene_id=i.split(";")[0]
            chr=i.split(";")[1].split(":")[1]
            strand=i.split(";")[1].split(":")[6]
            upstreamExonEnd=i.split(";")[1].split(":")[2].split("-")[0]
            upstreamExonStart=list(gtf[(gtf['gene_id']==gene_id) & (gtf['end']==int(upstreamExonEnd))]['start'])[0]
            firstExonStart=i.split(";")[1].split(":")[2].split("-")[1]
            firstExonEnd=i.split(";")[1].split(":")[3].split("-")[0]
            downstreamExonStart=i.split(";")[1].split(":")[3].split("-")[1]
            downstreamExonEnd=list(gtf[(gtf['gene_id']==gene_id) & (gtf['start']==int(downstreamExonStart))]['end'])[0]
            secondExonStart=i.split(";")[1].split(":")[4].split("-")[1]
            secondExonEnd=i.split(";")[1].split(":")[5].split("-")[0]
            
            chr_list.append(chr)            
            strand_list.append(strand)
            firstExonStart_list.append(firstExonStart)
            firstExonEnd_list.append(firstExonEnd)
            secondExonStart_list.append(secondExonStart)
            secondExonEnd_list.append(secondExonEnd)
            upstreamExonStart_list.append(upstreamExonStart)
            upstreamExonEnd_list.append(upstreamExonEnd)
            downstreamExonStart_list.append(downstreamExonStart)
            downstreamExonEnd_list.append(downstreamExonEnd)
            
        psi_data['chr']=chr_list
        psi_data['strand']=strand_list
        psi_data['1stExonstart']=firstExonStart_list
        psi_data['1stExonEnd']=firstExonEnd_list
        psi_data['2ndExonstart']=secondExonStart_list
        psi_data['2ndExonEnd']=secondExonEnd_list
        psi_data['upstreamExonStart']=upstreamExonStart_list
        psi_data['upstreamExonEnd']=upstreamExonEnd_list
        psi_data['downstreamExonStart']=downstreamExonStart_list
        psi_data['downstreamExonEnd']=downstreamExonEnd_list
        
    return psi_data




def psis_modal_df(adata: ad.AnnData, groupby: str, event_list: str, valid_cells: int = 25, groups: Union[str, list] = "all", pkl: Union[str, None] = None):
    """
    Calculate modal probabilities for events in a given AnnData object based on a pre-trained classifier.

    Parameters:
    adata (ad.AnnData): An AnnData object containing the data.
    groupby (str): The column in adata.obs used for grouping.
    event_list (str): List of events to calculate modal probabilities for.
    valid_cells (int): Minimum number of valid cells required to calculate modal probabilities.
    groups (Union[str, list]): Groups to consider for calculating modal probabilities. Defaults to "all".
    pkl (Union[str, None]): Path to a custom pre-trained classifier. Defaults to None.

    Returns:
    event_modal_data (pd.DataFrame): DataFrame with modal probabilities for each event.
    event_best_modal_data (pd.DataFrame): DataFrame with the best modal events identified.
    """
    
    # Load the pre-trained classifier
    if pkl is None:
        # Load the default pre-trained classifier
        clf = joblib.load(os.path.join(__root_dir__,"model","PSI_random_forest_model4.pkl"))

    else:
        # Load the specified classifier
        if not os.path.exists(pkl):
            logging.error(f"Data file {pkl} does not exist.")
            return None
        clf = joblib.load(pkl)

    # Check if the grouping column exists in adata.obs
    if groupby not in adata.obs.columns:
        raise ValueError(f"Grouping column '{groupby}' not found in adata.obs")

    # Get the categories for grouping
    
    cats = adata.obs[groupby].cat.categories
    if groups == "all":
        groups_order = cats
    else:
        groups_order = groups
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
        if not set(groups_order) <= set(cats):
            raise ValueError(f"Groups {groups} need to be one of the groupby categories {cats}")

    # Check if event_list is present in adata.var
    if len(set(event_list) & set(adata.var.index)) != len(set(event_list)):
        raise ValueError("Some events in event_list are not found in adata.var")

    # Filter data based on groups and event_list
    sub_adata = adata[adata.obs[groupby].isin(groups_order)]
    sub_adata = sub_adata[:,sub_adata.var.index.isin(event_list)]
    print(sub_adata)
    data = sub_adata.to_df()

    class_names = clf.classes_

    # Calculate modal probabilities for each event
    event_modal_dict = {}
    for i in data.columns:
        sub_data = data[[i]]
        sub_data = sub_data.dropna(subset=[i])
        if sub_data.shape[0] >= valid_cells:
            bins = np.linspace(0, 1, 11)
            sub_data['bin'] = pd.cut(sub_data[i], bins, include_lowest=True)
            bin_counts = sub_data['bin'].value_counts().sort_index()
            bin_proportions = bin_counts / bin_counts.sum()
            bin_proportions = np.array(bin_proportions).reshape(1, -1)
            probabilities = clf.predict_proba(bin_proportions)
            event_modal_dict[i] = probabilities[0].tolist()

    event_modal_data = pd.DataFrame(event_modal_dict)
    event_modal_data.index = class_names

    # Identify the best modal events
    best_event_modal = []
    for i in event_modal_data.columns:
        pro_data = event_modal_data[[i]]
        pro_data=pro_data.sort_values(i,ascending=False)
        if pro_data.iloc[0].values[0] > 0.4 and pro_data.iloc[0].values[0] > (1.5 * pro_data.iloc[1].values[0]):
            best_event_modal.append(pro_data.index[0])
        else:
            best_event_modal.append(None)

    event_best_modal_data = pd.DataFrame({"event": event_modal_data.columns.to_list(), "state": best_event_modal})

    return event_modal_data, event_best_modal_data