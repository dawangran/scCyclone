# -*- coding: utf-8 -*-
"""
@File    :   _iso_analysis.py
@Time    :   2024/09/01
@Author  :   Dawn
@Version :   1.0
@Desc    :   Isoform annotation for scCyclone
"""

import os
import logging
import pandas as pd
import anndata as ad


_feature_dict = {
    "cpc2": ["coding_probability","label"],
    "deeploc2": ["Localizations", "Signals"],
    "sqanti3": ["chrom", "strand", "length", "exons", "structural_category", "CDS_length", "CDS_start", "CDS_end", "CDS_genomic_start", "CDS_genomic_end", "predicted_NMD"]   
}



def add_sqanti3(
    adata: ad.AnnData, 
    sqanti3_result_path: str):
    
    """
    Add SQANTI3 information to AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    sqanti3_result_path (str): Path to the SQANTI3 result file.

    Returns:
    ad.AnnData: Updated AnnData object with SQANTI3 information.
    """

    # Check if the SQANTI3 result file exists
    if not os.path.exists(sqanti3_result_path):
        logging.error(f"Transcript info file {sqanti3_result_path} does not exist.")
    
    # Read SQANTI3 result file
    logging.info(f"Reading transcript info from {sqanti3_result_path}")
    t_info = pd.read_csv(sqanti3_result_path, sep="\t")
    
    # Define columns to merge from SQANTI3 result
    colnames = _feature_dict['sqanti3']
    
    # Merge each column from SQANTI3 result with isoform data
    for colname in colnames:
        adata.var[colname] = pd.merge(adata.var, t_info, on='isoform', how='left')[colname].tolist()
        
    return adata


def add_gtf(
    adata: ad.AnnData, 
    sqanti3_gtf_path: str):
    
    """
    Add SQANTI3 GTF to AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    sqanti3_gtf_path (str): Path to the SQANTI3 GTF file.

    Returns:
    ad.AnnData: Updated AnnData object with SQANTI3 GTF.
    """
    
    #Check if the SQANTI3 gtf file exists
    if not os.path.exists(sqanti3_gtf_path):
        logging.error(f"Transcript info file {sqanti3_gtf_path} does not exist.")
        
    data=pd.read_csv(sqanti3_gtf_path,sep="\t",header=None)
    data['isoform']=[i[1:-2] for i in data[8].str.split(" ",expand=True)[1]]
    data=data[data[2]=="exon"]
    data=data.rename(columns={3:"exon_start",4:"exon_end"})
    data=data.groupby("isoform").agg({"exon_start":list,"exon_end":list})
    data=data[data.index.isin(list(adata.var.index))]
    data=data.loc[list(adata.var.index)]
    exon_dict={}
    for n,s,e in zip(data.index,data['exon_start'],data['exon_end']):
        exon_dict[n]=[s,e]
        
    adata.uns["isoform_exon"]=exon_dict
    




def add_cpc2(
    adata: ad.AnnData, 
    cpc2_result_path: str):
    
    """
    Add CPC2 information to AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    cpc2_result_path (str): Path to the CPC2 result file.

    Returns:
    ad.AnnData: Updated AnnData object with CPC2 information.
    """

    # Check if the CPC2 result file exists
    if not os.path.exists(cpc2_result_path):
        logging.error(f"Transcript info file {cpc2_result_path} does not exist.")

    # Read CPC2 result file and set 'isoform' column
    t_info = pd.read_csv(cpc2_result_path, index_col=0)
    t_info['isoform'] = list(t_info.index)
    t_info = t_info.replace({"coding":True,"nocoding":False})


    # Merge each column from CPC2 result with isoform data
    
    colnames = _feature_dict['cpc2']
    
    for colname in colnames:
        adata.var[colname] = pd.merge(adata.var, t_info, on='isoform', how='left')[colname].tolist()


    return adata
  


def add_pfam(
    adata: ad.AnnData, 
    pfam_result_path: str):
    
    """
    Add PFAM domain information to AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    pfam_result_path (str): Path to the PFAM domain result file.

    Returns:
    ad.AnnData: Updated AnnData object with PFAM domain information.
    """

    # Check if the PFAM domain result file exists
    if not os.path.exists(pfam_result_path):
        logging.error(f"Transcript info file {pfam_result_path} does not exist.")

    # Read PFAM domain result file
    t_info = pd.read_csv(pfam_result_path, index_col=0)
    t_info['isoform'] = list(t_info.index)
    
    # Filter out "RSB_motif" from the PFAM results
    t_info = t_info[t_info['hmm_name'] != "RSB_motif"]
    
    # Group PFAM results by isoform to aggregate domain count and list of domains
    t_info_type = t_info.groupby("isoform").agg({"clan": "count", "hmm_list": list})
    
    # Merge domain count and domain list with isoform data in adata
    adata.var["domain_number"] = pd.merge(adata.var, t_info_type, on='isoform', how='left')["clan"].tolist()
    adata.var["hmm_list"] = pd.merge(adata.var, t_info_type, on='isoform', how='left')["hmm_list"].tolist()
    
    return adata


def add_deepLoc2(adata: ad.AnnData, deeploc2_result_path: str):
    """
    Add DeepLoc2 information to AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    deeploc2_result_path (str): Path to the DeepLoc2 result file.

    Returns:
    ad.AnnData: Updated AnnData object with DeepLoc2 information.
    """

    # Check if the DeepLoc2 result file exists
    if not os.path.exists(deeploc2_result_path):
        logging.error(f"Transcript info file {deeploc2_result_path} does not exist.")

    # Read DeepLoc2 result file
    t_info = pd.read_csv(deeploc2_result_path, index_col=0)
    t_info['isoform'] = list(t_info.index)

    # Define DeepLoc2 feature columns to merge
    colnames = _feature_dict['deeploc2']

    # Merge DeepLoc2 information with isoform data for each feature column
    for colname in colnames:
        adata.var[colname] = pd.merge(adata.var, t_info, on='isoform', how='left')[colname].tolist()

    return adata
    


def add_custom(adata: ad.AnnData, custom_result_path: str, feature_colname: str):
    """
    Add custom feature information to AnnData object.

    Parameters:
    adata (ad.AnnData): Annotated data object.
    custom_result_path (str): Path to the custom feature result file.
    feature_colname (str): Name of the custom feature column.

    Returns:
    ad.AnnData: Updated AnnData object with custom feature information.
    """

    # Check if the custom feature result file exists
    if not os.path.exists(custom_result_path):
        logging.error(f"Transcript info file {custom_result_path} does not exist.")
        return adata

    # Read custom feature result file
    t_info = pd.read_csv(custom_result_path, index_col=0)
    t_info['isoform'] = list(t_info.index)

    # Merge custom feature information with isoform data based on the specified column name
    adata.var[feature_colname] = pd.merge(adata.var, t_info, on='isoform', how='left')[feature_colname].tolist()

    return adata