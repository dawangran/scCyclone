# -*- coding: utf-8 -*-
"""
@File    :   reader.py
@Time    :   2024/05/29 20:55:55
@Author  :   Dawn
@Version :   1.0
@Desc    :   Generate adata for scCyclone
"""
       
        
import os
import numpy as np
import pandas as pd
import anndata as ad
import logging
from typing import Union





def generate_Iso_adata(
    data_path: str,
    ):
    """
    Process Cyclone isoform data and save as AnnData objects.
    
    Parameters:
        data_path (str): Path to the isoform count data CSV file (isoquant result).
    Returns:
        adata (anndata.AnnData): Processed AnnData object.
    """
    
    # Check if the data file exists
    if not os.path.exists(data_path):
        logging.error(f"Data file {data_path} does not exist.")
        return None

  
    # Read isoform count data
    logging.info(f"Reading isoform count data from {data_path}")
    data = pd.read_csv(data_path, index_col=0)
    
    # Create AnnData object
    adata = ad.AnnData(data.T)
    adata.obs['batch'] = [i.rsplit("_", 1)[0] for i in adata.obs_names]    
    adata.var['isoform'] = adata.var.index
   
    
    # Remove entries with NaN gene names

    # adata = adata[:, adata.var.dropna(subset=['gene_id']).index]
    
    logging.info("Processing completed successfully.")
    
    return adata
    


def generate_PSI_adata(
    adata: ad.AnnData,
    event_info_path: str
    ):
    """
    Read Cyclone PSI (Percent Spliced In) data and save as AnnData object.

    Parameters:
        adata (anndata.AnnData): Input AnnData object(iso adata).
        event_info_path (str): Path to the event information file(suppa2 result).

    Returns:
        adata_as (anndata.AnnData): Processed AnnData object.
    """

    
    # Check if the event info file exists
    if not os.path.exists(event_info_path):
        logging.error(f"Transcript info file {event_info_path} does not exist.")
        return None

    # Read event information
    data = pd.read_csv(event_info_path, sep="\t")
    data['type'] = data['event_id'].str.split(";", expand=True)[1].str.split(":", expand=True)[0]
    data = data[['event_id', 'gene_id', 'type', 'alternative_transcripts', 'total_transcripts']]

    # Prepare PSI values
    cell = adata.obs.index.to_list()
    feature = data['event_id'].to_list()
    data["acceptor"] = data['alternative_transcripts'].str.split(",")
    data['total'] = data['total_transcripts'].str.split(",")   
    df = pd.DataFrame(index=cell, columns=feature)
    data_list = []

    for i in range(0, df.shape[1], 1000):
        start_col = i
        end_col = i + 1000 if i + 1000 <= df.shape[1] else df.shape[1]
        df_subset = df.iloc[:, start_col:end_col]

        for column in df_subset.columns:
            alternative_transcripts = data[data['event_id'] == column]['acceptor'].values[0]
            total_transcripts = data[data['event_id'] == column]['total'].values[0]
            
            alternative_transcripts_number = adata[:, adata.var_names.isin(alternative_transcripts)].X.sum(axis=1)
            total_transcripts_number = adata[:, adata.var_names.isin(total_transcripts)].X.sum(axis=1)
            df_subset[column] = alternative_transcripts_number / total_transcripts_number

        data_list.append(df_subset)

    data_as = pd.concat(data_list, axis=1)
    adata_as = ad.AnnData(data_as)
    adata_as.obs = adata.obs
    adata_as.var['gene_id'] = data['gene_id'].to_list()
    adata_as.var['type'] = data['type'].to_list()
    adata_as.var['alternative_transcripts'] = data['alternative_transcripts'].to_list()
    adata_as.var['total_transcripts'] = data['total_transcripts'].to_list()
    
    logging.info("Processing completed successfully.")
    
    return adata_as

        

def generate_Gene_adata(
    adata: ad.AnnData,
    var_name: str ='gene_name'
    ):
    """
    Aggregate isoform expression data to gene level and save as AnnData object.

    Parameters:
        adata (anndata.AnnData): Input AnnData object containing isoform expression data.
        var_name (str): Name of the column in `adata.var` that contains.

    Returns:
        adata_gene (anndata.AnnData): AnnData object containing aggregated gene expression data.
    """
    # Check if input AnnData object is provided

    
    # Check if gene_name column exists in adata.var
    if var_name not in adata.var:
        raise ValueError(f"Column '{var_name}' not found in adata.var.")

    # Aggregate isoform expression data to gene level
    gene_data = adata.to_df().T
    gene_data[var_name] = adata.var[var_name].to_list()
    gene_data = gene_data.groupby(var_name).agg(sum)
    
    # Check if any genes are aggregated
    if len(gene_data) == 0:
        raise ValueError("No genes found in the input AnnData object.")

    # Create AnnData object for gene expression
    adata_gene = ad.AnnData(gene_data.T)
    adata_gene.obs = adata.obs

    
    logging.info("Processing completed successfully.")
    
    return adata_gene


def generate_IF_adata(
    adata: ad.AnnData, 
    var_name: str ='gene_name',
    bulk: bool =False,
    obs_name: Union[str, None] = None
    ):
    """
    Generate an AnnData object with isoform fraction data.

    Parameters:
        adata (anndata.AnnData): Input AnnData object.
        gene_name (str, optional): Name of the column in `adata.var` that contains. Defaults to 'gene_name'.
        bulk (bool): IF bulk is true, the obs_name command is used to merge the obs_name command to generate the IF of the bulk level.
        obs_name (str): Name of the column in `adata.obs` that contains. Defaults to  None.
        
    Returns:
        adata_IF (anndata.AnnData): AnnData object with isoform fraction data.
    """
    
    
    # Check if gene_name column exists in adata.var
    if var_name not in adata.var:
        raise ValueError(f"Column '{var_name}' not found in adata.var.")
    
    if bulk:
        if obs_name not in adata.obs:
            raise ValueError(f"Column '{obs_name}' not found in adata.obs.")
        
        data = adata.to_df()
        data[obs_name] = adata.obs[obs_name].to_list()
        data_iso = data.groupby(obs_name).agg("sum")
        data_iso = data_iso.T
        data_iso[var_name] = adata.var[var_name].to_list()
        data_gene = data_iso.iloc[:, :-1]
        data_gene[var_name] = adata.var[var_name].to_list()
        data_gene = data_gene.groupby(var_name).agg("sum")
        data_IF = pd.DataFrame()
        for i, j in enumerate(data_iso.index):
            sub_data = data_iso.iloc[i, :-1].values
            gene = data_iso.iloc[i, -1]
            sub_gene = data_gene.loc[gene].values
            data_IF[j] = list(sub_data / sub_gene)
            
            if i % 1000 == 0:
                print("process succes {}".format(i))
                
        data_IF = data_IF.T
        data_IF.columns = data_gene.columns
        adata_IF = ad.AnnData(data_IF.T)
        adata_IF.obs = adata_IF.obs.rename_axis("cell", axis="index")
        adata_IF.obs[obs_name] = adata_IF.obs.index
        
    else:
                                
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
    
    logging.info("Processing completed successfully.")
    
    return adata_IF
