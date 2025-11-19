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
from typing import Optional






def generate_Iso_adata(
    data_path: str,
    ):
    """
    Process Cyclone isoform data and save as AnnData objects.
    
    Parameters:
    ----------
    
    data_path (str): Path to the isoform count data CSV file (isoquant result).
        
    Returns:
    ----------
    
        adata (anndata.AnnData): Processed AnnData object.
    """
    
    # Check if the data file exists
    if not os.path.exists(data_path):
        logging.error(f"Data file {data_path} does not exist.")
        return None

  
    # Read isoform count data
    logging.info(f"Reading isoform count data from {data_path}")
    chunk_size = 10000  # 每块大小
    chunks = pd.read_csv(data_path, index_col=0, sep="\t", chunksize=chunk_size)

    # 合并分块数据
    data = pd.concat(chunks, axis=0)
    
    # Create AnnData object
    adata = ad.AnnData(data.T)
    adata.obs['batch'] = [i.rsplit("_", 1)[0] for i in adata.obs_names]    
    adata.var['isoform'] = adata.var.index
   
    
    # # Remove entries with NaN gene names

    # adata = adata[:, adata.var.dropna(subset=['gene_id']).index]
    
    logging.info("Processing completed successfully.")
    
    return adata
    




def generate_PSI_adata(adata: ad.AnnData, event_info_path: str) -> Optional[ad.AnnData]:
    """
    Read Cyclone PSI (Percent Spliced In) data and save as AnnData object.

    Parameters:
    ----------
    adata (anndata.AnnData): Input AnnData object (iso adata).
    event_info_path (str): Path to the event information file (suppa2 result).

    Returns:
    ----------
    adata_as (anndata.AnnData): Processed AnnData object.
    """
    # Check if the event info file exists
    if not os.path.exists(event_info_path):
        logging.error(f"Event info file {event_info_path} does not exist.")
        return None

    # Read event information
    logging.info(f"Reading event information from {event_info_path}")
    event_info = pd.read_csv(event_info_path, sep="\t")
    event_info['type'] = event_info['event_id'].str.split(";", expand=True)[1].str.split(":", expand=True)[0]
    event_info = event_info[['event_id', 'gene_id', 'type', 'alternative_transcripts', 'total_transcripts']]

    # Extract cell names from adata
    cell_names = adata.obs.index.to_list()

    # Initialize PSI DataFrame
    logging.info("Initializing PSI DataFrame")
    psi_df = pd.DataFrame(index=cell_names, columns=event_info['event_id'])

    # Vectorized computation of PSI values
    logging.info("Computing PSI values")
    for _, row in event_info.iterrows():
        event_id = row['event_id']
        alternative_transcripts = row['alternative_transcripts'].split(",")
        total_transcripts = row['total_transcripts'].split(",")

        # Extract relevant columns from adata
        alternative_counts = adata[:, adata.var_names.isin(alternative_transcripts)].X.sum(axis=1)
        total_counts = adata[:, adata.var_names.isin(total_transcripts)].X.sum(axis=1)

        # Compute PSI values
        psi_df[event_id] = alternative_counts / total_counts

    # Create AnnData object
    logging.info("Creating AnnData object")
    adata_as = ad.AnnData(psi_df)
    adata_as.obs = adata.obs  # Copy cell annotations
    adata_as.var = event_info.set_index('event_id')  # Add event annotations

    logging.info("Processing completed successfully.")
    return adata_as

        

def generate_Gene_adata(
    adata: ad.AnnData,
    var_name: str ='gene_name'
    ):
    """
    Aggregate isoform expression data to gene level and save as AnnData object.

    Parameters:
    ----------
    
    adata (anndata.AnnData): Input AnnData object containing isoform expression data.
    var_name (str): Name of the column in `adata.var` that contains.

    Returns:
    ----------
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
    var_name: str = 'gene_name',
    bulk: bool = False,
    obs_name: Union[str, None] = None
) -> ad.AnnData:
    """
    Generate an AnnData object with isoform fraction data.

    Parameters:
    ----------
    adata (anndata.AnnData): Input AnnData object.
    var_name (str, optional): Name of the column in `adata.var` that contains gene names. Defaults to 'gene_name'.
    bulk (bool): If True, compute bulk-level isoform fractions using `obs_name`.
    obs_name (str, optional): Name of the column in `adata.obs` to group by for bulk-level computation. Defaults to None.

    Returns:
    ----------
    adata_IF (anndata.AnnData): AnnData object with isoform fraction data.
    """
    # Check if required columns exist
    if var_name not in adata.var:
        raise ValueError(f"Column '{var_name}' not found in adata.var.")
    
    if bulk:
        if obs_name not in adata.obs:
            raise ValueError(f"Column '{obs_name}' not found in adata.obs.")
        
        # Group by obs_name and compute bulk-level isoform fractions
        logging.info("Computing bulk-level isoform fractions")
        data = adata.to_df()
        data[obs_name] = adata.obs[obs_name].to_list()
        data_iso = data.groupby(obs_name).sum().T
        data_iso[var_name] = adata.var[var_name].to_list()
        
        # Compute gene-level sums
        data_gene = data_iso.groupby(var_name).sum()
        
        # Broadcast gene-level sums to match isoform-level data
        gene_sums = data_gene.reindex(data_iso[var_name]).values
        
        # Compute isoform fractions
        data_IF = data_iso.iloc[:, :-1].div(gene_sums, axis=0)
        
        # Create AnnData object
        adata_IF = ad.AnnData(data_IF.T)
        adata_IF.obs = adata_IF.obs.rename_axis("cell", axis="index")
        adata_IF.obs[obs_name] = adata_IF.obs.index
    else:
        # Compute single-cell isoform fractions
        logging.info("Computing single-cell isoform fractions")
        data = adata.to_df().T
        data[var_name] = adata.var[var_name].to_list()
        
        # Compute gene-level sums
        data_gene = data.groupby(var_name).sum()
        
        # Broadcast gene-level sums to match isoform-level data
        gene_sums = data_gene.reindex(data[var_name]).values
        
        # Compute isoform fractions
        data_IF = data.iloc[:, :-1].div(gene_sums, axis=0)
        
        # Create AnnData object
        adata_IF = ad.AnnData(data_IF.T)
        adata_IF.obs = adata.obs
    
    # Copy variable annotations from input AnnData
    adata_IF.var = adata.var
    
    logging.info("Processing completed successfully.")
    return adata_IF



def generate_APA_adata(
    adata: ad.AnnData,
    gene_name: str = 'gene_name',
    polyA_dist_name: str = "polyA_dist"
):
    """
    Aggregate isoform expression data to APA (Alternative Polyadenylation) level and save as AnnData object.

    Parameters:
    ----------
    adata (anndata.AnnData): Input AnnData object containing isoform expression data.
    gene_name (str): Name of the column in `adata.var` that contains gene names.
    polyA_dist_name (str): Name of the column in `adata.var` that contains polyA distances.

    Returns:
    ----------
    adata_apa (anndata.AnnData): AnnData object containing aggregated APA expression data.
    """

    # Check if the specified columns exist in adata.var
    if gene_name not in adata.var:
        raise ValueError(f"Column '{gene_name}' not found in adata.var.")
    if polyA_dist_name not in adata.var:
        raise ValueError(f"Column '{polyA_dist_name}' not found in adata.var.")

    # Filter out rows where polyA_dist_name has missing values
    adata = adata[:, adata.var[polyA_dist_name].dropna().index]

    # Create a new column 'apa_id' in the format 'gene_name_polyA_dist'
    adata.var['apa_id'] = ["{}_{}".format(i, abs(int(j))) for i, j in zip(adata.var[gene_name], adata.var[polyA_dist_name])]

    # Convert the AnnData object to a DataFrame and add the 'apa_id' column
    data_iso = adata.to_df().T  # Transpose to make each isoform's expression data a row
    data_iso['apa_id'] = list(adata.var['apa_id'])  # Add the 'apa_id' column

    # Group by 'apa_id' and aggregate the data by summing up the values
    data_iso = data_iso.groupby("apa_id").agg(sum).T  # Aggregate and transpose back to the original format

    # Convert the aggregated data back to an AnnData object
    adata_apa = ad.AnnData(data_iso)

    # Extract gene names and lengths from the 'apa_id' and add them to the new AnnData object's var
    adata_apa.var[gene_name] = [i.rsplit("_",1)[0] for i in adata_apa.var.index]  # Extract gene names
    adata_apa.var['length'] = [int(i.rsplit("_",1)[1]) for i in adata_apa.var.index]  # Extract lengths
    adata_apa.var['apa'] = list(adata_apa.var.index)

    # Rank the lengths within each gene group in descending order and add the ranks to the var
    adata_apa.var['rank'] = adata_apa.var.groupby(gene_name)['length'].rank(method='dense', ascending=False).astype(int)

    # Sort the var by 'gene_name' and 'rank' to organize the data
    adata_apa.var = adata_apa.var.sort_values(by=[gene_name, 'rank'])
    
    # Filter out genes with only one APA site (since APA analysis requires multiple sites)
    adata_apa=adata_apa[:,adata_apa.var[gene_name].isin(adata_apa.var[adata_apa.var['rank']>1][gene_name].to_list())]
    
    
    # Compute apa fractions within each gene
    data = adata_apa.to_df().T
    data[gene_name] = adata_apa.var[gene_name].to_list()
    data_gene = data.groupby(gene_name).sum()
    gene_sums = data_gene.reindex(data[gene_name]).values
    data_apai = data.iloc[:, :-1].div(gene_sums, axis=0)
    adata_apa_IF = ad.AnnData(data_apai.T)
    adata_apa_IF.obs = adata.obs
    adata_apa_IF.var = adata_apa.var
    
    logging.info("Processing completed successfully.")
    
    return adata_apa_IF



# def generate_IF_adata(
#     adata: ad.AnnData, 
#     var_name: str ='gene_name',
#     bulk: bool =False,
#     obs_name: Union[str, None] = None
#     ):
#     """
#     Generate an AnnData object with isoform fraction data.

#     Parameters:
#     ----------

#     adata (anndata.AnnData): Input AnnData object.
#     gene_name (str, optional): Name of the column in `adata.var` that contains. Defaults to 'gene_name'.
#     bulk (bool): IF bulk is true, the obs_name command is used to merge the obs_name command to generate the IF of the bulk level.
#     obs_name (str): Name of the column in `adata.obs` that contains. Defaults to  None.
        
#     Returns:
#     ----------
#     adata_IF (anndata.AnnData): AnnData object with isoform fraction data.
#     """
    
    
#     # Check if gene_name column exists in adata.var
#     if var_name not in adata.var:
#         raise ValueError(f"Column '{var_name}' not found in adata.var.")
    
#     if bulk:
#         if obs_name not in adata.obs:
#             raise ValueError(f"Column '{obs_name}' not found in adata.obs.")
        
#         data = adata.to_df()
#         data[obs_name] = adata.obs[obs_name].to_list()
#         data_iso = data.groupby(obs_name).agg("sum")
#         data_iso = data_iso.T
#         data_iso[var_name] = adata.var[var_name].to_list()
#         data_gene = data_iso.iloc[:, :-1]
#         data_gene[var_name] = adata.var[var_name].to_list()
#         data_gene = data_gene.groupby(var_name).agg("sum")
#         data_IF = pd.DataFrame()
#         for i, j in enumerate(data_iso.index):
#             sub_data = data_iso.iloc[i, :-1].values
#             gene = data_iso.iloc[i, -1]
#             sub_gene = data_gene.loc[gene].values
#             data_IF[j] = list(sub_data / sub_gene)
            
#             if i % 1000 == 0:
#                 print("process succes {}".format(i))
                
#         data_IF = data_IF.T
#         data_IF.columns = data_gene.columns
#         adata_IF = ad.AnnData(data_IF.T)
#         adata_IF.obs = adata_IF.obs.rename_axis("cell", axis="index")
#         adata_IF.obs[obs_name] = adata_IF.obs.index
        
#     else:
                                
#         # Convert AnnData to DataFrame
#         data = adata.to_df().T
#         data[var_name] = adata.var[var_name].to_list()


#         # Group data by gene and aggregate by sum
#         data_gene = data.groupby(var_name).agg(sum)

#         # Initialize DataFrame for isoform fraction data
#         data_IF = pd.DataFrame()

#         # Calculate isoform fraction for each gene
#         for i, j in enumerate(data.index):
#             sub_data = data.iloc[i, :-1].values
#             gene = data.iloc[i, -1]
#             sub_gene = data_gene.loc[gene].values
#             data_IF[j] = list(sub_data / sub_gene)
            
#             if i % 1000 == 0:
#                 print("Process successful for {}".format(i))

#         # Transpose DataFrame and create AnnData object
#         data_IF = data_IF.T
#         data_IF.columns = data_gene.columns
#         adata_IF = ad.AnnData(data_IF.T)
#         adata_IF.obs = adata.obs
        
#     # Replace NaN values with 0
#     adata_IF.X = np.nan_to_num(adata_IF.X)
    
#     # Copy variable annotations from input AnnData
#     adata_IF.var = adata.var
    
#     logging.info("Processing completed successfully.")
    
#     return adata_IF


# def generate_PSI_adata(
#     adata: ad.AnnData,
#     event_info_path: str
#     ):
#     """
#     Read Cyclone PSI (Percent Spliced In) data and save as AnnData object.

#     Parameters:
#     ----------
    
#     adata (anndata.AnnData): Input AnnData object(iso adata).
#     event_info_path (str): Path to the event information file(suppa2 result).

#     Returns:
#     ----------
#         adata_as (anndata.AnnData): Processed AnnData object.
#     """

    
#     # Check if the event info file exists
#     if not os.path.exists(event_info_path):
#         logging.error(f"Transcript info file {event_info_path} does not exist.")
#         return None

#     # Read event information
#     data = pd.read_csv(event_info_path, sep="\t")
#     data['type'] = data['event_id'].str.split(";", expand=True)[1].str.split(":", expand=True)[0]
#     data = data[['event_id', 'gene_id', 'type', 'alternative_transcripts', 'total_transcripts']]

#     # Prepare PSI values
#     cell = adata.obs.index.to_list()
#     feature = data['event_id'].to_list()
#     data["acceptor"] = data['alternative_transcripts'].str.split(",")
#     data['total'] = data['total_transcripts'].str.split(",")   
#     df = pd.DataFrame(index=cell, columns=feature)
#     data_list = []

#     for i in range(0, df.shape[1], 1000):
#         start_col = i
#         end_col = i + 1000 if i + 1000 <= df.shape[1] else df.shape[1]
#         df_subset = df.iloc[:, start_col:end_col]

#         for column in df_subset.columns:
#             alternative_transcripts = data[data['event_id'] == column]['acceptor'].values[0]
#             total_transcripts = data[data['event_id'] == column]['total'].values[0]
            
#             alternative_transcripts_number = adata[:, adata.var_names.isin(alternative_transcripts)].X.sum(axis=1)
#             total_transcripts_number = adata[:, adata.var_names.isin(total_transcripts)].X.sum(axis=1)
#             df_subset[column] = alternative_transcripts_number / total_transcripts_number

#         data_list.append(df_subset)

#     data_as = pd.concat(data_list, axis=1)
#     adata_as = ad.AnnData(data_as)
#     adata_as.obs = adata.obs
#     adata_as.var['gene_id'] = data['gene_id'].to_list()
#     adata_as.var['type'] = data['type'].to_list()
#     adata_as.var['alternative_transcripts'] = data['alternative_transcripts'].to_list()
#     adata_as.var['total_transcripts'] = data['total_transcripts'].to_list()
    
#     logging.info("Processing completed successfully.")
    
#     return adata_as



# def generate_APA_adata(
#     adata: ad.AnnData,
#     gene_name: str = 'gene_name',
#     polyA_dist_name: str = "polyA_dist"
# ):
#     """
#     Aggregate isoform expression data to APA (Alternative Polyadenylation) level and save as AnnData object.

#     Parameters:
#     ----------
#     adata (anndata.AnnData): Input AnnData object containing isoform expression data.
#     gene_name (str): Name of the column in `adata.var` that contains gene names.
#     polyA_dist_name (str): Name of the column in `adata.var` that contains polyA distances.

#     Returns:
#     ----------
#     adata_apa (anndata.AnnData): AnnData object containing aggregated APA expression data.
#     """

#     # Check if the specified columns exist in adata.var
#     if gene_name not in adata.var:
#         raise ValueError(f"Column '{gene_name}' not found in adata.var.")
#     if polyA_dist_name not in adata.var:
#         raise ValueError(f"Column '{polyA_dist_name}' not found in adata.var.")

#     # Filter out rows where polyA_dist_name has missing values
#     adata = adata[:, adata.var[polyA_dist_name].dropna().index]

#     # Create a new column 'apa_id' in the format 'gene_name_polyA_dist'
#     adata.var['apa_id'] = ["{}_{}".format(i, abs(int(j))) for i, j in zip(adata.var[gene_name], adata.var[polyA_dist_name])]

#     # Convert the AnnData object to a DataFrame and add the 'apa_id' column
#     data_iso = adata.to_df().T  # Transpose to make each isoform's expression data a row
#     data_iso['apa_id'] = list(adata.var['apa_id'])  # Add the 'apa_id' column

#     # Group by 'apa_id' and aggregate the data by summing up the values
#     data_iso = data_iso.groupby("apa_id").agg(sum).T  # Aggregate and transpose back to the original format

#     # Convert the aggregated data back to an AnnData object
#     adata_apa = ad.AnnData(data_iso)

#     # Extract gene names and lengths from the 'apa_id' and add them to the new AnnData object's var
#     adata_apa.var[gene_name] = [i.split("_")[0] for i in adata_apa.var.index]  # Extract gene names
#     adata_apa.var['length'] = [i.split("_")[1] for i in adata_apa.var.index]  # Extract lengths

#     # Rank the lengths within each gene group in descending order and add the ranks to the var
#     adata_apa.var['rank'] = adata_apa.var.groupby(gene_name)['length'].rank(method='dense', ascending=False).astype(int)

#     # Sort the var by 'gene_name' and 'rank' to organize the data
#     adata_apa.var = adata_apa.var.sort_values(by=[gene_name, 'rank'])
    
#     logging.info("Processing completed successfully.")

#     return adata_apa