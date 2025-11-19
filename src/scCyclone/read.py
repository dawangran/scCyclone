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
import logging
from scipy.io import mmread
from scipy.sparse import csr_matrix


def generate_Iso_adata(
    data_path: str,
    mode: str,
    matrix_prefix: str = "OUT.discovered_transcript_grouped_counts",
    chunk_size: int = 10000,
):
    """
    Process Cyclone isoform data and return an AnnData object.

    Parameters
    ----------
    data_path : str
        - If `mode="csv"`: path to the isoform count TSV/CSV file (isoquant result),
          where rows = isoforms and columns = cells (or vice versa; this function
          assumes rows = isoforms and transposes).
        - If `mode="mtx"`: path to a directory containing:
            {matrix_prefix}.matrix.mtx
            {matrix_prefix}.barcodes.tsv
            {matrix_prefix}.features.tsv

    mode : {"csv", "mtx"}
        - "csv" : read as a dense TSV/CSV count matrix.
        - "mtx" : read as 10x/isoquant-style matrix.mtx + barcodes + features.

    matrix_prefix : str, default "OUT.discovered_transcript_grouped_counts"
        Prefix used for the mtx / barcodes / features files in `mode="mtx"`.

    chunk_size : int, default 10000
        Chunk size for reading large TSV/CSV files in `mode="csv"`.

    Returns
    -------
    adata : anndata.AnnData or None
        Processed AnnData object. Returns None if any required file is missing
        or if an error is encountered.
    """

    if mode not in {"csv", "mtx"}:
        raise ValueError(f"mode must be 'csv' or 'mtx', got '{mode}'")

    # ===================== CSV / TSV MODE =====================
    if mode == "csv":
        if not os.path.exists(data_path):
            logging.error(f"Data file {data_path} does not exist.")
            return None

        logging.info(f"Reading isoform count data (CSV/TSV) from {data_path}")

        # Read in chunks to avoid loading huge matrices into memory at once
        chunks = pd.read_csv(
            data_path,
            index_col=0,
            sep="\t",
            chunksize=chunk_size
        )
        data = pd.concat(chunks, axis=0)

        # AnnData expects cells in rows, features in columns
        # Here we assume rows = isoforms, so we transpose
        adata = ad.AnnData(data.T)

        # Store isoform IDs in .var column
        adata.var["isoform"] = adata.var.index.astype(str)

        logging.info("Processing in CSV/TSV mode completed successfully.")
        return adata

    # ===================== MTX MODE =====================
    if mode == "mtx":
        if not os.path.isdir(data_path):
            logging.error(f"Directory {data_path} does not exist.")
            return None

        matrix_path = os.path.join(
            data_path, f"{matrix_prefix}.matrix.mtx"
        )
        barcodes_path = os.path.join(
            data_path, f"{matrix_prefix}.barcodes.tsv"
        )
        features_path = os.path.join(
            data_path, f"{matrix_prefix}.features.tsv"
        )

        # Check that all required files exist
        for fp in [matrix_path, barcodes_path, features_path]:
            if not os.path.exists(fp):
                logging.error(f"Required file {fp} does not exist.")
                return None

        logging.info(f"Reading sparse matrix from {matrix_path}")
        # mmread returns (features × cells) by default for 10x-style matrices
        # We transpose to get shape (cells × features)
        mtx = csr_matrix(mmread(matrix_path).T)

        logging.info(f"Reading barcodes from {barcodes_path}")
        # barcodes.tsv: one barcode per row, use as index
        obs = pd.read_csv(barcodes_path, header=None, index_col=0, sep="\t")
        obs.index = obs.index.astype(str)
        obs.index.name = None

        logging.info(f"Reading features from {features_path}")
        # features.tsv: one feature (isoform) per row, use as index
        var_raw = pd.read_csv(features_path, header=None, index_col=0, sep="\t")
        var_raw.index = var_raw.index.astype(str)

        # *** KEY CHANGE: discard extra columns, keep only index ***
        var = pd.DataFrame(index=var_raw.index)
        var.index.name = None

        # Create AnnData object from sparse matrix and metadata
        adata = ad.AnnData(mtx, obs=obs, var=var)

        # If isoform column is missing in .var, create it from the index
        if "isoform" not in adata.var.columns:
            adata.var["isoform"] = adata.var.index

        logging.info("Processing in MTX mode completed successfully.")
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


