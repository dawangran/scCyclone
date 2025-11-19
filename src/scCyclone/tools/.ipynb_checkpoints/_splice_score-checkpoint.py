# -*- coding: utf-8 -*-
"""
@File    :   _splice_score.py
@Time    :   2024/09/04 
@Author  :   Dawn
@Version :   1.0
@Desc    :   splice score
"""



import numpy as np
import anndata as ad



def splice_score(
    adata_iso: ad.AnnData, 
    event_dict: dict,
    score_column: str = 'domain_number'
) -> ad.AnnData:
    """
    Compute splice scores for each event and add the total score to adata_iso.obs.

    Parameters:
    ----------
    event_dict (dict): Dictionary where keys are event IDs and values are lists of two lists:
                       - First list: Alternative transcripts (a)
                       - Second list: Total transcripts (b)
    adata_iso (anndata.AnnData): AnnData object containing isoform expression data.
                                 Expected columns in `adata_iso.var`:
                                 - score_column: Column containing scores (e.g., domain_number).
    score_column (str): Column name in `adata_iso.var` that contains the scores for each isoform. 
                        Defaults to 'domain_number'.

    Returns:
    ----------
    adata_iso (anndata.AnnData): Updated AnnData object with the total splice score added to `adata_iso.obs`.
    """
    # Check if required columns exist
    if score_column not in adata_iso.var:
        raise ValueError(f"adata_iso.var must contain the column '{score_column}'.")

    # Initialize score array
    score = np.zeros(adata_iso.n_obs)

    # Get valid isoforms in adata_iso.var
    valid_isoforms = adata_iso.var.index

    # Compute splice score for each event
    for k, v in event_dict.items():
        a, b = v

        # Filter out invalid isoforms
        a_valid = [iso for iso in a if iso in valid_isoforms]
        b_valid = [iso for iso in b if iso in valid_isoforms]

        if not a_valid or not b_valid:
            # Skip events with no valid isoforms
            continue

        a_score = np.array(adata_iso.var.loc[a_valid][score_column])
        b_score = np.array(adata_iso.var.loc[b_valid][score_column])

        # Compute dot products
        a_dot = np.dot(adata_iso[:, a_valid].X, a_score)
        b_dot = np.dot(adata_iso[:, b_valid].X, b_score)

        # Update score
        score += a_dot / b_dot

    # Add splice score to adata_iso.obs
    adata_iso.obs['splice_score'] = score

    return adata_iso