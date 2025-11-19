# -*- coding: utf-8 -*-
"""
@File    :   _rank_ifs_groups.py
@Time    :   2024/12/23
@Author  :   Dawn
@Version :   1.1
@Desc    :   DIF for scCyclone
             - Stable and reproducible permutation test (two-tailed)
             - Expression ratio denominator: if no NA → use total rows; if NA → use non-NA rows
             - var_name is written according to parameter, not hard-coded
"""

from typing import Union, List, Dict, Any
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import anndata as ad
from joblib import Parallel, delayed

# Kept for API compatibility (although not directly used here)
import scanpy as sc

from . import _utils


# ------------------------------
# Utility functions (local fallback for multiple testing correction)
# ------------------------------
def _bonferroni_safe(pvals: List[float]) -> List[float]:
    """If _utils.compute_pvalue_bonferroni is unavailable, use local Bonferroni."""
    n = len(pvals)
    if n == 0:
        return []
    arr = np.asarray([1.0 if (p is None or not np.isfinite(p)) else float(p) for p in pvals], dtype=float)
    adj = np.minimum(arr * n, 1.0)
    return adj.tolist()


def _two_sided_pvalue(dif_shuffle: List[float], dif_observed: float) -> float:
    """
    Two-tailed permutation p-value.  
    Returns 1.0 if no permutations or invalid observation.  
    No rounding here to avoid collapsing small effects into ties.
    """
    if not np.isfinite(dif_observed) or len(dif_shuffle) == 0:
        return 1.0
    null = np.asarray(dif_shuffle, dtype=float)
    k = int((np.abs(null) >= abs(dif_observed)).sum())
    n = null.size
    # Classical proportion; to be more conservative use (k+1)/(n+1)
    return k / n


def _compute_dif(
    e: str,
    data_ref: pd.DataFrame,
    data_tgt: pd.DataFrame,
    valid_cells: int,
    n_bins: int,
    rng: np.random.Generator,
) -> Dict[str, Any]:
    """
    Compute DIF statistics for a single isoform.

    Conventions:
    - Denominator for expression ratio (rpr/tpr):
        If the column contains no NA → use total rows;
        If the column contains NA → use non-NA rows.
    - Permutation is done by shuffling indices (no per-loop seeding to avoid dependency).
    - Permutation performed only if both groups have ≥ valid_cells and n_bins > 0.
    """
    # Drop NA for mean and numerator calculation
    sub_ref = data_ref[[e]].dropna()
    sub_tgt = data_tgt[[e]].dropna()

    # Denominator rule: no NA → full rows; otherwise → non-NA rows
    ref_den = data_ref.shape[0] if data_ref[e].isna().sum() == 0 else max(sub_ref.shape[0], 1)
    tgt_den = data_tgt.shape[0] if data_tgt[e].isna().sum() == 0 else max(sub_tgt.shape[0], 1)

    # Expression ratio (>0 cell proportion)
    rpr_value = float(sub_ref[e].gt(0).sum()) / ref_den if ref_den > 0 else 0.0
    tpr_value = float(sub_tgt[e].gt(0).sum()) / tgt_den if tgt_den > 0 else 0.0
    dpr_value = tpr_value - rpr_value

    # Group means (sub_* already drop NA)
    rif_value = float(sub_ref[e].mean()) if sub_ref.shape[0] else np.nan
    tif_value = float(sub_tgt[e].mean()) if sub_tgt.shape[0] else np.nan
    dif_observed = float(tif_value - rif_value) if np.isfinite(rif_value) and np.isfinite(tif_value) else np.nan

    # Permutation
    dif_shuffle: List[float] = []
    if (
        sub_ref.shape[0] >= valid_cells
        and sub_tgt.shape[0] >= valid_cells
        and n_bins > 0
        and np.isfinite(dif_observed)
    ):
        # Concatenate values; each permutation shuffles indices
        vals = pd.concat([sub_ref[e], sub_tgt[e]], axis=0).to_numpy()
        n_ref = sub_ref.shape[0]
        n_tot = vals.shape[0]

        for _ in range(n_bins):
            perm = rng.permutation(n_tot)
            ref_idx = perm[:n_ref]
            tgt_idx = perm[n_ref:]
            dif_shuffle.append(float(vals[tgt_idx].mean() - vals[ref_idx].mean()))

    return {
        "iso": e,
        "dif_observed": dif_observed,
        "dif_shuffle": dif_shuffle,
        "rif_value": rif_value,
        "tif_value": tif_value,
        "rpr_value": rpr_value,
        "tpr_value": tpr_value,
        "dpr_value": dpr_value,
    }


def rank_ifs_groups(
    adata: ad.AnnData,
    groupby: str,
    groups: Union[str, List[str]] = "all",
    reference: str = "rest",
    key_added: Union[str, None] = None,
    valid_cells: int = 50,
    n_bins: int = 100,
    random_seed: int = 22,
    var_name: str = "gene_name",
) -> ad.AnnData:
    """
    Compute DIF across groups and rank isoforms.

    Parameters
    ----------
    adata : AnnData
    groupby : str
        Column name in adata.obs for grouping.
    groups : list or "all"
        Target groups; "all" means all groups.
    reference : str
        Comparison reference: 'rest' means merge all remaining groups,
        or specify an explicit group name.
    key_added : str or None
        Where to store results in adata.uns[key_added].
    valid_cells : int
        Minimum number of valid (non-NA) cells per group to run permutation.
    n_bins : int
        Number of permutations per isoform; 0 = no permutation (p=1.0).
    random_seed : int
        Random seed for reproducibility.
    var_name : str
        Column name from adata.var to attach as annotation.

    Returns
    -------
    AnnData (results added to adata.uns[key_added])
    """
    # Order of groups (utility from your helper functions)
    groups_order = _utils.check_groups(adata, groupby, groups, reference)
    if key_added is None:
        key_added = "rank_ifs_groups"

    adata.uns[key_added] = {
        "params": {
            "groupby": groupby,
            "reference": reference,
            "valid_cells": valid_cells,
            "n_bins": n_bins,
            "random_seed": random_seed,
            "var_name": var_name,
        }
    }

    iso_index = list(adata.var.index)

    data_iso_dict: Dict[str, List[str]] = {}
    data_dif_observed_dict: Dict[str, List[float]] = {}
    data_pval_dict: Dict[str, List[float]] = {}
    data_pval_adj_dict: Dict[str, List[float]] = {}
    data_dpr_dict: Dict[str, List[float]] = {}
    data_rpr_dict: Dict[str, List[float]] = {}
    data_tpr_dict: Dict[str, List[float]] = {}
    data_rif_dict: Dict[str, List[float]] = {}
    data_tif_dict: Dict[str, List[float]] = {}
    var_name_dict: Dict[str, List[Any]] = {}

    # Master RNG (used to generate per-isoform sub-seeds; ensures parallel reproducibility)
    master_rng = np.random.default_rng(random_seed)

    for g in groups_order:
        if g == reference:
            continue

        print(f"Group {g} start!")
        adata_tgt = adata[adata.obs[groupby] == g]
        if reference == "rest":
            adata_ref = adata[adata.obs[groupby] != g]
        else:
            adata_ref = adata[adata.obs[groupby].isin([reference])]

        df_tgt = adata_tgt.to_df()
        df_ref = adata_ref.to_df()

        # Per-isoform seeds
        seeds = master_rng.integers(0, np.iinfo(np.int32).max, size=len(iso_index))

        # Parallel DIF computation
        print("Compute DIF...")
        results = Parallel(n_jobs=-1, prefer="threads")(
            delayed(_compute_dif)(
                e,
                df_ref,
                df_tgt,
                valid_cells,
                n_bins,
                np.random.default_rng(int(seeds[i])),
            )
            for i, e in enumerate(iso_index)
        )
        res = pd.DataFrame(results)

        # Two-sided p-values
        print("Compute p-values...")
        res["pvals"] = [_two_sided_pvalue(dsh, dob) for dsh, dob in zip(res["dif_shuffle"], res["dif_observed"])]

        # Ranking: |effect| desc, p-value asc
        res["_abs_dif"] = res["dif_observed"].abs().fillna(-np.inf)
        res = res.sort_values(by=["_abs_dif", "pvals"], ascending=[False, True]).drop(columns=["_abs_dif"])

        # Collect results (rounding only for display)
        data_iso_dict[g] = res["iso"].tolist()
        data_dif_observed_dict[g] = [None if pd.isna(x) else float(np.round(x, 3)) for x in res["dif_observed"]]
        data_rif_dict[g] = [None if pd.isna(x) else float(np.round(x, 3)) for x in res["rif_value"]]
        data_tif_dict[g] = [None if pd.isna(x) else float(np.round(x, 3)) for x in res["tif_value"]]
        data_pval_dict[g] = [1.0 if pd.isna(x) else float(x) for x in res["pvals"]]
        data_dpr_dict[g] = [None if pd.isna(x) else float(np.round(x, 4)) for x in res["dpr_value"]]
        data_rpr_dict[g] = [None if pd.isna(x) else float(np.round(x, 4)) for x in res["rpr_value"]]
        data_tpr_dict[g] = [None if pd.isna(x) else float(np.round(x, 4)) for x in res["tpr_value"]]

        # Multiple testing correction
        try:
            data_pval_adj_dict[g] = _utils.compute_pvalue_bonferroni(data_pval_dict[g])
        except Exception:
            data_pval_adj_dict[g] = _bonferroni_safe(data_pval_dict[g])

        # Attach var_name annotation
        if var_name in adata.var.columns:
            var_series = adata_tgt[:, res["iso"].tolist()].var[var_name]
            var_name_dict[g] = var_series.tolist()
        else:
            var_name_dict[g] = res["iso"].tolist()

        print(f"Group {g} complete!")
        print("-----------------------------------------")

    # Convert to recarray (scanpy-compatible)
    def _to_records(d: Dict[str, List[Any]]):
        return pd.DataFrame(d).to_records(index=False) if d else pd.DataFrame().to_records(index=False)

    adata.uns[key_added]["names"]      = _to_records(data_iso_dict)
    adata.uns[key_added]["dif"]        = _to_records(data_dif_observed_dict)
    adata.uns[key_added]["pvals"]      = _to_records(data_pval_dict)
    adata.uns[key_added]["pvals_adj"]  = _to_records(data_pval_adj_dict)
    adata.uns[key_added]["dpr"]        = _to_records(data_dpr_dict)
    adata.uns[key_added]["rpr"]        = _to_records(data_rpr_dict)
    adata.uns[key_added]["tpr"]        = _to_records(data_tpr_dict)
    adata.uns[key_added]["rif"]        = _to_records(data_rif_dict)
    adata.uns[key_added]["tif"]        = _to_records(data_tif_dict)
    # Stored under user-specified var_name instead of hard-coded "gene_name"
    adata.uns[key_added][var_name]     = _to_records(var_name_dict)

    return adata
