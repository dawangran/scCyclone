# -*- coding: utf-8 -*-
"""
@File    :   _rank_psis_groups.py
@Time    :   2024/09/01
@Author  :   Dawn
@Version :   1.3
@Desc    :   DPSI for scCyclone
            - Keep zeros, drop NaN-only (filter by non-NaN fraction)
            - Actually use filtered events
            - Fast p-values: approx / numba-permutation / python-permutation
            - Numba RNG fix (no RandomState; deterministic Fisher–Yates)
            - Robust RNG (no per-iter reseeding), safer outputs
"""

import warnings
warnings.filterwarnings('ignore')

from typing import Union, List, Dict, Any

import math
import numpy as np
import pandas as pd
import anndata as ad
from joblib import Parallel, delayed

# optional: numba for fast permutation; will gracefully fall back if unavailable
try:
    import numba as nb
    _HAS_NUMBA = True
except Exception:
    _HAS_NUMBA = False

# external utils in the same package (must provide: check_groups, compute_pvalue_bonferroni)
from . import _utils


# =============================================================================
# Utilities
# =============================================================================

def _dict_to_rec(d: Dict[str, List[Any]]) -> np.recarray:
    """Convert dict of lists to a structured recarray (safe on empty/misaligned)."""
    if not d:
        return np.rec.array([])
    lengths = [len(v) for v in d.values()]
    if len(set(lengths)) > 1:
        min_len = min(lengths)
        d = {k: v[:min_len] for k, v in d.items()}
    return pd.DataFrame(d).to_records(index=False)


def _filter_event(
    adata: ad.AnnData,
    groupby: str,
    groups: List[str],
    percent: float,
) -> List[str]:
    """
    Filter events by requiring (in at least one group in `groups`):
    fraction of non-NaN cells for the event >= `percent`.

    Notes
    -----
    - 0 values are meaningful and kept.
    - Only NaN are considered missing.
    """
    df = adata.to_df()  # dense (cells x events), preserves NaN
    keep: set = set()

    for g in groups:
        mask = (adata.obs[groupby] == g).values
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue
        sub = df.loc[mask]
        valid_frac = 1.0 - sub.isna().sum(axis=0) / float(n_cells)
        cols = valid_frac[valid_frac >= percent].index
        keep.update(cols.tolist())

    return list(keep)


def _obs_stats_median_mean(x_ref: np.ndarray, x_tgt: np.ndarray):
    """
    Return observed stats based on medians plus valid-cell proportions.

    Returns
    -------
    dpsi_obs, (rpsi, tpsi, rpr, tpr, dpr)
    """
    r_valid = ~np.isnan(x_ref)
    t_valid = ~np.isnan(x_tgt)
    r = x_ref[r_valid]
    t = x_tgt[t_valid]

    rpr = r.size / x_ref.size if x_ref.size else np.nan
    tpr = t.size / x_tgt.size if x_tgt.size else np.nan
    dpr = (tpr - rpr) if (np.isfinite(rpr) and np.isfinite(tpr)) else np.nan

    if r.size == 0 or t.size == 0:
        return np.nan, (np.nan, np.nan, rpr, tpr, dpr)

    rpsi = float(np.round(np.median(r), 3))
    tpsi = float(np.round(np.median(t), 3))
    dpsi_obs = float(np.round(tpsi - rpsi, 3))
    return dpsi_obs, (rpsi, tpsi, rpr, tpr, dpr)


# =============================================================================
# Permutation engines
# =============================================================================

def _perm_pval_py(x_ref: np.ndarray, x_tgt: np.ndarray, n_bins: int, seed: int) -> float:
    """
    Pure Python permutation on pooled valid cells.
    Statistic = mean(target) - mean(ref) (kept for compatibility with older code).
    P-value is one-tailed sign-aware by comparing against observed median-based sign.
    """
    r_valid = x_ref[~np.isnan(x_ref)]
    t_valid = x_tgt[~np.isnan(x_tgt)]
    n_r = r_valid.size
    n_t = t_valid.size
    if n_r == 0 or n_t == 0 or n_bins <= 1:
        return 1.0

    obs = float(np.round(np.median(t_valid) - np.median(r_valid), 3))
    x = np.concatenate([r_valid, t_valid], axis=0)
    n = x.size
    k = n_r
    rng = np.random.default_rng(seed)

    count = 0
    for _ in range(n_bins):
        idx_r = rng.choice(n, size=k, replace=False)
        mask_r = np.zeros(n, dtype=bool)
        mask_r[idx_r] = True
        mr = x[mask_r].mean()
        mt = x[~mask_r].mean()
        stat = mt - mr
        if obs > 0 and stat > obs:
            count += 1
        elif obs < 0 and stat < obs:
            count += 1
        elif obs == 0 and stat != 0:
            count += 1

    return count / n_bins


# ----- Numba-compatible engine (deterministic Fisher–Yates; no RandomState) -----
if _HAS_NUMBA:
    @nb.njit(cache=True, fastmath=True)
    def _fy_shuffle_inplace(idx, seed):
        """In-place Fisher–Yates shuffle using Numba-friendly RNG."""
        np.random.seed(seed)
        n = idx.shape[0]
        for i in range(n - 1, 0, -1):
            j = np.random.randint(0, i + 1)  # [0, i]
            tmp = idx[i]
            idx[i] = idx[j]
            idx[j] = tmp

    @nb.njit(cache=True, fastmath=True)
    def _perm_count_numba(x, n_r, n_bins, seed, obs):
        """
        x: 1D float64 array (no NaN; concatenated ref+target)
        n_r: ref count
        statistic: mean(target) - mean(ref)
        one-tailed, sign-aware (match legacy behavior)
        """
        n = x.shape[0]
        idx = np.arange(n)
        count = 0

        for it in range(n_bins):
            # ensure determinism: vary seed each iteration
            _fy_shuffle_inplace(idx, seed + it)

            # first n_r -> ref, rest -> target
            s_r = 0.0
            for j in range(n_r):
                s_r += x[idx[j]]
            s_t = 0.0
            for j in range(n_r, n):
                s_t += x[idx[j]]

            mr = s_r / n_r
            mt = s_t / (n - n_r)
            stat = mt - mr

            if obs > 0.0 and stat > obs:
                count += 1
            elif obs < 0.0 and stat < obs:
                count += 1
            elif obs == 0.0 and stat != 0.0:
                count += 1

        return count

    def _perm_pval_numba(x_ref: np.ndarray, x_tgt: np.ndarray, n_bins: int, seed: int) -> float:
        r_valid = x_ref[~np.isnan(x_ref)]
        t_valid = x_tgt[~np.isnan(x_tgt)]
        n_r = r_valid.size
        n_t = t_valid.size
        if n_r == 0 or n_t == 0 or n_bins <= 1:
            return 1.0
        obs = float(np.round(np.median(t_valid) - np.median(r_valid), 3))
        x = np.concatenate((r_valid, t_valid)).astype(np.float64)
        cnt = _perm_count_numba(x, n_r, n_bins, seed, obs)
        return cnt / n_bins


def _pval_approx_normal(x_ref: np.ndarray, x_tgt: np.ndarray, two_tailed: bool) -> float:
    """
    Fast normal approximation:
      var( mean_t - mean_r ) ≈ pooled_var * (1/n_r + 1/n_t)
    We compute p for mean-difference; observed sign from median diff can be applied if one-tailed.
    """
    r = x_ref[~np.isnan(x_ref)]
    t = x_tgt[~np.isnan(x_tgt)]
    n_r, n_t = r.size, t.size
    if n_r == 0 or n_t == 0:
        return 1.0

    # mean difference as quick stat
    d_mean = t.mean() - r.mean()
    # pooled variance
    s2_r = r.var(ddof=1) if n_r > 1 else 0.0
    s2_t = t.var(ddof=1) if n_t > 1 else 0.0
    s2 = (s2_r + s2_t) / 2.0
    var = s2 * (1.0 / n_r + 1.0 / n_t)
    if var <= 0:
        return 1.0
    z = d_mean / math.sqrt(var)

    from math import erf, sqrt
    if two_tailed:
        p = 2.0 * (1.0 - 0.5 * (1.0 + erf(abs(z) / sqrt(2.0))))
    else:
        # one-tailed
        p = 1.0 - 0.5 * (1.0 + erf(z / sqrt(2.0)))
        p = float(max(min(p, 1.0), 0.0))
    return float(p)


# =============================================================================
# Core per-event observed statistic (no permutation here)
# =============================================================================

def _compute_dpsi_observed(
    e: str,
    data_a: pd.DataFrame,
    data_b: pd.DataFrame,
) -> dict:
    """
    Compute observed dPSI and descriptive stats for one event.
    This function avoids permutations for speed; p-values are computed later.
    """
    x_ref = data_a[e].to_numpy(dtype=float, copy=False)
    x_tgt = data_b[e].to_numpy(dtype=float, copy=False)

    dpsi_observed, (rpsi_value, tpsi_value, rpr_value, tpr_value, dpr_value) = \
        _obs_stats_median_mean(x_ref, x_tgt)

    return {
        "event": e,
        "dpsi_shuffle": None,         # placeholder for backward compatibility
        "dpsi_observed": dpsi_observed,
        "rpsi_value": rpsi_value,
        "tpsi_value": tpsi_value,
        "rpr_value": rpr_value,
        "tpr_value": tpr_value,
        "dpr_value": dpr_value,
        "_x_ref": x_ref,              # keep arrays for fast p-value step
        "_x_tgt": x_tgt,
    }


# =============================================================================
# Public API
# =============================================================================

def rank_psis_groups(
    adata: ad.AnnData,
    groupby: str,
    groups: Union[str, List[str]] = "all",
    reference: str = "rest",
    key_added: Union[str, None] = None,
    percent: float = 0.1,
    valid_cells: int = 50,
    n_bins: int = 100,
    random_seed: int = 22,
    var_name: str = "gene_name",
    two_tailed: bool = False,
    method: str = "perm_numba",   # "perm_numba" | "perm_py" | "approx"
) -> ad.AnnData:
    """
    Rank PSI-like events for characterizing groups (target vs reference/rest).

    Parameters
    ----------
    adata : AnnData
    groupby : str
        Key in adata.obs for grouping.
    groups : 'all' or list[str]
        Target groups to compare; if 'all', use all observed categories.
    reference : str
        'rest' or a specific group id to be used as control.
    key_added : str or None
        adata.uns key to save results. Default: "rank_psis_groups".
    percent : float
        Minimal fraction of non-NaN cells within a group to keep an event (0 kept; NaN dropped).
    valid_cells : int
        Minimal valid cells per group to enable p-value estimation.
    n_bins : int
        Number of permutations (when a permutation method is used).
    random_seed : int
        RNG seed for permutations.
    var_name : str
        Column in adata.var to attach as "gene_name".
    two_tailed : bool
        Two-tailed p for approx, or interpreted for one-tailed when False.
    method : str
        "perm_numba" (fastest, if numba available), "perm_py" (pure Python), or "approx" (normal approx).

    Returns
    -------
    AnnData with results in adata.uns[key_added].
    """
    groups_order = _utils.check_groups(adata, groupby, groups, reference)
    print(f"Groups order: {groups_order}")

    # Filter events by NaN fraction rule (keep zeros)
    event_list = _filter_event(adata, groupby=groupby, groups=groups_order, percent=percent)
    print(f"Filter event: removed {adata.n_vars - len(event_list)} / {adata.n_vars}")

    if key_added is None:
        key_added = "rank_psis_groups"

    adata.uns[key_added] = {}
    adata.uns[key_added]["params"] = {
        "groupby": groupby,
        "reference": reference,
        "percent": percent,
        "valid_cells": valid_cells,
        "n_bins": n_bins,
        "random_seed": random_seed,
        "two_tailed": two_tailed,
        "var_name": var_name,
        "used_events": len(event_list),
        "method": method,
        "numba": _HAS_NUMBA,
    }

    # Prepare containers
    data_event_dict: Dict[str, List[str]] = {}
    data_dpsi_observed_dict: Dict[str, List[float]] = {}
    data_pval_dict: Dict[str, List[float]] = {}
    data_pval_adj_dict: Dict[str, List[float]] = {}
    data_dpr_dict: Dict[str, List[float]] = {}
    data_rpr_dict: Dict[str, List[float]] = {}
    data_tpr_dict: Dict[str, List[float]] = {}
    data_rpsi_dict: Dict[str, List[float]] = {}
    data_tpsi_dict: Dict[str, List[float]] = {}
    var_name_dict: Dict[str, List[str]] = {}

    if len(event_list) == 0:
        print("No events passed filtering; writing empty results.")
    else:
        # One-time dense dataframe (to keep NaN)
        df_all = adata.to_df()

        for i in groups_order:
            if i == reference:
                continue

            print(f"Group {i} start!")

            mask_t = (adata.obs[groupby] == i).values
            if reference == "rest":
                mask_r = (adata.obs[groupby] != i).values
            else:
                mask_r = (adata.obs[groupby].isin([reference])).values

            data_target = df_all.loc[mask_t, :]
            data_ref    = df_all.loc[mask_r, :]

            # 1) Observed stats only (parallel) — fast and low overhead
            print("Compute dpsi (observed only)...")
            rs = Parallel(n_jobs=-1, prefer="threads", batch_size=64)(
                delayed(_compute_dpsi_observed)(e, data_ref, data_target)
                for e in event_list
            )
            result = pd.DataFrame(rs)

            # 2) p-values (fast path)
            print("Compute pvalue (fast)...")
            pvals: List[float] = []
            for x_ref, x_tgt, obs in zip(result["_x_ref"], result["_x_tgt"], result["dpsi_observed"]):
                if not np.isfinite(obs):
                    pvals.append(1.0)
                    continue
                # enforce minimal valid cells rule
                if (np.sum(~np.isnan(x_ref)) < valid_cells) or (np.sum(~np.isnan(x_tgt)) < valid_cells):
                    pvals.append(1.0)
                    continue

                if method == "approx":
                    p = _pval_approx_normal(x_ref, x_tgt, two_tailed=two_tailed)
                elif method == "perm_numba" and _HAS_NUMBA:
                    p = _perm_pval_numba(x_ref, x_tgt, n_bins, random_seed)
                elif method == "perm_py":
                    p = _perm_pval_py(x_ref, x_tgt, n_bins, random_seed)
                else:
                    # fallback: no numba -> pure python
                    p = _perm_pval_py(x_ref, x_tgt, n_bins, random_seed)

                pvals.append(float(p))

            result["pvals"] = pvals

            # drop temp arrays to save memory/serialization
            result.drop(columns=["_x_ref", "_x_tgt"], inplace=True)

            # 3) Sort: by effect size (desc) then p-value (asc)
            result = result.sort_values(by=["dpsi_observed", "pvals"], ascending=[False, True])

            # 4) Save into dicts
            evts = result["event"].astype(str).to_list()
            data_event_dict[i] = evts
            data_dpsi_observed_dict[i] = result["dpsi_observed"].astype(float).to_list()
            data_pval_dict[i] = result["pvals"].astype(float).to_list()
            data_pval_adj_dict[i] = _utils.compute_pvalue_bonferroni(data_pval_dict[i])
            data_rpsi_dict[i] = result["rpsi_value"].astype(float).to_list()
            data_tpsi_dict[i] = result["tpsi_value"].astype(float).to_list()
            data_dpr_dict[i] = result["dpr_value"].astype(float).to_list()
            data_rpr_dict[i] = result["rpr_value"].astype(float).to_list()
            data_tpr_dict[i] = result["tpr_value"].astype(float).to_list()

            # var_name (fallback to event id)
            if var_name in adata.var.columns:
                var_vals = adata[:, evts].var[var_name].astype(str).to_list()
            else:
                var_vals = evts
            var_name_dict[i] = var_vals

            print(f"Group {i} complete!")
            print("-----------------------------------------")

    # Convert to structured arrays
    name_data   = _dict_to_rec(data_event_dict)
    dpsi_data   = _dict_to_rec(data_dpsi_observed_dict)
    pval_data   = _dict_to_rec(data_pval_dict)
    pval_adj    = _dict_to_rec(data_pval_adj_dict)
    dpr_data    = _dict_to_rec(data_dpr_dict)
    rpr_data    = _dict_to_rec(data_rpr_dict)
    tpr_data    = _dict_to_rec(data_tpr_dict)
    rpsi_data   = _dict_to_rec(data_rpsi_dict)
    tpsi_data   = _dict_to_rec(data_tpsi_dict)
    vname_data  = _dict_to_rec(var_name_dict)

    # Save to adata.uns
    adata.uns[key_added] = {
        **adata.uns[key_added],
        "names":     name_data,
        "dpsi":      dpsi_data,
        "pvals":     pval_data,
        "pvals_adj": pval_adj,
        "dpr":       dpr_data,
        "rpr":       rpr_data,
        "tpr":       tpr_data,
        "rpsi":      rpsi_data,
        "tpsi":      tpsi_data,
        "gene_name": vname_data,
    }

    return adata

