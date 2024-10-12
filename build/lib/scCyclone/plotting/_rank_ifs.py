# -*- coding: utf-8 -*-
"""
@File    :   _rank_ifs.py
@Time    :   2024/10/10
@Author  :   Dawn
@Version :   1.0
@Desc    :   bug
"""


from collections.abc import Iterable, Sequence, Mapping
from typing import Literal
import pandas as pd
import numpy as np

from anndata import AnnData

from ..get import get

from scanpy._utils import savefig_or_show
from scanpy._settings import settings


def _fig_show_save_or_axes(plot_obj, return_fig, show, save):
    """
    Decides what to return
    """
    if return_fig:
        return plot_obj
    plot_obj.make_figure()
    savefig_or_show(plot_obj.DEFAULT_SAVE_PREFIX, show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return plot_obj.get_axes()


def _get_values_to_plot(
    adata,
    values_to_plot: Literal[
        "dif",
        "pvals",
        "pvals_adj",
        "dpr",
    ],
    gene_names: Sequence[str],
    *,
    groups: Sequence[str] | None = None,
    key: str | None = "rank_ifs_groups",
    gene_symbols: str | None = None,
):
    """
    If rank_genes_groups has been called, this function
    prepares a dataframe containing scores, pvalues, logfoldchange etc to be plotted
    as dotplot or matrixplot.

    The dataframe index are the given groups and the columns are the gene_names

    used by rank_genes_groups_dotplot

    Parameters
    ----------
    adata
    values_to_plot
        name of the value to plot
    gene_names
        gene names
    groups
        groupby categories
    key
        adata.uns key where the rank_genes_groups is stored.
        By default 'rank_genes_groups'
    gene_symbols
        Key for field in .var that stores gene symbols.
    Returns
    -------
    pandas DataFrame index=groups, columns=gene_names

    """
    valid_options = [
        "dif",
        "pvals",
        "pvals_adj",
        "dpr",
    ]
    if values_to_plot not in valid_options:
        raise ValueError(
            f"given value_to_plot: '{values_to_plot}' is not valid. Valid options are {valid_options}"
        )

    values_df = None
    check_done = False
    if groups is None:
        groups = adata.uns[key]["names"].dtype.names
    if values_to_plot is not None:
        df_list = []
        for group in groups:
            df = get.rank_ifs_groups_df(adata, group, key=key)
            if gene_symbols is not None:
                df["names"] = df[gene_symbols]
            # check that all genes are present in the df as sc.tl.rank_genes_groups
            # can be called with only top genes
            if not check_done:
                if df.shape[0] < adata.shape[1]:
                    message = (
                        "Please run `sc.tl.rank_genes_groups` with "
                        "'n_genes=adata.shape[1]' to save all gene "
                        f"scores. Currently, only {df.shape[0]} "
                        "are found"
                    )
                    raise ValueError(message)
            df["group"] = group
            df_list.append(df)

        values_df = pd.concat(df_list)
       
        column = values_to_plot
        values_df = pd.pivot(
            values_df, index="names", columns="group", values=column
        ).fillna(1)

        values_df = values_df.loc[gene_names].T

    return values_df




def _rank_ifs_groups_plot(
    adata: AnnData,
    plot_type: str = "heatmap",
    *,
    groups: str | Sequence[str] | None = None,
    n_genes: int | None = None,
    groupby: str | None = None,
    values_to_plot: str | None = None,
    var_names: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
    min_logfoldchange: float | None = None,
    key: str | None = None,
    show: bool | None = None,
    save: bool | None = None,
    return_fig: bool | None = False,
    gene_symbols: str | None = None,
    **kwds,
):
    """\
    Common function to call the different rank_genes_groups_* plots
    """
    if var_names is not None and n_genes is not None:
        raise ValueError(
            "The arguments n_genes and var_names are mutually exclusive. Please "
            "select only one."
        )

    if var_names is None and n_genes is None:
        # set n_genes = 10 as default when none of the options is given
        n_genes = 10

    if key is None:
        key = "rank_genes_groups"

    if groupby is None:
        groupby = str(adata.uns[key]["params"]["groupby"])
    group_names = adata.uns[key]["names"].dtype.names if groups is None else groups

    if var_names is not None:
        if isinstance(var_names, Mapping):
            # get a single list of all gene names in the dictionary
            var_names_list = sum([list(x) for x in var_names.values()], [])
        elif isinstance(var_names, str):
            var_names_list = [var_names]
        else:
            var_names_list = var_names
    else:
        # dict in which each group is the key and the n_genes are the values
        var_names = {}
        var_names_list = []
        for group in group_names:
            df = get.rank_ifs_groups_df(
                adata,
                group,
                key=key,
                gene_symbols=gene_symbols,
                log2fc_min=min_logfoldchange,
            )

            if gene_symbols is not None:
                df["names"] = df[gene_symbols]

            genes_list = df.names[df.names.notnull()].tolist()

            if len(genes_list) == 0:
                raise ValueError(f"No genes found for group {group}")
                continue
            if n_genes < 0:
                genes_list = genes_list[n_genes:]
            else:
                genes_list = genes_list[:n_genes]
            var_names[group] = genes_list
            var_names_list.extend(genes_list)

    # by default add dendrogram to plots
    kwds.setdefault("dendrogram", True)

    if plot_type in ["dotplot", "matrixplot"]:
        # these two types of plots can also
        # show score, logfoldchange and pvalues, in general any value from rank
        # genes groups
        title = None
        values_df = None
        if values_to_plot is not None:
            values_df = _get_values_to_plot(
                adata,
                values_to_plot,
                var_names_list,
                key=key,
                gene_symbols=gene_symbols,
            )
            title = values_to_plot
            if values_to_plot == "logfoldchanges":
                title = "log fold change"
            else:
                title = values_to_plot.replace("_", " ").replace("pvals", "p-value")

        if plot_type == "dotplot":
            from scanpy.plotting._dotplot import dotplot

            _pl = dotplot(
                adata,
                var_names,
                groupby,
                dot_color_df=values_df,
                return_fig=True,
                gene_symbols=gene_symbols,
                **kwds,
            )
            if title is not None and "colorbar_title" not in kwds:
                _pl.legend(colorbar_title=title)
        elif plot_type == "matrixplot":
            from scanpy.plotting._matrixplot import matrixplot

            _pl = matrixplot(
                adata,
                var_names,
                groupby,
                values_df=values_df,
                return_fig=True,
                gene_symbols=gene_symbols,
                **kwds,
            )

            if title is not None and "colorbar_title" not in kwds:
                _pl.legend(title=title)

        return _fig_show_save_or_axes(_pl, return_fig, show, save)

    elif plot_type == "stacked_violin":
        from scanpy.plotting._stacked_violin import stacked_violin

        _pl = stacked_violin(
            adata,
            var_names,
            groupby,
            return_fig=True,
            gene_symbols=gene_symbols,
            **kwds,
        )
        return _fig_show_save_or_axes(_pl, return_fig, show, save)
    elif plot_type == "heatmap":
        from scanpy.plotting._anndata import heatmap

        return heatmap(
            adata,
            var_names,
            groupby,
            show=show,
            save=save,
            gene_symbols=gene_symbols,
            **kwds,
        )

    elif plot_type == "tracksplot":
        
        from scanpy.plotting._anndata import tracksplot

        return tracksplot(
            adata,
            var_names,
            groupby,
            show=show,
            save=save,
            gene_symbols=gene_symbols,
            **kwds,
        )
