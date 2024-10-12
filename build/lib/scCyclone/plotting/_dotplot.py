# -*- coding: utf-8 -*-
"""
@File    :   _dotplot.py
@Time    :   2024/10/10
@Author  :   Dawn, ziqi
@Version :   1.0
@Desc    :   Plot dotplot switch consequences 
"""

import dotplot

from ..get import get


def dotplot_switch_consequences(adata, 
                 key='rank_switch_consequences_groups', 
                 colcmap='coolwarm',
                 vmin: float = -0.5, 
                 vmax: float = 0.5,
                 dot_title='Fraction of cells\n in group',
                 colorbar_title='log fold change',
                 cluster_order=None,
                 **kwargs):
    """
    Plots a dot plot for the rank switch consequences in the AnnData object.

    Parameters:
    ----------
    adata: AnnData
        The AnnData object containing the data.
    key: str (default: 'rank_switchs_consequences_groups')
        The key in adata.uns to access the rank switch consequences data.
    colcmap: str (default: 'coolwarm')
        The colormap to use for the dot plot.
    vmin: float (default: -0.5)
        The minimum value for the color range.
    vmax: float (default: 0.5)
        The maximum value for the color range.
    dot_title: str (default: 'Fraction of cells\n in group')
        The title for the dot plot.
    colorbar_title: str (default: 'log fold change')
        The title for the color bar.
    cluster_order: list (default: None)
        The order of the clusters to plot.
    **kwargs:
        Additional keyword arguments to pass to the dot plot function.
    """

    # 获取转换后果的DataFrame
    df = get.rank_switch_consequences_groups_df(adata, key=key)
    
    # 根据'group'和'feature'聚合数据，计算每个组中细胞的百分比
    df_sizes = df.pivot(index='group', columns='feature', values="percentage")
    
    # 根据'group'和'feature'聚合数据，计算每个组的log2 fold change
    df_colors = df.pivot(index='group', columns='feature', values="log2fc")
    
    # 如果提供了cluster_order，则按照指定的顺序对数据进行重新排序
    if cluster_order:
        df_sizes = df_sizes.loc[cluster_order]
        df_colors = df_colors.loc[cluster_order]
    
    # 创建DotPlot对象
    dp = dotplot.DotPlot(df_size=df_sizes, df_color=df_colors)
    
    # 绘制点图，设置大小因子、颜色映射、标题、颜色范围和颜色条标题
    _ = dp.plot(size_factor=100, cmap=colcmap, dot_title=dot_title,
                vmin=vmin, vmax=vmax, colorbar_title=colorbar_title, **kwargs)
    
    


