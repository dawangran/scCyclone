# -*- coding: utf-8 -*-
"""
@File    :   _isoform.py
@Time    :   2024/10/10
@Author  :   Dawn, ziqi
@Version :   1.0
@Desc    :   Plot transcript structure 
"""


import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from tkinter import *
import anndata as ad



def transcript_structure(
    adata: ad.AnnData, 
    gene_name: str,
    var_name:str,
    output_path=str
    ):
    
    """
    Plots a transcript diagram for the given gene name.

    Parameters:
    ----------
    adata: AnnData
        The AnnData object containing the transcript information.
    gene_name: str
        The name of the gene to plot.
    var_name: str
        The variable name in adata.var to match the gene name.
    output_path: str
        The path to save the plot.
    """
    
    transcript_list=list(adata.var[adata.var[var_name]==gene_name].index)
    strand=adata.var[adata.var['gene_name']==gene_name]['strand'][0]
    if strand=="+":
        arr = '-|>'
    else:
        arr = '<|-'
        
    transcript_num = len(transcript_list)
    line_width = 10
    num = 0
    tmp_colors = ['#2166ac', 'indigo']
    names_tmp_colors = ['gene', 'exon']
    colors_legend_name = ['gene', 'exon']
    color_dict = dict(zip(names_tmp_colors, tmp_colors))

    exon_dict=adata.uns['isoform_exon']

    filter_exon_dict = {key: exon_dict[key] for key in transcript_list if key in exon_dict}


    filter_exon_values = list(filter_exon_dict.values())
    flattened_list = [item for sublist1 in filter_exon_values for sublist2 in sublist1 for item in sublist2]

    gene_start=min(flattened_list)
    gene_end=max(flattened_list)

    fig = plt.figure(1)
    ax = fig.add_axes([0.2, 0.2, 0.5, 0.6])
    arrow = mpatches.FancyArrowPatch(
    (int(gene_start), 0.1),
    (int(gene_end), 0.1),
    arrowstyle=arr,
    mutation_scale=25, lw=2, color='#2166ac', antialiased=True)
    ax.add_patch(arrow)
    ax.set_xlim(int(gene_start), int(gene_end))
    ax.set_ylim(-0.5, transcript_num + 1)
    ax.set_xticks(np.linspace(int(gene_start), int(gene_end), 5))
    ax.set_yticks([0.1] + list(range(1, transcript_num + 1)))
    ax.set_yticklabels([gene_name]+transcript_list)
    ax.set_xticklabels([str(j) for j in [int(i) for i in np.linspace
        (int(gene_start), int(gene_end),5)]])

    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.get_xaxis().set_tick_params(direction='out')
    ax.tick_params(axis=u'y', which=u'both', length=0)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)

    for i in transcript_list:
        num += 1
        isoform_info=filter_exon_dict[i]
        t_start=min(isoform_info[0])
        t_end=max(isoform_info[1])
        line1 = [(int(t_start), num), (int(t_end), num)]
        (line1_xs, line1_ys) = zip(*line1)
        arrow = mpatches.FancyArrowPatch(
        (line1_xs[0], line1_ys[0]),
        (line1_xs[1], line1_ys[1]),
        arrowstyle=arr,
        mutation_scale=10, lw=2, color='black', antialiased=True)
        ax.add_patch(arrow)
        
        for s,e in zip(isoform_info[0],isoform_info[1]):
            line2 = [(int(s) - 0.5, num), (int(e) + 0.5, num)]
            (line2_xs, line2_ys) = zip(*line2)
            ax.add_line(lines.Line2D(line2_xs, line2_ys,
                                    solid_capstyle='butt', solid_joinstyle='miter',
                                    linewidth=int(line_width), alpha=1,
                                    color=color_dict["exon"],
                                    antialiased=False))

    ax_legend = fig.add_axes([0.75, 0.2, 0.2, 0.4])
    for i in range(len(colors_legend_name)):
        line3 = [(0, (9 - i) * 0.1), (0.1, (9 - i) * 0.1)]
        (line3_xs, line3_ys) = zip(*line3)
        ax_legend.add_line(lines.Line2D(line3_xs, line3_ys, linewidth=5,
                                        color=color_dict[names_tmp_colors[i]],
                                        solid_capstyle='butt', solid_joinstyle='miter',
                                        antialiased=False))
        ax_legend.text(0.2, (8.9 - i) * 0.1, colors_legend_name[i], fontsize=6)
        ax_legend.set_axis_off()
            
    if output_path:
        fig.savefig(output_path, dpi=800)
        
