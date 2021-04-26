import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from typing import List, Union, Dict, Collection, Iterable
from pathlib import Path, PosixPath, WindowsPath
import pandas as pd
from itertools import combinations
from crispr_tools.tools import drop_nonumeric, size_factor_normalise
from crispr_tools.tools import list_not_str

#todo test plot_clonal_X
#todo clonal lfcs should be multi-indexed for other uses
#  like selecting control reps that match condition x and treat reps that match y
# todo qc plots should include intial abundance vs lfc dispersion, as obtained from clonal
# todo use the same validation of



def plot_clonal_counts(count:pd.DataFrame, sample_reps:Dict[str, List[str]], file_fmt_str='',
                       title_fmt_str="Clonal counts, {}",):
    """plot the abundance of every replicate in the same sample against each other.

    Args:
        count: DataFrame of abundances. Non-numeric columns are ignored.
        sample_reps: Dictionary mapping samples to replicates.
        file_fmt_str: a string path containing {} into which will be inserted
            the sample name when saving the chart. If empty chart is not saved.
        show_plt: If true plt.show() will be called for each chart.

    Returns:
        None
    """

    if file_fmt_str:
        if '{}' not in file_fmt_str:
            raise ValueError("Save file string doesn't contain required '{}', "+str(file_fmt_str))

    count = drop_nonumeric(count)
    for samp_name, reps in sample_reps.items():
        combs = [c for c in combinations(reps, 2)]
        n = len(combs)
        fig, axes = plt.subplots(1, n, figsize=(n * 5, 5))
        if n == 1:
            axes = [axes]
        print(axes)
        for i, ((ra, rb), ax) in enumerate(zip(combs, axes)):
            if (ra not in count.columns) or (rb not in count.columns):
                print(f"Key missing, on or both of {ra}, {rb}")
                continue

            x, y = count[ra], count[rb]
            # plt.sca(ax)

            ax.hexbin(x, y, bins='log', gridsize=40)
            ax.set_xlabel(ra)
            ax.set_ylabel(rb)
            if title_fmt_str:
                ax.set_title(title_fmt_str.format(samp_name))
        plt.tight_layout()
        if file_fmt_str:
            plt.savefig(file_fmt_str.format(samp_name), dpi=150)
        else:
            plt.show()


def clustermap(df:pd.DataFrame, z_score:int=None, xlabel='',
               quantile:float=0, cmkwargs:dict=None):
    """
    not tested
    df should have samples as columns.

    args:
        z_score: 0 or 1 applies the transformation to rows (0) or columns (1),
        none uses raw data.

        quantile: use on the the most variable quantile in the heat map, 0 uses
        all data and 0.99 the top 1% most variable.

        cmkwargs: additional kwargs passed to sns.clustermap(**cmkwargs)
    """

    #todo: fix 0 = middle of the cmap
    #todo: put the color bar inthe top left
    #todo: use col_colors to identify treats and clones by default. (probably rename the function)

    if cmkwargs is None:
        cmkwargs = {}
    guide_var = df.var(1).sort_values(ascending=False)
    # most_var = guide_var.head(300).index
    most_var = guide_var > guide_var.quantile(quantile)

    # do the cluster
    cm = sns.clustermap(df.loc[most_var], cmap='coolwarm', z_score=z_score, **cmkwargs)
    # cm.savefig(f'charts/heatmap dendrogram of LFCs by clone {grp} most var.png', dpi=150)
    hm = cm.ax_heatmap.get_position()
    plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
    cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width, hm.height * 1.5])
    new_hm = cm.ax_heatmap.get_position()
    col = cm.ax_col_dendrogram.get_position()

    # set the col dendrogram y0 to the heatmap's y1
    cm.ax_col_dendrogram.set_position([col.x0, new_hm.y1 + 0.12, col.width, col.height])
    # stretch out the row dend
    rowd = cm.ax_row_dendrogram.get_position()
    cm.ax_row_dendrogram.set_position([rowd.x0, rowd.y0, rowd.width, new_hm.height])
    cm.ax_heatmap.xaxis.tick_top()
    for tick in cm.ax_heatmap.get_xticklabels():
        tick.set_rotation(90)
    cm.ax_heatmap.set_xlabel(f"{xlabel} {quantile*100}%")
    # cm.savefig(f'charts/heatmap dendrogram of LFCs by clone {grp} {quantile*100} percentile variance.png', dpi=150)



def plot_clonal_lfc(count:pd.DataFrame, sample_reps:Dict[str, List[str]],
                    ctrl_maps:Dict[str, Dict[str, List[str]]],
                    title_fmt_str="Clonal LFC, {}->{}", file_fmt_str='', show_plt=False):
    """plot the lfc of every replicate calculating lfc from paired reps
    with matching indicies in the sample_reps. LFCs calculated between
    samples specified by ctrl_maps.

    Args:
        count: DataFrame of abundances. Non-numeric columns are ignored.
            Should probably be normalised log2.
        sample_reps:  Dictionary mapping samples to replicates.
        ctrl_maps:  A dict containing control groups {ctrla:[treat1, ...], ...}
        file_fmt_str:  A string path containing 2 {} into which will be inserted
            the control and treament sample names when saving the chart.
            If empty chart is not saved.
        title_fmt_str: Used for plot title, control and treat names mapped to {}{}
        show_plt: If true plt.show() will be called for each chart.

    Returns:
        None
    """

    if file_fmt_str:
        assert file_fmt_str.count('{}') == 2
        if file_fmt_str.count('{}') != 2:
            raise ValueError("Save file string doesn't contain required '{}', "+str(file_fmt_str))
    for grp, ctrl_dict in ctrl_maps.items():
        for ctrl_samp, treat_samples in ctrl_dict.items():
            for trt_samp in treat_samples:
                # we shall assume that the clones are int he same order in sample_reps
                treat_clone_pairs = zip(sample_reps[ctrl_samp], sample_reps[trt_samp])
                # Do every pair of clones in these samples, if there are two clones there's one comb
                xy_pairs = [thing for thing in combinations(treat_clone_pairs, 2)]
                n = len(xy_pairs)
                fig, axes = plt.subplots(1, n, figsize=(n * 5, 5))
                if n == 1:
                    axes = [axes]
                for i, (x_pair, y_pair) in enumerate(xy_pairs):
                    plt.sca(axes[i])

                    # subtract log2 ctrl from log2 samp
                    lfcx = count[x_pair[1]] - count[x_pair[0]]
                    lfcy = count[y_pair[1]] - count[y_pair[0]]

                    plt.hexbin(lfcx, lfcy, bins='log', gridsize=40)
                    plt.xlabel(x_pair)
                    plt.ylabel(y_pair)
                    if title_fmt_str:
                        plt.title(title_fmt_str.format(ctrl_samp, trt_samp))

                    plt.plot([-1, 1], [-1, 1], 'w--')

                plt.tight_layout()
                plt.savefig(file_fmt_str.format(ctrl_samp, trt_samp), dpi=150)
                if show_plt:
                    plt.show()

