import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from typing import List, Union, Dict, Collection, Iterable
from pathlib import Path, PosixPath, WindowsPath
import pandas as pd
from itertools import combinations
from crispr_tools.tools import drop_nonumeric

#todo test plot_clonal_X

def plot_clonal_counts(count:pd.DataFrame, sample_reps:Dict[str, List[str]], file_fmt_str='',
                       title_fmt_str="Clonal counts, {}", show_plt=False):
    """plot the abundance of every replicate in the same sample against each other.

    Args:
        count: DataFrame of abundances. Non-numeric columns are ignored.
            Should probably be log2.
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
        for i, (ra, rb) in enumerate(combs):
            plt.sca(axes[i])
            #plt.figure(figsize=(5, 5))
            plt.hexbin(count[ra], count[rb], bins='log', gridsize=40)
            plt.xlabel(ra)
            plt.ylabel(rb)
            if title_fmt_str:
                plt.title(title_fmt_str.format(samp_name))
            plt.tight_layout()
            if file_fmt_str:
                plt.savefig(file_fmt_str.format(samp_name), dpi=150)
            if show_plt:
                plt.show()
            else:
                plt.close()


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
                    plt.xlabel(x_pair);
                    plt.ylabel(y_pair)
                    if title_fmt_str:
                        plt.title(title_fmt_str.format(ctrl_samp, trt_samp))

                    plt.plot([-1, 1], [-1, 1], 'w--')

                plt.tight_layout()
                plt.savefig(file_fmt_str.format(ctrl_samp, trt_samp), dpi=150)
                plt.show()