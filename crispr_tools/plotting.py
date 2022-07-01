import matplotlib.pyplot as plt
from crispr_tools.tools import size_factor_normalise, ROC_values, clonal_lfcs
import numpy as np
import pandas as pd
import seaborn as sns

import itertools
from typing import Tuple, List, Dict, Iterable, Collection

from scipy import stats

try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    def adjust_text(*args, **kwargs):
        pass

def plot_read_violins(tab, samp_per_row='all', column_labs=None, log=True, size_norm=False, ax=None):
    """Takes counts table, does log2 (optionally size factor normalises)
    and produces violin density plots.

    figsize is (3*n, 6)

    Returns a plt.Axes

    per_row specifies the number of violins on each row"""

    n_samp = tab.shape[1]

    if samp_per_row == 'all':
        samp_per_row = n_samp

    rows = n_samp // samp_per_row + int((n_samp % samp_per_row) > 0)

    if ax is None:
        fig, axes = plt.subplots(rows, 1, figsize=(3 * samp_per_row, 6 * rows))
    else:
        axes = ax

    if rows == 1:
        axes = [axes]

    # remove nonnumeric columns
    tab = tab.copy()
    # drop any text columns
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )
    # transform the data
    if size_norm:
        tab = size_factor_normalise(tab, log=log)
    elif log:
        tab = tab.__add__(1).apply(np.log2)

    # apply specified column labels
    if column_labs is None:
        column_labs = tab.columns

    for rowi, ax in enumerate(axes):
        inds = slice(rowi * samp_per_row, (1 + rowi) * samp_per_row)
        stab = tab.iloc[:, inds]

        # this is poo
        if stab.shape[0] == 0:
            break
        ax.violinplot(stab, widths=0.8)
        ax.boxplot(stab, widths=0.2)

        ax.set_xticks(range(1, samp_per_row + 1))
        clabs = list(column_labs[inds])
        while len(clabs) < samp_per_row:
            clabs.append('')
        ax.set_xticklabels(clabs, rotation=40)
    plt.tight_layout()
    if rows == 1:
        return axes[0]
    else:
        return axes


def plot_ROC(tab, things_oi, label=None, ax = None):
    #todo: make things_oi work as any kind of iterable (currently pd.Series)
    """Tab needs to be a series or dataframe containing values used to
    order the tab index."""

    if type(things_oi) is list:
        things_oi = pd.Series(things_oi, index=things_oi)

    if ax is None:
        _, ax = plt.subplots(1,1, figsize=(4,4))
    else:
        plt.sca(ax)
    tab = tab.copy()
    things_oi = things_oi[things_oi.isin(tab.index)]
    #print(things_oi)
    if len(tab.shape) > 1:
        for i in range(tab.shape[1]):
            if label is not None:
                lab = label[i]
            else:
                lab = tab.columns[i]
            # get ordered list of the values and subset things of interest to exclude nan values.
            col = tab.iloc[:, i].dropna().sort_values()
            toi = things_oi[things_oi.isin(col.index)]
            plt.plot(np.arange(col.shape[0])/col.shape[0],
                     np.cumsum(col.index.isin(toi))/len(toi),
                     label=lab)
    else:
        tab = tab.sort_values()
        plt.plot(np.arange(tab.shape[0])/tab.shape[0], np.cumsum(tab.index.isin(things_oi))/len(things_oi))
    #     plt.plot(np.arange(kin.shape[0])/kin.shape[0], np.cumsum(kin.index.isin(gn))/len(gn),
    #             label='RPE1 screen essentials')
    plt.plot([0,1], [0,1], 'k--', alpha=0.4)
    plt.xlabel('Not essential')
    plt.ylabel('Essential')
    plt.title('ROC')
    plt.legend()
    plt.tight_layout()
    return ax


def plot_volcano_from_mageck(tab, title='', label_genes=None, outfn='', ax=None, source='mageck'):
    """Take a mageck table and do a volcano plot. """

    dep = tab.loc[:, 'pos|lfc'] > 0
    enr = tab.loc[:, 'pos|lfc'] < 0
    tab.loc[dep, 'ml10_fdr'] = -np.log10(tab.loc[dep, 'pos|fdr'])
    tab.loc[enr, 'ml10_fdr'] = -np.log10(tab.loc[enr, 'neg|fdr'])
    plot_volcano(tab['pos|lfc'], tab['ml10_fdr'], title=title, other_labels=label_genes,
                 outfn=outfn, ax=ax)

#POINT = 0

def plot_volcano(lfc, fdr, tab=None, title='', label_deplet=0, label_enrich=0,
                 other_labels=None, p_thresh=0.05, outfn='', ax=None,
                 exclude_labs=('NonT', 'Rando'), plot_kw: dict = None):
    """Draw a volcano plot of lfc vs fdr. assumes fdr is -log10.

    :param lfc: str giving tab[lfc] or series with gene names as index
    :param fdr: str giving tab[fdr] or series with gene names as index
    :param tab: DataFrame containing lfc and fdr
    :param title: optional plot title
    :param label_enrich: int giving top n enriched genes to label
    :param label_deplet: int giving top n depleted genes to label
    :param other_labels: other genes to label
    :param p_thresh:
    :param outfn: if provided a .png will be written
    :param ax: optional plt.Axes instance to use
    :return: plt.Axes
    """
    # print('TP ', POINT)

    # this is silly
    lfc_lab, fdr_lab = None, None
    if tab is not None:
        lfc_lab = lfc
        fdr_lab = fdr
        lfc = tab[lfc]
        fdr = tab[fdr]

    sctkw = dict(marker='o', linestyle='none', alpha=0.4)
    if plot_kw is not None:
        sctkw.update(plot_kw)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 12))

    # todo update to use scatter and set minimum matplotlib version
    ax.plot(lfc, fdr, **sctkw)
    # plt.yscale('log')
    # plt.gca().invert_yaxis()
    ax.set_title(title)

    xmin, xmax, ymin, ymax = ax.axis()

    p_thresh = -np.log10(p_thresh)
    ax.plot([xmin, xmax], [p_thresh, p_thresh], 'k--')
    ax.plot([0, 0], [ymin, ymax], color='silver')
    # label the top and bottom most, if specified
    texts = []
    texts_done = []

    # get subtables
    # filter out excluded labels by getting Series containing only included labs,
    # turn that into a mask of the original series and combining that with the dep/enr masks.
    filtered_lfc = lfc.copy()
    if exclude_labs:
        for exclude in exclude_labs:
            filtered_lfc = filtered_lfc.loc[~filtered_lfc.index.str.contains(exclude)]

        included = lfc.index.isin(filtered_lfc.index)
        depmask = (lfc < 0) & included
        enrmask = (lfc > 0) & included
    else:
        depmask = (lfc < 0)
        enrmask = (lfc > 0)

    # get tails from ascending order
    if tab is not None:
        dep = tab.loc[depmask, :].sort_values(lfc_lab, ascending=False).sort_values(fdr_lab).tail(label_deplet)[fdr_lab]
        enr = tab.loc[enrmask, :].sort_values([fdr_lab, lfc_lab]).tail(label_enrich)[fdr_lab]
        # label the tails
        for end in dep, enr:
            for lab, an_fdr in end.items():
                if an_fdr < p_thresh:
                    continue
                texts_done.append(lab)
                texts.append(plt.text(lfc[lab], fdr[lab], lab))
    # label additional genes
    if other_labels is not None:
        for lab in other_labels:
            if lab in texts_done:
                continue
            texts.append(
                ax.text(lfc[lab], fdr[lab], lab)
            )

    if texts and adjust_text:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), ax=ax)

    if lfc_lab is None:
        ax.set_xlabel('Log$_2$ Fold Change')
    else:
        ax.set_xlabel(lfc_lab)

    ax.set_ylabel('-log$_{10}$(FDR)')

    if outfn:
        plt.tight_layout()
        plt.savefig(outfn, dpi=150)
    return ax


def scores_scatter_plotly(x, y, fdr,
                          fdr_thresholds=(0.05,),
                          fdr_colours=('#CC7700',),
                          gene_set_masks: List[Tuple] = None, ):
    """gene_set_masks, [('label1', mask1), ('label2', mask2), ...]. Mask here either a
    pd.Series boolean mask, or list of genes that can be passed to .loc[]"""
    from plotly import graph_objects as go
    fig = go.Figure()
    fig.update_layout(template='plotly_white')
    Xall = x
    Yall = y

    # These will be populated with values and passed to go.Scatter
    xys = []
    marker_labels = []
    markers = []
    symbols = []

    # get list of genes between FDR thresholds
    #  Sort thresholds and colours, largest to smallest
    sranks = stats.rankdata(fdr_thresholds, 'ordinal')
    fdr_thresholds = sorted(fdr_thresholds, key=lambda x: sranks[fdr_thresholds.index(x)], reverse=True)
    fdr_colours = sorted(fdr_colours, key=lambda x: sranks[fdr_colours.index(x)], reverse=True)

    # rank big to small
    sranks = stats.rankdata(fdr_thresholds, 'ordinal')
    fdr_thresholds = [1] + sorted(fdr_thresholds, key=lambda x: sranks[fdr_thresholds.index(x)], reverse=True)
    fdr_colours = ['#d9d9d9'] + sorted(fdr_colours, key=lambda x: sranks[fdr_colours.index(x)], reverse=True)

    # get masks giving a gene's fdr bracket
    sig_masks = []
    genes_in_mask = []
    for i in range(len(fdr_thresholds) - 1):
        m = (fdr_thresholds[i + 1] < fdr) & (fdr <= fdr_thresholds[i])
        m = m[m].index
        sig_masks.append(m)
        genes_in_mask.extend(m)
    m = fdr < fdr_thresholds[-1]
    sig_masks.append(m[m].index)

    # populate the fdr related Scatter values
    for m, c in zip(sig_masks, fdr_colours):
        xys.append((Xall.loc[m], Yall.loc[m]))
        symbols.append('circle-open')
        markers.append(dict(color=c, size=7))

    marker_labels.append(f"FDR > {fdr_thresholds[1]}")
    for thresh in fdr_thresholds[1:]:
        marker_labels.append(f"FDR < {thresh}")
    # finished with FDR

    # populate gene set scatter values
    if gene_set_masks is not None:
        for lab, m in gene_set_masks:
            xys.append((Xall.loc[m], Yall.loc[m]))
            symbols.append('circle')
            markers.append(dict(size=4))
            marker_labels.append(lab)

    # plot it finally
    for (x, y), data_label, marker, symbol in zip(xys, marker_labels, markers, symbols):
        fig.add_trace(
            go.Scatter(
                x=x, y=y,
                name=data_label,
                mode='markers',
                marker_symbol=symbol,
                marker=marker,
                text=x.index,
                customdata=fdr.loc[x.index],
                hovertemplate='%{text}:<br>  %{x:.2f},%{y:.2f}<br>  FDR=%{customdata:.2f}',
            )
        )
    return fig

def write_plotly_html(results_table:pd.DataFrame,
                      xy_key:str,
                      samplePairs:Iterable=None,
                      gene_set_masks:Tuple[str, pd.Series]=None,
                      fdr_thresholds=(0.05,),
                      fdr_colours=('#CC7700',),
                      out_fmt_str='',
                      samp_labels=Dict[str, str],
                      samp_pair_key_str='{}-{}',
                      show_fig=False):
    """
    Produces set of HTMLs applying scores_scatter_plotly using a table, but I
    don't really get the format of the input table and need to go back and look
    at where I originally used this.
    TODO: how do you use write_plotly_html

    Args:
        results_table: DF with multiindex columns, (sample, stat). Must include
            neg_fdr & pos_fdr stat columns
        xy_key: The stat name (in results_table) that will be used as x/y values
        samplePairs: All X/Y samples that will have html produced, all pairs by
            default.
        gene_set_masks: List of labelling strings and bool masks that indicate
            subsets of genes that will be plotted separately as small dots within
            the main points.
        out_fmt_str: Path to which files will be written, containing one {} to
            which the comparison will be written.
        samp_labels: Dictionary of sample labels
        samp_pair_key_str:
    """

    if samplePairs is None:
        if type(results_table) is pd.DataFrame:
            samps = results_table.columns.levels[0]
        else:
            samps = results_table.keys()
        samplePairs = itertools.permutations(samps, 2)

    for rsamp, psamp in samplePairs:

        smpk = samp_pair_key_str.format(rsamp, psamp)

        Xall, Yall = results_table[smpk][xy_key]

        enr, dep = [results_table[smpk][pk] for pk in ('pos_fdr', 'neg_fdr')]

        # get the lowest fdr for hovertext
        fdr_lowest = enr.copy()
        # wherever dep is smaller, overwrite the enr fdr value
        m = dep < enr
        fdr_lowest[m] = dep[m]

        axis_labels = [samp_labels[s] if samp_labels else s for s in (rsamp, psamp)]

        fig = scores_scatter_plotly(Xall, Yall, fdr_lowest, fdr_thresholds, fdr_colours, gene_set_masks)

        # comp = f"{rsamp}-{psamp}"
        if out_fmt_str:
            fn = out_fmt_str.format(rsamp, psamp)
            fig.write_html(
                fn,
                include_plotlyjs='directory',
            )
        if show_fig:
            fig.show()


def plot_rank_vs_score(score:pd.Series, sig:pd.Series, n_labels, sig_threshold=0.05, step=0.015):
    """Score is the thing they will be ranked on,
    sig series is just used for drawing the significance line
    todo: make sig optional."""

    def get_rank_of_thresholds(score, sig, threshold=0.05):

        # first dropouts, then enrichments
        threshold_yvalue = []
        for endmask in (score < 0), (score > 0):
            yend = sig.loc[endmask]
            thresh_rank = (yend < threshold).sum()
            threshold_yvalue.append(
                score.loc[endmask].abs().sort_values(ascending=False).iloc[thresh_rank]
            )

        threshold_yvalue[0] = -threshold_yvalue[0]
        return threshold_yvalue

    plt.figure(figsize=(7, 10))
    #tab = res_drugz[k]
    normz_rank = score.rank()
    plt.scatter(normz_rank, score)
    low, hi = get_rank_of_thresholds(score, sig, sig_threshold)
    for thresh in (low, hi):
        plt.plot([0, score.shape[0]], [thresh, thresh], 'k--', alpha=0.6)
        plt.text(score.shape[0] / 2, thresh + 0.1, f'{sig_threshold * 100}% FDR')

    for isneg, endmask in (True, (score < 0)), (False, (score > 0)):
        ax = plt.gca()
        axis_to_data = (ax.transAxes + ax.transData.inverted()).transform
        scoreend = score.loc[endmask]

        top_genes = scoreend.loc[(sig.loc[endmask] < sig_threshold)].abs().sort_values(ascending=False).head(n_labels).index
        try:
            extreme_value = scoreend.loc[top_genes[0]]
        except IndexError:
            continue

        for gi, gn in enumerate(top_genes):
            x = normz_rank[gn]
            y = score[gn]


            if isneg:
                text_x, text_y = axis_to_data([0.25, 0.05 + gi * step])
                alignment = 'left'

            else:
                text_x, text_y = axis_to_data([0.75, 1 - (gi * step)])
                alignment = 'right'

            plt.text(
                text_x, text_y, gn,
                fontsize=9,
                bbox={'facecolor': 'white', 'alpha': 1, 'pad': 0.1},
                ha=alignment,
            )
            plt.plot([x + 50, text_x, ], [y, text_y], 'k-')

def pca_grid(pca, hue_deet, style_deet, max_components=5, also_return_fig=False):
    """Plot components against each other in a grid of scatter plots.
    Deets MUST be in the same order as the columns used for the PCA.
    Args:
        pca: a sklearn.decomposition.PCA() object that has been fit
        hue_deet: grouping variables used for setting colours, by sample
        style_deet: as hue_deet but marker style
        max_components: components to plot, max_comp*(max_comp-1) scatter plots will be produced
        also_return_fig: When false only axes are returned, set to true to also return the Figures"""


    thing = sns.scatterplot(x=pca.components_[0], y=pca.components_[0],
                            hue=hue_deet, style=style_deet,
                            s=150)

    leg = thing.get_legend_handles_labels()
    plt.close()

    max_b = max_components - 1
    fig, axes = plt.subplots(max_components, max_b, figsize=(4.5 * max_b, 3.85 * max_components))
    # pc ind also used for subplots
    for pc_a in range(max_components):
        for pc_b in range(max_b):
            plt.sca(axes[pc_a][pc_b])

            if pc_a == pc_b:
                plt.legend(*leg)
                continue

            # pc_a/b swapped to get column PC on the x and row PC on the y
            sns.scatterplot(x=pca.components_[pc_b], y=pca.components_[pc_a],
                            hue=hue_deet, style=style_deet,
                            s=150, legend=False)
            plt.xlabel(f"PC {pc_b+1} (variance explained: {pca.explained_variance_ratio_[pc_b]*100:.3}%)")
            plt.ylabel(f"PC {pc_a+1} (variance explained: {pca.explained_variance_ratio_[pc_a]*100:.3}%)")
    plt.tight_layout()
    if also_return_fig:
        return (fig, axes)
    else:
        return axes

def get_req_infection(libcnt, minimum=100, proportion = 0.95):
    """Find the minimum number of infected cells required to get a specified
    proportion of guides above a minimum abundance.

    Defaults to 95% of guides with abundance > 100."""
    libcnt = libcnt.sort_values().copy()
    x = libcnt.quantile(1-proportion)
    if x == 0:
        return np.inf
    return sum(libcnt)/(x/minimum)

def plot_req_inf(counts, reps, qrange=(0.99,0.90), moi=1):

    plt.figure(figsize=(6, 10))
    for n in reps:
        ys=[]
        xs = np.linspace(*qrange, 100)
        for quantile in xs:
            ys.append(get_req_infection(counts, n, proportion = quantile))
        ys=pd.Series(ys)
        plt.plot((1-xs)*100, ys/moi, label=str(n))
        plt.xlabel("% guides > X ")
        plt.ylabel("required cell infections")
    plt.legend(title='X (guide rep)')