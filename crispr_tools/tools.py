import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

from pathlib import Path, PosixPath, WindowsPath
import pandas as pd

__version__ = 'v1.3.0'

#v1.3.0
# Adding Mageck tab
# fixed plot volcanos enrichment labeling

#todo: sometimes lfc tables not writte, make it optional
#todo: boxpot on the violins
#todo: pass through non numeric columns in SFN
#todo; pass plot_volcano a filen string and it loads the table


def drop_nonumeric(tab):
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )
    return tab

def size_factor_normalise(cnt_tab, log=True):
    """The number by which MAGeCK uses to normalise its reads.
    Returns factor for each column.

    Uses +1 counts, returns a copy of the table.
    log=True returns np.log2 of the counts."""

    cnt_tab = cnt_tab.copy()

    cnt_tab = drop_nonumeric(cnt_tab)

    cnt_tab = cnt_tab + 1

    # cnt_tab = cnt_tab.apply(np.log2)
    # get geometric means, by guide
    gm = stats.gmean(cnt_tab, axis=1)
    # divide the counts by the gmean for each guide
    tab = cnt_tab.T / gm
    # get the median gmean for each experiment and normalise by that
    norm_tab = cnt_tab / tab.T.apply(np.median)
    #norm_tab = norm_tab.apply(round)
    if log:
        norm_tab = norm_tab.apply(np.log2)
    return norm_tab


def plot_read_violins(tab, column_labs=None, log=True, size_norm=False, ax=None):
    """Takes counts table, does log2 (optionally size factor normalises)
    and produces violin density plots.

    Returns a plt.Axes"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(3 * len(tab.columns), 6))

    tab = tab.copy()
    # drop any text columns
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )
    if size_norm:
        tab = size_factor_normalise(tab, log=log)
    elif log:
        tab = tab.__add__(1).apply(np.log2)

    if column_labs is None:
        column_labs = tab.columns

    ax.violinplot(tab.T)
    ax.boxplot(tab.T, widths=0.2)
    ax.set_xticks(range(1, len(column_labs) + 1))
    ax.set_xticklabels(column_labs, rotation=40)
    return ax


def plot_ROC(tab, things_oi, label=None):
    """Tab needs to be a series or dataframe containing values used to
    order the tab index."""
    plt.figure(figsize=(4,4))
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
        plt.plot(np.arange(tab.shape[0])/tab.shape[0], np.cumsum(tab.index.isin(things_oi))/len(things_oi))
#     plt.plot(np.arange(kin.shape[0])/kin.shape[0], np.cumsum(kin.index.isin(gn))/len(gn),
#             label='RPE1 screen essentials')
    plt.plot([0,1], [0,1], 'k--', alpha=0.4)
    plt.xlabel('Not essential')
    plt.ylabel('Essential')
    plt.title('ROC')
    plt.legend()
    plt.tight_layout()


def plot_volcano_from_mageck(tab, title='', label_genes=None, outfn='', ax=None, source='mageck'):
    """Take a mageck table and do a volcano plot. """

    dep = tab.loc[:, 'pos|lfc'] > 0
    enr = tab.loc[:, 'pos|lfc'] < 0
    tab.loc[dep, 'ml10_fdr'] = -np.log10(tab.loc[dep, 'pos|fdr'])
    tab.loc[enr, 'ml10_fdr'] = -np.log10(tab.loc[enr, 'neg|fdr'])
    plot_volcano(tab['pos|lfc'], tab['ml10_fdr'], title=title, other_labels=label_genes,
                 outfn=outfn, ax=ax)


def plot_volcano(lfc, fdr, tab=None, title='', label_deplet=0, label_enrich=0,
                 other_labels=None, p_thresh=0.05, outfn='', ax=None, exclude_labs = ('NonT',)):
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

    if tab is not None:
        lfc = tab[lfc]
        fdr = tab[fdr]

    sctkw = dict(marker='o', linestyle='none')
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 12))
    ax.plot(lfc, fdr, alpha=0.4, **sctkw)
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

    filtered_lfc = lfc.copy()
    for exclude in exclude_labs:
        filtered_lfc = filtered_lfc.loc[~filtered_lfc.index.str.contains(exclude)]

    depmask = lfc < 0
    enrmask = lfc > 0
    # ascending order by default
    dep = fdr[depmask].sort_values().tail(label_deplet)
    enr = fdr[enrmask].sort_values().tail(label_enrich)

    for end in dep, enr:
        for lab, an_fdr in end.items():
            if an_fdr < p_thresh:
                continue
            texts_done.append(lab)
            texts.append(plt.text(lfc[lab], fdr[lab], lab))
    # label additional genes
    if other_labels:
        for lab in other_labels:
            if lab in texts_done:
                continue
            texts.append(
                plt.text(lfc[lab], fdr[lab], lab)
            )

    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

    ax.set_xlabel('Log$_2$ Fold Change')
    ax.set_ylabel('-log$_{10}$(FDR)')

    if outfn:
        plt.tight_layout()
        plt.savefig(outfn, dpi=150)
    return ax


revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])


def tabulate_mageck(prefix):
    """
    :param prefix: Input file prefix, including path
    :return: pd.DataFrame
    """
    prefix = Path(prefix)
    tables = {}
    tab = None
    for fn in os.listdir(prefix.parent):
        if not fn.endswith('.gene_summary.txt') or \
                prefix.parts[-1] not in fn:
            continue
        #mtab from the mageck output, reformatted into tab
        mtab = pd.read_csv(prefix.parent / fn, '\t', index_col=0)
        tab = pd.DataFrame(index=mtab.index)
        tab.loc[:, 'lfc'] = mtab.loc[:, 'neg|lfc']
        # turn sep pos|neg columns into one giving only the appropriate LFC/FDR
        pos = mtab['pos|lfc'] > 0
        tab.loc[pos, 'fdr'] = mtab.loc[pos, 'pos|fdr']
        tab.loc[~pos, 'fdr'] = mtab.loc[~pos, 'neg|fdr']
        tab.loc[:, 'fdr_log10'] = tab.fdr.apply(lambda x: -np.log10(x))

        sampnm = fn.split(prefix.stem)[1].split('.gene_s')[0]
        tables[sampnm] = tab
    if tab is None:
        raise FileNotFoundError('Failed to find any .gene_summary.txt files with prefix '+str(prefix))
    tbcolumns = pd.MultiIndex.from_product([sorted(tables.keys()), ['lfc', 'fdr', 'fdr_log10']],
                                           1)
    table = pd.DataFrame(index=tab.index, columns=tbcolumns)
    for exp, tab in tables.items():
        table[exp] = tab
    return table
