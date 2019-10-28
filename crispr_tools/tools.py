import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    def adjust_text(*args, **kwargs):
        pass
from typing import List, Union, Dict, Collection, Iterable
from pathlib import Path, PosixPath, WindowsPath
import pandas as pd
#import multiprocessing as mp

#hart_path = pkg_resources.resource_filename('crispr_tools', 'files/Hart2017_TableS2_core_genes.txt')
#if os.path.isfile(hart_path):
#    hart_list = pd.read_csv(hart_path, '\t', header=None)[0].values
#    hart_list = pd.Series(hart_list, index=hart_list)
#else:
#    print('Hart list not found')
from crispr_tools import hart_list

__version__ = 'v1.4.0'

# import logging
# slowvolcLOG = logging.getLogger('slow_volc')
# slowvolcLOG.setLevel(logging.DEBUG)

#v1.3.0
# Adding Mageck tab
# fixed plot volcanos enrichment labeling


#todo; pass plot_volcano a filen string and it loads the table
# todo standardised tables for all results that allow you to pull out x, y wiht same sample and score keys
import pkg_resources

def list_not_str(thing):
    """If thing is a string, wrap it in a list."""
    if type(thing) is str:
        return [thing]
    return thing


def drop_nonumeric(tab):
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )
    return tab

def split_nonumeric(tab):
    """drop columns from DF (tab) where [0, col] type is string,
    and return as a tuple with dropped columns as a sperate DF"""
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    texttab = tab.loc[:, nonumeric].copy()
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )

    return tab, texttab

def size_factor_normalise(cnt_tab, log=True):
    """The number by which MAGeCK uses to normalise its reads.
    Returns factor for each column.

    Uses +1 counts, returns a copy of the table.
    log=True returns np.log2 of the counts."""

    cnt_tab = cnt_tab.copy()
    out_tab = cnt_tab.copy()
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
    out_tab.loc[:, norm_tab.columns] = norm_tab
    return out_tab


def plot_read_violins(tab, rows=1, column_labs=None, log=True, size_norm=False, ax=None):
    """Takes counts table, does log2 (optionally size factor normalises)
    and produces violin density plots.

    figsize is (3*n, 6)

    Returns a plt.Axes

    High nrows sometimes results in empty rows."""

    samp_per_row = tab.shape[1] // rows
    # if the samp per row does not divide cleanly into the nsamps
    #   an extra row is required for the remainders
    if tab.shape[0] % samp_per_row:
        samp_per_row += 1

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
        stab = tab.iloc[:, inds].T
        # this is poo
        if stab.shape[0] == 0:
            break
        ax.violinplot(stab, widths=0.8)
        ax.boxplot(stab, widths=0.2)
        ax.set_xticks(range(1, samp_per_row + 1))
        ax.set_xticklabels(column_labs[inds], rotation=40)
    plt.tight_layout()
    if rows == 1:
        return axes[0]
    else:
        return axes


def plot_ROC(tab, things_oi, label=None, ax = None):
    #todo: make things_oi work as any kind of iterable (currently pd.Series)
    """Tab needs to be a series or dataframe containing values used to
    order the tab index."""
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
                 exclude_labs = ('NonT', 'Rando'), plot_kw:dict=None):
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
    #print('TP ', POINT)

    #this is silly
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
    #todo update to use scatter and set minimum matplotlib version
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
                plt.text(lfc[lab], fdr[lab], lab)
            )

    if texts and adjust_text:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

    if lfc_lab is None:
        ax.set_xlabel('Log$_2$ Fold Change')
    else:
        ax.set_xlabel(lfc_lab)

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
        # select correct files, the '._' are mac files that get caught on unix filesystem
        if not fn.endswith('.gene_summary.txt') or \
                prefix.parts[-1] not in fn or \
                fn.startswith('._'):
            continue
        #mtab from the mageck output, reformatted into tab
        mtab = pd.read_csv(prefix.parent / fn, '\t', index_col=0)
        tab = pd.DataFrame(index=mtab.index)

        # keep the pos and neg p values with better names
        tab.loc[:, 'neg_p'] = mtab['neg|p-value']
        tab.loc[:, 'pos_p'] = mtab['pos|p-value']

        tab.loc[:, 'lfc'] = mtab.loc[:, 'neg|lfc']
        # turn sep pos|neg columns into one giving only the appropriate LFC/FDR
        pos = mtab['pos|lfc'] > 0
        for stat in 'fdr', 'p-value':
            tabk = stat.replace('-value', '')
            tab.loc[pos, tabk] = mtab.loc[pos, f'pos|{stat}']
            tab.loc[~pos, tabk] = mtab.loc[~pos, f'neg|{stat}']
            tab.loc[:, f'{tabk}_log10'] = tab[tabk].apply(lambda x: -np.log10(x))

        sampnm = fn.split(prefix.stem)[1].split('.gene_s')[0]
        tables[sampnm] = tab
    if tab is None:
        raise FileNotFoundError('Failed to find any .gene_summary.txt files with prefix '+str(prefix))

    # create multiindex using the sample names and stats
    tbcolumns = pd.MultiIndex.from_product(
        [sorted(tables.keys()), ['lfc', 'fdr', 'fdr_log10', 'p', 'p_log10', 'pos_p', 'neg_p']],
        1
    )
    table = pd.DataFrame(index=tab.index, columns=tbcolumns)
    for exp, tab in tables.items():
        table[exp] = tab
    return table


def pca_grid(pca, hue_deet, style_deet, max_components=5, also_return_fig=False):
    """"""
    thing = sns.scatterplot(pca.components_[0], pca.components_[0],
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
            sns.scatterplot(pca.components_[pc_b], pca.components_[pc_a],
                            hue=hue_deet, style=style_deet,
                            s=150, legend=False)
            plt.xlabel(f"PC {pc_b+1} (variance explained: {pca.explained_variance_ratio_[pc_b]*100:.3}%)")
            plt.ylabel(f"PC {pc_a+1} (variance explained: {pca.explained_variance_ratio_[pc_a]*100:.3}%)")
    plt.tight_layout()
    if also_return_fig:
        return (fig, axes)
    else:
        return axes


def write_repmap(sample_reps:Dict[str, list], ctrlmap:Dict[str, list], repmap_fn:str):
    """Write a JACKS format replicate file.

    Args:
        sample_reps: maps given sample names to replicate names in the count file.
        ctrlmap: specifies which samples are controls for which other samples.
        repmap_fn: path to file that will be written.
    Returns:
        None
        """
    reptab = [['rep', 'samp', 'ctrl']]
    for ctrl_samp, treat_samps in ctrlmap.items():
        for ctrl_rep in sample_reps[ctrl_samp]:
            reptab.append([ctrl_rep, ctrl_samp, ctrl_samp])
        for treat in treat_samps:
            for rep in sample_reps[treat]:
                reptab.append([rep, treat, ctrl_samp])
                # print([rep, samp, ctrl])

    with open(repmap_fn, 'w') as f:
        for line in reptab:
            f.write('\t'.join(line) + '\n')