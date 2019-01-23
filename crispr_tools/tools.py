import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

__version__ = 'v1.2.1'

#todo: sometimes lfc tables not writte, make it optional
#todo: boxpot on the violins
#todo: pass through non numeric columns in SFN
#todo; pass plot_volcano a filen string and it loads the table
#todo make plot volcano work with both

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


def plot_volcano(tab, title='', label_genes=None, outfn='', ax=None, source='mageck'):
    """Take a mageck table and do a volcano plot. """

    dep = tab.loc[:, 'pos|lfc'] > 0
    enr = tab.loc[:, 'pos|lfc'] < 0
    tab.loc[dep, 'ml10_fdr'] = -np.log10(tab.loc[dep, 'pos|fdr'])
    tab.loc[enr, 'ml10_fdr'] = -np.log10(tab.loc[enr, 'neg|fdr'])
    sctkw = dict(marker='o', linestyle='none')
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 12))
    ax.plot(tab['pos|lfc'], tab['ml10_fdr'], alpha=0.4, **sctkw)
    # plt.yscale('log')
    # plt.gca().invert_yaxis()
    ax.set_title(title)

    xmin, xmax, ymin, ymax = ax.axis()
    ax.plot([xmin, xmax], [-np.log10(0.05), -np.log10(0.05)], 'k--')
    ax.plot([0, 0], [ymin, ymax], color='silver')
    if label_genes is not None:

        if type(label_genes) is int:
            ngoi = label_genes
            label_genes = list(tab.sort_values(['neg|fdr', 'neg|lfc']).index[:ngoi]) \
                          + list(tab.sort_values(['pos|fdr', 'pos|lfc'], ascending=[True, False]).index[:ngoi])
            label_genes = [gn for gn in label_genes if tab.loc[gn, 'neg|fdr'] < 0.05 or tab.loc[gn, 'pos|fdr'] < 0.05]

        texts = []
        if adjust_text is not None:
            for gene in label_genes:
                goix, goiy = tab.loc[gene, 'pos|lfc'], tab.loc[gene, 'ml10_fdr']
                tx = ax.text(goix, goiy, gene)
                ax.plot(goix, goiy, color='orange', **sctkw)
                texts.append(tx)
            adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

    ax.set_xlabel('Log$_2$ Fold Change')
    ax.set_ylabel('-log$_{10}$(FDR)')

    if not outfn:
        return ax
    if outfn:
        plt.savefig(outfn)

revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])

