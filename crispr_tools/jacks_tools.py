import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from typing import Union
#import pickle
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    def adjust_text(*args, **kwargs):
        pass

__version__ = '0.6'

#todo: sepfucntions for each table, one function that calls all, option to only return scores
#todo subclass scores_table from DF, add functions as methods, plots, mahal, scores etc
#todo: put eff and fold change on the same table (currently nans)
#todo pass ax to plot2d
#todo bootstrapping?
#todo: write decent readme and documention for github
#todo make scores scatterplot work with other distance methods

# in 0.5
# refactor "esstab" to "score_table"

def plot_volcano_from_scoretable(score_table, savefn=None, ax=None,
                 label_deplet = 0, label_enrich = 0, other_labels=None,
                 p_thresh = 0.05):
    """Supply pandas table with 'jacks_score' and 'fdr_pos/neg' columns.
    Returns fig, ax if no ax is supplied."""
    if ax is None:
        fig, ax_ = plt.subplots(1, 1, figsize=(10, 10))
    else:
        plt.sca(ax)

    score_table = score_table.copy()

    pos = score_table['jacks_score'] > 0
    neg = ~pos

    # get lowest not zero and set zeroes to 1/10 that value
    min_pos = min(score_table.loc[score_table['fdr_pos'] > 0, 'fdr_pos'])
    min_neg = min(score_table.loc[score_table['fdr_neg'] > 0, 'fdr_neg'])
    #print(min_neg, min_pos)
    for fdri, fdr in enumerate(['fdr_pos', 'fdr_neg']):
        score_table.loc[score_table[fdr] == 0, fdr] = (min_pos, min_neg)[fdri] / 10

    for mask, fdr in (pos, 'fdr_pos'), (neg, 'fdr_neg'):
        score_table.loc[mask, 'fdr'] = score_table.loc[mask, fdr]

    score_table.loc[:, 'fdr'] = score_table['fdr'].apply(lambda x: -np.log10(x))

    faces = dict(facecolors='none', edgecolors='b')
    # plt.scatter(score_table.loc[pos, 'jacks_score'], score_table.loc[pos, 'fdr_pos'], **faces)
    # plt.scatter(score_table.loc[neg, 'jacks_score'], score_table.loc[neg, 'fdr_neg'], **faces)
    plt.scatter(score_table.jacks_score, score_table.fdr, **faces)

    p = -np.log10(p_thresh)

    plt.plot([min(score_table['jacks_score']), max(score_table['jacks_score'])],
             [p, p], 'k--')

    # label the top and bottom most, if specified
    texts = []
    texts_done = []
    # get subtables
    dep = score_table.sort_values('jacks_score').head(label_deplet)
    enr = score_table.sort_values('jacks_score').tail(label_enrich)
    for stab in dep, enr:
        for lab, row in stab.iterrows():
            if row.fdr < p:
                continue
            texts_done.append(lab)
            texts.append(plt.text(row.jacks_score, row.fdr, lab))
    # label additional genes
    if other_labels:
        for lab, row in score_table.loc[other_labels, :].iterrows():
            if lab in texts_done:
                continue
            texts.append(
                plt.text(score_table.loc[lab, 'jacks_score'],
                         score_table.loc[lab, 'fdr'],
                         lab)
            )
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

    plt.xlabel('JACKS score')
    plt.ylabel('-log10(FDR)')

    if savefn :
        plt.savefig(savefn)
    if ax is None:
        return fig, ax_
    else:
        return ax


def scores_scatterplot(x, y, table=None, min_mahal:float=False,
                       label_mahal:Union[int, bool]=False,
                        labels=None, formatters=None, ):
    """Produce biplot of 2 essentiality series.
    args:
        x, y:
            str that point to column labels, or
            pd.Series with shared index that gives x and y values
            to be plotted.

        table:
            Score table from e.g. .tabulate_scores() (optional)

        labels:
            list of gene names to be labeled on the plot

        min_mahal:
            minimum distance to colour points

        label_mahal:
            True to label all points > min_mahal, pass an int to label the
            top N furthest points.

        formatters:
            list of tuples as below. format_dict is passed to
            plt.plot(..., **format_dict).
                [(list_of_genes, format_dict), ...]

    Formatters"""

    fig, ax = plt.subplots(1, 1, figsize=(8,8))

    if table is None:
        x_score = x
        y_score = y
    else:
        x_score = table[x].jacks_score
        y_score = table[y].jacks_score

    if min_mahal is not False:
        _, _, mahal = mahal_nocov(
            FakeOLS(1,0), x_score, table[x].stdev, y_score, table[y].stdev
        )

        # for this we don't care about pos/neg
        mahal = mahal.abs()
        max_mahal = max(mahal)
        #print(max_mahal)
        # normalised mahaladonis for color scaling
        # Anything > 0 will go on labels_maybe
        nmahal = (mahal-min_mahal)/(max_mahal-min_mahal)
        colrs = cm.viridis(nmahal) # loses the gene name indicies
        labels_maybe = []
        # RGBA
        grey = (0.68, 0.68, 0.68, 1)
        for i, m in enumerate(nmahal):
            # m == mahal-minmahal
            if m < 0:
                colrs[i] = grey
            else:
                labels_maybe.append(nmahal.index[i])
    else:
        colrs = None

    if label_mahal:
        if type(label_mahal) is int:
            topdist = mahal.sort_values().tail(label_mahal)
            labels_maybe = [l for l in labels_maybe if l in topdist.index]
        if labels:
            labels.extend(labels_maybe)
        else:
            labels = labels_maybe

    if not formatters:
        # all genes with empty format spec
        formatters = [(x_score.index, {})]
    for genes, formats in formatters:
        if 'marker' not in formats:
            formats['marker'] ='o'
        plt.scatter(x_score.loc[genes], y_score.loc[genes], c=colrs, **formats)

    # plot a square Axes with orientation lines
    lim = ax.axis()
    minlim = min(lim[0], lim[2])
    maxlim = max(lim[1], lim[3])
    plt.plot([minlim, maxlim], [minlim, maxlim], 'k--', alpha=0.2)
    plt.plot([0, 0], [minlim, maxlim], 'g--')
    plt.plot([minlim, maxlim], [0, 0], 'g--')

    if colrs is not None:
        sm = plt.cm.ScalarMappable(cmap=cm.get_cmap('viridis'))
        sm.set_array([])
        cb = plt.colorbar(sm, fraction=0.03, pad=0.01, aspect=5)

        cb.set_ticks([0.0, 1.0])
        cb.set_ticklabels([str(round(min_mahal, 2)),
                           str(round(max_mahal, 2))])

    # labels
    txt = []
    # label significants and control
    if labels:
        for lab in set(labels):
            txt.append(plt.text(x_score[lab], y_score[lab], lab))

    # label axes
    if table is not None:
        plt.xlabel(x+' JACKS score')
        plt.ylabel(y + ' JACKS score')

    adjust_text(txt, arrowprops=dict(arrowstyle='->'))


    return fig, ax


def tablulate_score(prefix):
    """Return a multiindexed DF of JACKS results.

    Table columns are sample names as given in the repmap at level 0,
    and then 'jacks_score, fdr_log10, fdr_neg, fdr_pos, stdev' at level 1"""
    # othertab = pd.DataFrame(columns=("IC10","IC90","D14"), index=essen['D14'].index)
    # Tables produced by jacks have columns that are the groups
    genes = pd.read_table(prefix + '_gene_JACKS_results.txt', sep='\t', index_col=0)
    genes_index = sorted(genes.index)
    genes = genes.reindex(genes_index)
    genesstd = pd.read_table(prefix + '_gene_std_JACKS_results.txt', sep='\t', index_col=0)
    genesstd = genesstd.reindex(genes_index)
    ps = genes / genesstd
    ps = ps.apply(norm.cdf)

    # multiindex DF for each experiment giving results
    sig_cols = pd.MultiIndex.from_product((ps.columns, ['jacks_score', 'fdr_pos', 'fdr_neg', 'fdr_log10', 'stdev']),
                                          names=('exp', 'stat'))
    sig_df = pd.DataFrame(index=genes_index, columns=sig_cols)

    for exp in ps.columns:
        sig_df.loc[:, (exp, 'fdr_neg')] = multipletests(ps[exp], method='fdr_bh')[1]
        sig_df.loc[:, (exp, 'fdr_pos')] = multipletests(1 - ps[exp], method='fdr_bh')[1]
        sig_df.loc[:, (exp, 'jacks_score')] = genes[exp]
        sig_df.loc[:, (exp, 'stdev')] = genesstd[exp]
        score_table = sig_df[exp].copy()
        # get one FDR
        pos = score_table['jacks_score'] > 0
        neg = ~pos
        # get lowest not zero and set zeroes to 1/10 that value
        min_pos = min(score_table.loc[score_table['fdr_pos'] > 0, 'fdr_pos'])
        min_neg = min(score_table.loc[score_table['fdr_neg'] > 0, 'fdr_neg'])
        #print(min_neg, min_pos)
        for fdri, fdr in enumerate(['fdr_pos', 'fdr_neg']):
            score_table.loc[score_table[fdr] == 0, fdr] = (min_pos, min_neg)[fdri] / 10

        for mask, fdr in (pos, 'fdr_pos'), (neg, 'fdr_neg'):
            score_table.loc[mask, 'fdr_log10'] = score_table.loc[mask, fdr]

        sig_df.loc[:, (exp, 'fdr_log10')] = score_table['fdr_log10'].apply(lambda x: -np.log10(x))
        # put the columns in a good order
        sig_df = sig_df.reindex(['jacks_score', 'fdr_log10', 'stdev'],
                       axis=1, level=1, )
    return sig_df

def tabulate_jacks(prefix):
    """Return 3 tables giving the:
        1. jacks_score|fdr_pos|fdr_neg|std,
        2. guide efficacy data,
    and 3. fold changes for each gene.

    Tables are multiindexed by sample name and then results columns for those
    samples.

    fdr in the

    Prefix is the used to identify the results files. So prefix
    should contain the path to the files if they aren't in os.getcwd()"""

    kwtab = dict(sep='\t', index_col=0)

    sig_df = tablulate_score(prefix)
    samples = sig_df.columns.levels[0]
    # get guide data, foldchange and efficacies
    guide_cols = pd.MultiIndex.from_product((samples, ['foldchange', 'fold_std', 'eff', 'eff_std']),
                                            names=['exp', 'stat'])
    fchange_df = pd.DataFrame(columns=guide_cols)
    foldchange = pd.read_table(prefix + '_logfoldchange_means.txt', **kwtab)
    foldstd = pd.read_table(prefix + '_logfoldchange_std.txt', **kwtab)
    eff_tab = pd.read_table(prefix + '_grna_JACKS_results.txt', **kwtab)

    for exp in samples:
        fchange_df.loc[:, (exp, 'lfc')] = foldchange[exp]
        fchange_df.loc[:, (exp, 'fold_std')] = foldstd[exp]
    fchange_df.loc[:, 'gene'] = foldchange['gene']

    efficacies = pd.DataFrame(columns=('eff', 'eff_std'))
    efficacies.loc[:, 'eff'] = eff_tab['X1']
    efficacies.loc[:, 'eff_std'] = (eff_tab['X2'] - eff_tab['X1'] ** 2) ** 0.5
    efficacies.loc[:, 'gene'] = fchange_df['gene']

    return sig_df, efficacies, fchange_df


def dist_to_line(ols, p):
    """ols: sm.OLS results obj
    p: a tuple giving xy.
    doesn't expect endogenous constant"""
    p = np.array(p)
    x = p[0]
    y = p[1]
    # get the points that define the line
    a = np.array([x, ols.predict(x)])
    b = np.array([y, ols.predict(y)])

    d = np.cross(b - a, p - a) / np.linalg.norm(b - a)

    return d


class FakeOLS:
    def __init__(self, slope, intcpt):
        self.slope = slope
        self.intcpt = intcpt

    def predict(self, x):
        # assert len(x) == 2
        # x = x[1]
        return x * self.slope + self.intcpt


def closest_point(ols, xs, ys, verbose=False):
    from shapely.geometry import LineString, Point
    # xs and ys are the gene essentialities for ctrl & exp
    # ols is the statsmodel.OLS object giving the trend line for xs&ys
    # FakeOLS can be used to compare to x=y slope
    closest_x = []
    closest_y = []
    mins = min(min(xs), min(ys))
    maxs = max(max(xs), max(ys))

    for x, y in zip(xs, ys):
        line = LineString([(mins, ols.predict(mins)), (maxs, ols.predict(maxs))])
        p = Point(x, y)
        pj = line.project(p)
        cp = line.interpolate(pj)
        closest_x.append(cp.x)
        closest_y.append(cp.y)
        if verbose:
            print(p, '->', cp)
    return np.array(closest_x), np.array(closest_y)


def mahal_nocov(ols, xs, xs_sd, ys, ys_sd, **cp_kwargs):

    # get closest points, calc significance
    cxs, cys = closest_point(ols, xs, ys, **cp_kwargs)
    mahal_dist = np.sqrt((cxs - xs) ** 2 / xs_sd + (cys - ys) ** 2 / ys_sd)
    # p: a tuple giving xy. a, b, start and stop of line
    minx = min(min(xs), min(ys))
    maxx = max(max(xs), max(ys))
    a, b = [(minx, ols.predict(minx)), (maxx, ols.predict(maxx))]
    a, b = np.array(a), np.array(b)
    points = [np.array(p) for p in zip(xs, ys)]
    euc_dist = [np.cross(b - a, p - a) / np.linalg.norm(b - a) for p in points]
    euc_dist = np.array(euc_dist)
    # get significance
    ps = 1 - pd.Series(mahal_dist).apply(norm.cdf)

    # return negative distance for the depleted genes
    depleted = euc_dist < 0
    mahal_dist[depleted] = 0 - mahal_dist[depleted]

    return euc_dist, ps, pd.Series(mahal_dist)





