import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from typing import Union, List, Dict, Tuple
#import pickle
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    def adjust_text(*args, **kwargs):
        pass
from copy import copy
__version__ = '0.7'

#todo: sepfucntions for each table, one function that calls all, option to only return scores
#todo subclass scores_table from DF, add functions as methods, plots, mahal, scores etc
#todo: put eff and fold change on the same table (currently nans)
#todo bootstrapping?
#todo: write decent readme and documention for github

# in 0.5
# refactor "esstab" to "score_table"


def scores_scatterplot(x, y, table=None, distance_gradient=True, label_pos=0, label_neg=0,
                       distance:pd.Series=None, min_label_dist=0.5, min_label_diff=0.5,
                        labels=None, formatters=None, dist_name = None, gradient=1, ax=None) -> plt.Axes:
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

        mahal_gradient:
            if True points will be colored by mahal value.

        label_pos, label_neg:
            Pass an int or True to label the most distant top or bottom points.

        minlabel:
            minimum abs(mahal) for labels to be applied by label_pos/neg

        formatters:
            list of tuples as below. format_dict is passed to
            plt.plot(..., **format_dict).
                [(list_of_genes, format_dict), ...]

    returns fig, ax if ax=None, otherwise the supplied ax is returnd."""

    # required changes;
    #   distance optionally supplied
    #   optional filter by abs(x-y)>threshold (use minlabel...?)
    #   remove refs to jacks_score


    fig=None
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8,8))

    if table is None:
        x_score = x
        y_score = y
    else:
        cols = table[x].columns
        for k in ('jacks_score', 'lfc'):
            if k in cols:
                break
        else:
            raise ValueError('No score column found in '+str(cols))
        x_score = table[x][k]
        y_score = table[y][k]

    if labels is None:
        labels = []
    else:
        labels = copy(labels)

    # deal with any nans
    nans = (x_score.isna() | y_score.isna())
    x_score, y_score = x_score.loc[~nans], y_score.loc[~nans]

    pos_genes, neg_genes = [], []
    if distance_gradient:
        if distance is None:
            _, _, distance = mahal_nocov(
                x_score, table[x].stdev, y_score, table[y].stdev, line_ig=(0, gradient)
            )
        # no nans
        distance = distance.loc[~distance.isna()]
        # get the genes to be labelled first
        distance = distance.sort_values()

        # so that most distant points are plotted on top, order x and y by distance
        x_score = x_score[distance.index]
        y_score = y_score[distance.index]

        # get most distance in the pos and neg directions
        # get most distance in the pos and neg directions
        if label_pos is True:
            label_pos = np.inf
        if label_neg is True:
            label_neg = np.inf

        # get top N labels, filtered for those passing distance and x-y difference thresholds
        diff_mask = abs(x_score - y_score) > min_label_diff
        pos_genes = distance.loc[(distance > min_label_dist) & diff_mask].tail(label_pos).index
        neg_genes = distance.loc[(distance < 0 - min_label_dist) & diff_mask].head(label_neg).index


        # for gradient we don't care about pos/neg
        # but we do want the order to match x_score as the index will be lost
        distance = distance.reindex(x_score.index).abs()
        max_dist = max(distance)

        # normalised mahaladonis for color scaling
        if distance_gradient is True:
            min_dist = 0
        else:
            min_dist = distance_gradient
        norm_dist = (distance-min_dist)/(max_dist-min_dist)
        colrs = cm.viridis(norm_dist) # loses the gene name indicies
        # RGBA
        grey = (0.68, 0.68, 0.68, 1)
        for i, m in enumerate(norm_dist):
            # m == mahal-minmahal
            if m < 0:
                colrs[i] = grey
    else:
        colrs = None

    labels.extend(pos_genes)
    labels.extend(neg_genes)
    #print(labels)
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
        cb.set_ticklabels([str(round(min_dist, 2)),
                           str(round(max_dist, 2))])
        if dist_name is not None:
            cb.ax.set_ylabel(dist_name)

    # labels
    txt = []
    # label significants and control
    if labels:
        for lab in set(labels):
            #txt.append(plt.text(x_score[lab], y_score[lab], lab))
            txt.append(
                plt.annotate(
                    lab,
                    (x_score[lab], y_score[lab]),
                    #bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.3),
                    arrowprops=dict(arrowstyle='->', color='#ff8533' )
                )
            )

    # label axes
    if table is not None:
        plt.xlabel(x+' '+k)
        plt.ylabel(y +' '+k)

    avoid_points = True if distance_gradient is True else False
    if avoid_points:
        adjust_text(txt, x_score, y_score, force_points=(0.2,0.25))
    else:
        adjust_text(txt, expand_text=(1.1, 1.2), )

    if fig is not None:
        return ax
    else:
        return fig, ax


def tabulate_score(prefix):
    """Return a multiindexed DF of JACKS results.

    Table columns are sample names as given in the repmap at level 0,
    and then 'jacks_score, fdr_log10, fdr_neg, fdr_pos, stdev' at level 1"""
    # othertab = pd.DataFrame(columns=("IC10","IC90","D14"), index=essen['D14'].index)
    # Tables produced by jacks have columns that are the groups
    genes = pd.read_csv(prefix + '_gene_JACKS_results.txt', sep='\t', index_col=0)
    genes_index = sorted(genes.index)
    genes = genes.reindex(genes_index)
    genesstd = pd.read_csv(prefix + '_gene_std_JACKS_results.txt', sep='\t', index_col=0)
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
        sig_df = sig_df.reindex(['jacks_score', 'fdr_log10', 'stdev'], axis=1, level=1, )
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

    sig_df = tabulate_score(prefix)
    samples = sig_df.columns.levels[0]
    # get guide data, foldchange and efficacies
    guide_cols = pd.MultiIndex.from_product((samples, ['foldchange', 'fold_std', 'eff', 'eff_std']),
                                            names=['exp', 'stat'])
    fchange_df = pd.DataFrame(columns=guide_cols)
    foldchange = pd.read_csv(prefix + '_logfoldchange_means.txt', **kwtab)
    foldstd = pd.read_csv(prefix + '_logfoldchange_std.txt', **kwtab)
    eff_tab = pd.read_csv(prefix + '_grna_JACKS_results.txt', **kwtab)

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
    def __init__(self, intcpt, slope):
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



def mahal_nocov(xs, xs_sd, ys, ys_sd, line_ig = (0,1), **cp_kwargs):
    """Pass x, y values and stdev. Returns (euclidian distance, p-values of
    mahal distance, and mahal_distance from the x=y line.

    Other lines can be specified using line_gi giving gradient and intercept
    (y=ax+b)"""

    ols = FakeOLS(*line_ig)
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



# from crispr_tools.tools import tabulate_mageck
# import os
# os.chdir('/Users/johnc.thomas/Dropbox/crispr/screens_analysis')
# tab = tabulate_mageck('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/scattertst/tst.')
#
# scores_scatterplot('B1-vs-B0', 'B2-vs-B0', tab,True, 5,5, distance=tab['B1-vs-B2'].fdr_log10, min_label_diff=0)
# plt.show()


