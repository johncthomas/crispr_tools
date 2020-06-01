mport numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from typing import Union, List, Dict, Tuple
#import pickle
from scipy import stats
from scipy.stats import norm
import statsmodels.api as sm


from statsmodels.stats.multitest import multipletests
try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    def adjust_text(*args, **kwargs):
        pass
from copy import copy
__version__ = '0.7.1'

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
    plt.plot([0, 0], [minlim, maxlim], 'g--', zorder=0)
    plt.plot([minlim, maxlim], [0, 0], 'g--', zorder=0)

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


def tabulate_score(prefix, return_ps=False):
    """Return a multiindexed DF of JACKS results.

    DF has multiindex, samples on level zero.

    Level one, if return_ps is false: 'jacks_score', 'fdr_log10', 'stdev'
    if True 'p' and 'p_log10' also included. This is for legacy reasons
    no reason to not get p values."""
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
        sig_df.loc[:, (exp, 'p_neg')] = ps[exp]
        sig_df.loc[:, (exp, 'p_pos')] = 1-ps[exp]
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

        for fdri, fdr in enumerate(['fdr_pos', 'fdr_neg']):
            score_table.loc[score_table[fdr] == 0, fdr] = (min_pos, min_neg)[fdri] / 10

        for mask, posneg in (pos,'pos'), (neg, 'neg'):
            score_table.loc[mask, 'fdr'] = score_table.loc[mask, 'fdr_'+posneg]
            score_table.loc[mask, 'p'] = score_table.loc[mask, 'p_' + posneg]

        sig_df.loc[:, (exp, 'fdr_log10')] = score_table['fdr'].apply(lambda x: -np.log10(x))
        sig_df.loc[:, (exp, 'p_log10')] = score_table['p'].apply(lambda x: -np.log10(x))
        # put the columns in a good order
        if not return_ps:
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


def mahalanobis_distance(table, xsample, ysample, line_ig = (0,1), **cp_kwargs):
    ols = FakeOLS(*line_ig)
    xs = table[xsample]['jacks_score']
    ys = table[ysample]['jacks_score']
    xs_sd = table[xsample]['stdev']
    ys_sd = table[ysample]['stdev']
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

    # return negative distance for the depleted genes
    depleted = euc_dist < 0
    mahal_dist[depleted] = 0 - mahal_dist[depleted]

    return pd.DataFrame({'Euclidean_dist':euc_dist, 'Mahal_dist':pd.Series(mahal_dist)})


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


def get_jacks_stats(df, sampx, sampy, negative_sig=True):
    # duplicated from crispr_screen_viewer.scatter
    """returns DF of  p-value, FDR, Difference, calculated between X
    and Y distributions. Significance values are -log10

    Args:
        df: a multiindex DF, levels[0] are sample names, levels[1] including
            jacks_score and stdev.
        sampx: name of column giving the stats of the sample to be plotted on x
        sampy: as above
        negative_sig: if True negative values returned for genes with neg diff

    Returns:
        a dataframe with columns: 'p-value', 'FDR', 'Difference'
        """
    # if you add more, update STAT_KEYS
    X = df[sampx]
    Y = df[sampy]
    # We can assume the subcol names i guess
    # todo make it return non-absolute diff and non-negative sig values by default, and update the charts
    diff = abs(Y.jacks_score - X.jacks_score)
    sumstd = np.sqrt(X.stdev**2 + Y.stdev**2)
    # get bool array, index will be
    negative = (Y.jacks_score < X.jacks_score).values
    # cdf value at zero
    ps = stats.norm(diff.abs(), sumstd).cdf(0)
    fdr = sm.stats.multipletests(ps, 0.1, 'fdr_bh')[1]
    ps, fdr = [-np.log10(s) for s in (ps, fdr)]
    if negative_sig:
        ps[negative] = 0-ps[negative]
        fdr[negative] = 0-fdr[negative]

    DISTANCES = {
        'p-value': ps,
        'FDR': fdr,
        'Difference': diff
    }

    return pd.DataFrame(DISTANCES)


def get_t(row, ak, bk, score_k='jacks_score', sterr_k='stdev'):
    """T-statistic for A>B. Written for use with DF.apply(get_t, axis=1, ...),
    but could work with properly structured obj, a dict of dict for example.
    Uses mean and stdev parameters to calculate t-stats.

    Args:
        row: Object with keys ak & bk, values are obj with keys score_k & sterr_k
        ak: key (column name) for sample A.
        bk: key (column name) for sample B.
        score_k: the key to access mean parameter.
        sterr_k: the key to access standard error param."""
    # t-statistic:
    #   abs(x-y)/(x.std+y.std)

    A, B = row[ak], row[bk]
    diff = A[score_k] - B[score_k]
    return diff / np.sqrt(A[sterr_k]**2 + B[sterr_k]**2)


def bootstrap_significance_from_control_genes(
        jacksRes, ctrlGeneMask, ctrl_k, treat_k,
        score_k='jacks_score', sterr_k='stdev'):
    """Get the significance using t-statistics of control genes (ideally
    resampled control guides to produce lots of genes). Minimum p will be
    1/ctrlGeneMask.sum().

    Returns:
        dict with keys "[enriched/depleted]_[p/fdr]" and values pd.Series
        with index containing only non-control genes

    Args:
        jacksRes: pd.DataFrame with multiindex columns
                (ctrl_k & treat_k, score_k & sterr_k)
        ctrlGeneMask: bool mask with same row index as jacksRes, where True
                specifies a control gene.
        ctrl_k, treat_k: sample names in jacksRes. Enrichment is defined as higher value
                in treatment compared to control.
        score_k, sterr_k: keys to access mean and stdev parameters to calculate T-statistics.

    """
    # test if control is higher
    t_depletion = jacksRes.apply(get_t, axis=1, ak=ctrl_k, bk=treat_k, score_k=score_k,sterr_k=sterr_k)
    # test treat
    t_enrichment = jacksRes.apply(get_t, axis=1, ak=treat_k, bk=ctrl_k, score_k=score_k,sterr_k=sterr_k )

    nFake = ctrlGeneMask.sum()
    res = {}
    for ts, label in (t_enrichment, 'enrich'), (t_depletion, 'deplete'):
        ctrl_ts = ts.loc[ctrlGeneMask]
        # get proportion of fake gene ts that are greater than each t
        res[label] = ts.apply(lambda x: (ctrl_ts > x).sum() / nFake)

    enrich_p = res['enrich'].loc[~ctrlGeneMask]
    deplet_p = res['deplete'].loc[~ctrlGeneMask]

    # set minimum p value, as we only know that p < 1/nFake when p == 0
    for ps in enrich_p, deplet_p:
        ps.loc[ps == 0] = 1/nFake

    enrich_fdr = pd.Series(
        sm.stats.multipletests(enrich_p, method='fdr_bh')[1],
        index=jacksRes.loc[~ctrlGeneMask].index
    )

    deplete_fdr = pd.Series(
        sm.stats.multipletests(deplet_p, method='fdr_bh')[1],
        index=jacksRes.loc[~ctrlGeneMask].index
    )
    return {'enriched_fdr': enrich_fdr, 'depleted_fdr': deplete_fdr,
            'enriched_p': enrich_p, 'depleted_p': deplet_p, }

try:
    from jacks.jacks_io import *
    def run_pseudo_genes(countfile, replicatefile, guidemappingfile,
                  rep_hdr=REP_HDR_DEFAULT, sample_hdr=SAMPLE_HDR_DEFAULT, common_ctrl_sample=COMMON_CTRL_SAMPLE_DEFAULT,
                  ctrl_sample_hdr=None, sgrna_hdr=SGRNA_HDR_DEFAULT, gene_hdr=GENE_HDR_DEFAULT,
                  apply_w_hp=APPLY_W_HP_DEFAULT, norm_type=NORM_TYPE_DEFAULT,
                  ignore_blank_genes=False, ctrl_genes=None, reffile=None, n_pseudo=0, count_prior=32):

        sample_spec, ctrl_spec, gene_spec, x_ref = preprocess(countfile, replicatefile, guidemappingfile,
                                                                    rep_hdr, sample_hdr, common_ctrl_sample,
                                                                    ctrl_sample_hdr, sgrna_hdr, gene_hdr,
                                                                    ignore_blank_genes, reffile)

        ctrl_geneset = ctrl_genes

        data, meta, sample_ids, genes, gene_index = loadDataAndPreprocess(sample_spec, gene_spec, ctrl_spec=ctrl_spec,
                                                                          normtype=norm_type, ctrl_geneset=ctrl_geneset,
                                                                          prior=count_prior)


        testdata, ctrldata, test_sample_idxs = collateTestControlSamples(data, sample_ids, ctrl_spec)


        # Add a set of pseudo genes, created by randomly sampling from guides targeting genes in the control set
        LOG.info('Running JACKS inference on %d pseudogenes' % n_pseudo)
        pseudo_gene_index = createPseudoNonessGenes(gene_index, ctrl_geneset, n_pseudo)
        # structure of results: results[gene] = (y, tau, x1, x2, w1, w2)
        jacks_pseudo_results = inferJACKS(pseudo_gene_index, testdata, ctrldata, apply_w_hp=apply_w_hp)

        # get W1 values
        pseudo_scores = [(np.nanmean(jacks_pseudo_results[gene][4]), gene) for gene in jacks_pseudo_results]

        return pseudo_scores


except ImportError:
    pass

