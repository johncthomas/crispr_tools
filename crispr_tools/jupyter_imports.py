"""Imports things commonly used by the author. Might not be of any use to anyone else."""

import platform, pathlib
import os
from typing import Dict, List, Union, Tuple
import typing
import itertools
from itertools import combinations, combinations_with_replacement
import attrdictionary
import collections

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid')
import scipy.stats as stats
import yaml
from functools import partial
revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])
import logging


pltlogger = logging.getLogger('matplotlib')
pltlogger.setLevel(logging.WARNING)


from crispr_tools import tools, jacks_tools, crispr_pipeline, drugz
from IPython.display import display

from pathlib import Path

try:
    from jacks.jacks_io import runJACKS
    import jacks
    jacks.jacks_io.LOG.setLevel(logging.WARNING)
except:
    pass

from .tools import (
    hart_list,
    load_analyses_via_expd,
    write_stats_workbook,
    clonal_lfcs,
    load_counts,
    write_counts,
    abundance_normalise
)

from .plotting import pca_grid, plot_read_violins, plot_ROC

from sklearn.metrics import precision_recall_curve
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import statsmodels.api as sm

OLS = sm.regression.linear_model.OLS


csv_good_encoding = 'utf-8-sig'

def hxbin(x,y, **kwarghole):
    plt.hexbin(x,y, gridsize=40, bins='log')

is_olfactory = lambda x: x.startswith("OR") and x[2].isnumeric()

comp_name = platform.node()
if comp_name == 'POS':
    DROPBOX = '/mnt/d/Dropbox/'
elif platform.system() == 'Darwin':
    # Darwin is apparently the name of the mac OS
    DROPBOX = '/Users/johnc.thomas/Dropbox/'
else:
    DROPBOX = '/usr/johnc.thomas/Dropbox/'

if os.path.isdir(DROPBOX):
    DROPBOX = pathlib.Path(DROPBOX)
else:
    print(f"DROPBOX location not set, platform '{comp_name}' not recognised")
    del DROPBOX

#TODO move these non-crispr specifc tools to another package

log2p1 = lambda x: np.log2(x+1)

from crispr_tools.tools import list_not_str

ARROW = '→'

def multipletests_fdr(ps, method='fdr_bh', **kwargs):
    """Calls statsmodels.stats.multipletests and returns the corrected
    p-values only, with default of fdr_bh, rather than a FWER method.
    """
    kwargs['method'] = method
    try:
        if ps.isna().any():
            raise RuntimeError('ps contains NAN, this will break it')
    except AttributeError:
        pass
    qs = sm.stats.multipletests(ps, **kwargs)[1]
    if type(ps) == pd.core.series.Series:
        return pd.Series(qs, index=ps.index)
    return qs


# just putting this here as it's only for jupyter really
def plt_labels(xlab='', ylab='', title=''):
    """Add axes labels and figure title to the current matplotlib Figure"""
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)

class working_dir:
    """Context manager to temporarily alter the working directory.
    Use as part of a `with` statment.

    E.g:

    os.chdir('/path/one')
    with working_dir('/path/two/'):
        print(os.getcwd())
    print(os.getcwd())

    out:
        [1] /path/two
        [2] /path/one"""

    def __init__(self, new_dir):
        self.new_dir = new_dir
        self.old_dir = os.getcwd()
    def __enter__(self):
        os.chdir(self.new_dir)
        return os.getcwd()
    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.old_dir)
        if exc_type:
            return False
        return True


def write_excel_text(sheets: Dict[str, pd.DataFrame],
                     filename: str,
                     text_columns: Union[List, Dict[str, List[str]]] = 'ALL',
                     index_as_text=False,
                     **to_excel_kwargs):
    """Write XLSX with columns optionally formated as text.
    Args:
        sheets: Dataframes to be written as sheets, dict keyed by sheetname,
            values DF
        filename: For the written XLSX
        text_columns: Specify columns (by name) to be formated as text.
            Default 'ALL' formats all columns as string.
            Lists of column names (including index name) can specify columns.
            A single list applies to all sheets, use a dictionary to specify
            different columns per sheet. Can't specify the index
        index_as_text: Since you can't just pass the name of the index (for
            tedious reasons), set this to True to write the index as text.
            Index name will be written as header. Won't work with multiindex
            columns for some reason.


        **to_excel_kwargs: Additional kwargs are passed to pd.DataFrame.to_excel()

    Doesn't work work with multiindex columns.
    """
    # tests 'all', text_col is list or dict

    text_columns = list_not_str(text_columns)

    for df in sheets.values():
        if type(df.columns) == pd.MultiIndex:
            raise RuntimeError('MultiIndex columns not supported')

    writer = pd.ExcelWriter(filename,
                            engine='xlsxwriter')

    # Go through the DF
    for sheet_name, df in sheets.items():
        if index_as_text:
            df = df.reset_index(drop=False)
        df.to_excel(writer, sheet_name=sheet_name, **to_excel_kwargs)

        worksheet = writer.sheets[sheet_name]

        workbook = writer.book
        txt = workbook.add_format({'num_format': '@'})

        # get the columns to be textualised

        if text_columns == ['ALL']:
            txtcols = df.columns
        elif type(text_columns) is Dict:
            txtcols = text_columns[sheet_name]
            if txtcols == 'ALL':
                txtcols = df.columns
        elif (text_columns == None) or (text_columns == False):
            txtcols = []
        else:
            txtcols = list(text_columns)

        # Note: index formats get written over by pd, and any sensible automated
        #   solution to work around this uses index=False - but that's incompatible with
        #   multiindex columns for some fucking reason, so there is no good
        #   solution for dealing with the index until Pandas impliments
        #   a stable way of disabling automatic formating of the index

        index_offset = 1
        if 'index' in to_excel_kwargs:
            if not to_excel_kwargs['index']:
                index_offset = 0

        # get indicies of named columns
        col_i = []
        col_i.extend([list(df.columns).index(c) + index_offset for c in txtcols])

        # set the format of selected columns
        for i in col_i:
            worksheet.set_column(i, i, cell_format=txt)

    writer.save()
    return writer


def underscore_columns(df):
    """replace spaces with underscores in all columns of a pandas DF.
    In place, but returns the DF"""
    df.columns = df.columns.map(lambda x: x.replace(' ', '_'))
    return df

def hixbin(x,y, ax=None, **hxbn_kwargs):
    if ax is None:
        plt.figure(figsize=(4,4))
        ax = plt.gca()
    kwargs = dict({'bins':'log', 'gridsize':40}, **hxbn_kwargs)
    ax.hexbin(x,y, **kwargs)


def prep_data_for_upsetplot(bool_table):
    """UpSetPlot (for displaying intersection sizes) requires data to be
    formated in a very specific way. This does that.

    Expects a table of bools with columns being the samples"""

    bool_counts = {}

    # get table with each combination of bools, and their occurance
    for row in bool_table.values:
        # row = ''.join([str(int(b)) for b in row])
        row = tuple(row)
        try:
            bool_counts[row] += 1
        except KeyError:
            bool_counts[row] = 1

    b = (True, False)
    multindx = pd.MultiIndex.from_product([b for _ in bool_table.columns],
                                          names=bool_table.columns)

    # duplicates the data in 2 dimensions, and the sample names are on the wrong dimension
    usp_data = pd.DataFrame(bool_counts, index=multindx)
    # so transpose a row
    usp_data.iloc[:, 0] = usp_data.iloc[0, :]
    # and select that column
    usp_data = usp_data.iloc[:, 0]

    return usp_data


def dir2(obj):
    return [c for c in dir(obj) if c[0] != '_']


def do_plot_regression(x, y, ax, colour='black'):
    if all([type(xy) is pd.Series for xy in (x,y)]):
        if x.isna().any() or y.isna().any():
            raise RuntimeError('x and/or y contain NANs.')

    ols = OLS(y, sm.add_constant(x)).fit()
    xs = [x.min() * 0.8, x.max() * 0.8]
    ys = ols.predict(sm.add_constant(xs))

    ax.plot(xs, ys, color=colour, )


def why_FileNotFound(p):
    """Prints the left-most part of a file path that doesn't exist"""
    p = os.path.normpath(p)
    if os.path.isfile(p):
        print(p, 'is a file')
        return p
    elif os.path.isdir(p):
        print(p, 'is a dir')
        return p
    else:
        target = ''
        p = p.split('/')

        for bit in p:
            if not os.path.isdir(target+'/'):
                print(target, 'is not found')
                return
            target = target+'/'+bit


def index_of_true(S:pd.Series) -> pd.Series:
    return S[S].index

def keys_crawler(d, level=0):
    """Accepts a dictionary, gets the type of values for each key,
    if it's a list the type of the FIRST item obtained. Lists and dicts
    encountered are explored recursively.

    Output gives the structure of the object as a string, levels of
    indentation show the hierarchy. Lists indicated with [...]"""
    out = ''
    for k in d.keys():
        out += '  '*level+k+':\n'
        if type(d[k]) is dict:
            out += keys_crawler(d[k], level=level+1)
        elif type(d[k]) is list:
            out += '  '*(level+1)+'[\n'
            if type(d[k][0]) is dict:
                out += keys_crawler(d[k][0], level=level+2)
            out += '  '*(level+1)+']\n'
        else:# type(d[k]) is str:
            out = out[:-1] + f'  {type(d[k])}\n'
    return out

def minminmaxmax(x,y):
    nn = min([min(x), min(y)])
    xx = max([max(x), max(y)])
    return (nn, xx)


def nlprint(things:typing.Collection[str], sort=False):
    if sorted:
        things = sorted(things)
    print('\n'.join(things))

data = np.random.randn(100)



width = 0.8     # the maximum width of each 'row' in the scatter plot


def jitter_density_plot(
        data: Union[pd.DataFrame, Dict],
        width=0.8,
        equal_width=True,
        ax=None,
        tick_labels=None,
        per_column_kwargs: List[Dict] = None,
        **scatter_kwargs):
    # The other thing I could do here is generate just a list of jittered x values
    """Jitter plots scaled with width scaled to density.
    Data is series of values in an obj with keys."""

    if tick_labels is None:
        try:
            tick_labels = data.columns
        except:
            try:
                tick_labels = list(data.keys())
            except:
                raise TypeError(f"Don't know how to get tick labels from object type {type(data)}. "
                                f"Please supply tick_labels or use a different object.")
    if type(data) is dict:
        data = pd.Series(data)
    if scatter_kwargs is None:
        scatter_kwargs = {}
    if ax is None:
        ax = plt.gca()

    for xpos, col in enumerate(data):
        if per_column_kwargs is not None:
            col_kwargs = per_column_kwargs[xpos]
        else:
            col_kwargs = {}
        y = data[col]

        y = pd.Series(y).dropna()
        kde = stats.gaussian_kde(y)
        # estimate the local density at each datapoint
        density = kde(y)
        if equal_width:
            density = density / max(density) * 0.8
        # generate some random jitter between 0 and 1
        jitter = np.random.rand(len(y)) - 0.5

        # scale the jitter by the KDE estimate and add it to the centre x-coordinate
        xvals = xpos + (density * jitter * width)

        ax.scatter(xvals, y, **(scatter_kwargs|col_kwargs))

    ax.set_xticks(range(0, data.shape[1]), tick_labels)


def string_to_df(s, sep='\t'):
    return pd.DataFrame([l.split(sep) for l in s.split('\n')])