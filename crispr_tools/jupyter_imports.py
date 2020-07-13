"""Imports things commonly used by the author."""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
import yaml
from functools import partial
revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])
import logging

pltlogger = logging.getLogger('matplotlib')
pltlogger.setLevel(logging.WARNING)

import os
from typing import Dict, List, Union, Tuple

from crispr_tools.qc import get_clonal_lfcs
from crispr_tools import *
from IPython.display import display
from jacks.jacks_io import runJACKS
import jacks
jacks.jacks_io.LOG.setLevel(logging.WARNING)

from .tools import hart_list

from sklearn.metrics import precision_recall_curve

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import itertools

from itertools import combinations, combinations_with_replacement
import statsmodels.api as sm

OLS = sm.regression.linear_model.OLS

def multipletests_fdr(ps, method='fdr_bh', **kwargs):
    """Calls statsmodels.stats.multipletests and returns the corrected
    p-values only, with default of fdr_bh, rather than a FWER method."""
    return sm.stats.multipletests(ps, method=method, **kwargs)[1]


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
            different columns per sheet.
            Note: if you can't specify the index, so copy that column to the DF.
            you can also set index=False to avoid the duplicate, but that doesn't
            work with multi-index columns.
        **to_excel_kwargs: Additional kwargs are passed to pd.DataFrame.to_excel()
    """
    # tests 'all', text_col is list or dict

    writer = pd.ExcelWriter(filename,
                            engine='xlsxwriter')

    # Go through the DF
    for sheet_name, df in sheets.items():
        df.to_excel(writer, sheet_name=sheet_name, **to_excel_kwargs)

        worksheet = writer.sheets[sheet_name]

        workbook = writer.book
        txt = workbook.add_format({'num_format': '@'})

        # get the columns to be textualised
        if text_columns:
            if text_columns == 'ALL':
                txtcols = df.columns
            # if we aren't specifying per sheet
            elif type(text_columns) is list:
                txtcols = text_columns
            else:
                txtcols = text_columns[sheet_name]
                if txtcols == 'ALL':
                    txtcols = df.columns

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


