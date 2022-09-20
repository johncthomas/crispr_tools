#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 17:00:24 2022

@author: simonlam, John C. Thomas
contact: sl681@cam.ac.uk, jcthomas000@gmail.com
"""

import argparse
import typing
from typing import Union, Any, Tuple, Dict

import pandas as pd
import re
import os

from crispr_tools.dataset import (
    analysis_short_names,
    analysis_score_names,
    get_treatment_str,
    AnalysisWorkbook
)

import logging
logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel('INFO')

ListOfPaths = typing.List[typing.Union[str, os.PathLike]]

"""
Take an output DrugZ or MAGeCK table created using crispr_pipeline.py,
extract relevant columns for plotting scores and fdrs on DDRcs, and join
to existing tables with scores and fdrs.

"""
__version__ = "1.1"
"""
Version history:
    1.0 - Initial release
    1.1 - Integrate into crispr_tools
"""

def get_statistic_table(analysis_path:str,
                        analysis_type:str, experiment_id:str, outdir:str=None, append_to:str=None,
                        symbol_map:pd.Series=None) -> typing.Dict[str, pd.DataFrame]:
    """Write gene/comparison stats files for use with CRISPR Screen Viewer,
    from a cirspr_pipeline analysis. Optionally join to existing tables.

    Args:
        analysis_path: Path to an analysis table.
        analysis_type: Type of infile, or infer from analysis_path.
        experiment_id: ID of the experiment, must be alphanumeric ("_" is acceptable).
            Prefixed to output column headers.
        outdir: Directory to write score and FDR tables to. Return DataFrame if None.
        append_to: Directory of existing scores and lfc tables.
        symbol_map: Mapping object (dict, Series, func) giving gene names
            to be written to file."""

    if type(symbol_map) is dict:
        symbol_map = pd.Series(symbol_map)

    # Check whether inputs are valid
    if not os.path.exists(analysis_path):
        raise Exception('File not found: ' + analysis_path) # Fail if infile not found

    # experiment ID can contain dash... might fail if it contains a dot?
    #todo test what is a valid expid, can it have a dot

    # if experiment_id:
    #     if re.search("[^A-Za-z0-9_]", experiment_id) is not None:
    #         raise Exception('Invalid --prefix in ' + experiment_id + '. Must be alphanumeric + "_".') # Fail if invalid prefix

    analysis_type = analysis_type.lower() # Lowercase kind
    if analysis_type == "infer": # Infer kind from infile name
        for analysis_name in analysis_short_names.keys():
            if analysis_name in analysis_path:
                analysis_type = analysis_name
    if analysis_type not in analysis_short_names.keys():
        # Fail if --kind is not one of the options
        raise Exception(f'Unexpected value for --kind. Got "{analysis_type}", expected one of: '
                        f'{", ".join(analysis_short_names.keys())}, or "infer".')

    short_name = analysis_short_names[analysis_type]
    score_name = analysis_score_names[analysis_type]

    if append_to is not None: # Check if files exist if append_to is specified
        append_score = os.path.join(append_to, f"/{short_name}_score.csv")
        append_fdr = os.path.join(append_to, f"/{short_name}_fdr.csv")

        if not os.path.exists(append_score):
            raise FileNotFoundError('File not found: ' + append_score + '.')
        elif not os.path.exists(append_fdr):
            raise FileNotFoundError('File not found: ' + append_fdr + '.')

    data = pd.read_csv(analysis_path, header=[0, 1], index_col=0)
    if symbol_map is not None:
        new_names = data.index.map(symbol_map)
        if new_names.isna().any():
            LOG.warning(f"NaN found in gene names, from original names: {data.index[new_names.isna()]}")
        data.index = new_names
        gene_dupes = data.index.duplicated()
        if gene_dupes.any():
            LOG.warning("Duplicate gene names found after mapping gene symbols: \n"
                        f"{', '.join(set(data.index[gene_dupes].dropna()))}")

    if experiment_id:
        experiment_id = experiment_id + '.'

    tables = {}
    for statcol, outstat in zip((score_name, 'fdr'), ('score', 'fdr')):
        stats = data.xs(statcol, axis=1, level=1)

        stats.columns = [f"{experiment_id}{x}" for x in stats.columns]
        outfile = f"{short_name}_{outstat}.csv"

        tables[outfile] = stats

        if outdir is not None:
            if not os.path.exists(outdir):
                os.makedirs(os.path.dirname(outdir), exist_ok=True) # Make the outdir if it doesn't exist already
            if append_to is not None:
                append_data = pd.read_csv(os.path.join(append_to, outfile), index_col=0)
                stats = pd.concat([append_data, stats], axis=1)
            outpath = os.path.join(outdir, outfile)
            stats.to_csv(outpath)

            LOG.info("Successfully wrote to file: " + outpath)
    return tables



def get_data_tables(
        analysis:AnalysisWorkbook, results_root_dir:str,
        analysis_version=None, results_prefix=None,
) -> typing.Dict[str, pd.DataFrame]:
    """Pull all the data tables for all ran methods specified by an
    AnalysisWorkbook."""

    all_results = {}

    # we are just assuming that group does not form part of the file names
    expid = analysis.expd['experiment_id']
    if analysis_version is None:
        analysis_version = analysis.expd['analysis_version']
    if results_prefix is None:
        results_prefix = analysis.expd['file_prefix']

    lib = analysis.wb['Experiment details']['Library']
    symb_map = pd.read_csv(f'symbol_mapping/{lib}.csv',
                           index_col=0).Symbol

    methods = set()
    for _, row in analysis.wb['Analyses'].iterrows():
        for meth in analysis.safesplit(row['Method']):
            methods.add(meth)
    for meth in methods:
        res = get_statistic_table(
            os.path.join(
                results_root_dir,
                expid,
                analysis_version,
                'tables',
                f"{results_prefix}.{meth}_table.csv"
            ),
            meth,
            expid,
            symbol_map=symb_map,
        )

        for k, tab in res.items():
            all_results[k] = tab

    return all_results


def get_comparisons_metadata(analysis:AnalysisWorkbook, ) -> pd.DataFrame:
    """Tabulate the comparisons from an analysis Excel workbook."""
    safesplit = analysis.safesplit

    # get methods used for each control group
    groups_methods = {}
    for _, row in analysis.wb['Analyses'].iterrows():
        for grp in safesplit(row['Control group']):
            for meth in safesplit(row['Method']):
                try:
                    groups_methods[grp].append(meth)
                except KeyError:
                    groups_methods[grp] = [meth]
    groups_methods = {k:'|'.join(m) for k, m in groups_methods.items()}

    # this will become the returned table
    comparisons_metadata = []

    # iterate through comparisons, pull the sample data for controls/treatments and build the table
    for group_name, (ctrl, treat) in analysis.iter_comps():
        comp_row = {}
        #print(ctrl, treat)
        comp_row['Ctrl samp'] = ctrl
        comp_row['Treat samp'] = treat
        comp_row['Treatment'] = get_treatment_str(
            analysis.wb['Sample details'],
            ctrl, treat
        )

        for k in ('Dose', 'Growth inhibition %', 'Days grown',
                  'Cell line',  'KO', 'Notes'):
            comp_row[k] = analysis.wb['Sample details'].loc[treat, k]
        exp = analysis.wb['Experiment details']
        comp_row['Experiment ID'] = exp_id = exp['Experiment name']
        comp_row['Library'] = exp['Library']
        comp_row['Comparison ID'] = f"{exp_id}.{ctrl}-{treat}"
        if pd.isna(comp_row['KO']) or (comp_row['KO'] == ''):
            comp_row['KO'] = 'WT'
        comp_row['Timepoint'] = group_name.split('_')[0]
        comp_row['Control group'] = group_name
        #comp_row['Source'] = source_dict[exp_id]
        comp_row['Available analyses'] = groups_methods[group_name]

        comparisons_metadata.append(comp_row)

    return pd.DataFrame(comparisons_metadata)



def write_screens_data_for_viewer(
        analysis_xlsx:ListOfPaths, results_root_dir:str, app_data_dir, append_to=None,
        results_analysis_version=None, results_file_prefix=None,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    """Convert the outputs of crispr_pipeline.py into tables of data formated
    for crispr_screen_viewer package. Write tables.

    :param analysis_xlsx: List of paths to experiment XLSX files.
    :param results_root_dir: Source directory for results.
    :param app_data_dir: Output directory.
    :param append_to: Data from here will be appended to the end of generated tables.
    :param results_analysis_version: --analysis-version used when writing pipeline output files.
        By default, the value in the experiment xlsx is used; this over-rides that. Optional
    :param results_file_prefix: --file-prefix used when writing pipeline output files.
        By default, the value in the experiment xlsx is used; this over-rides that.Optional

    """
    data_tables = {}
    comparisons_metadata = []
    experiment_metadata = []

    os.makedirs(app_data_dir, exist_ok=True)
    filestr = "\n\t".join(analysis_xlsx)
    LOG.info(f'Parsing results for input experiment books: {filestr}')

    for fn in analysis_xlsx:
        LOG.info('\n\n')
        analysis = AnalysisWorkbook(fn)

        # data tables from the drugz/etc analyses
        tables = get_data_tables(
            analysis, results_root_dir, results_analysis_version,
            results_file_prefix
        )
        for k, tab in tables.items():
            try:
                data_tables[k].append(tab)
            except KeyError:
                data_tables[k] = [tab]

        comparisons_metadata.append(get_comparisons_metadata(analysis))
        experiment_metadata.append(analysis.wb['Experiment details'])

    comparisons_metadata = pd.concat(comparisons_metadata).set_index('Comparison ID')
    compdupes = comparisons_metadata.index.duplicated()
    if compdupes.any():
        print(f'duplicate comparisons found, they will be removed: {comparisons_metadata.index[compdupes]}')
        comparisons_metadata = comparisons_metadata.loc[~compdupes]

    # The experiments metadata should just be taken verbatim from details sheet
    # but some of the column names are verbose or have changed. And may change
    # in the future.
    column_renamer = {
        'Experiment description (a few sentances)':'Experiment description',
        'Experiment description (a few sentences)':'Experiment description',
        'Experiment ID':'Experiment name',
    }
    experiment_metadata = pd.DataFrame(experiment_metadata)
    colmap = {k:k for k in experiment_metadata.columns}
    for old, new in column_renamer.items():
        if old in colmap.keys():
            colmap[old] = new
    experiment_metadata.columns = experiment_metadata.columns.map(colmap)
    experiment_metadata.set_index('Experiment name', inplace=True)

    single_tables = {}
    for tabfn, tables in data_tables.items():
    # At the moment it's possible for their to be duplicate genes (in the
    #   TKO libraries specifically)
        tab = pd.concat([t.loc[~t.index.duplicated()] for t in tables], axis=1)
        single_tables[tabfn] = tab

    metatabs = {'comparisons_metadata.csv':comparisons_metadata,
     'experiments_metadata.csv':experiment_metadata}
    for axis, tabledict in [(1, single_tables), (0, metatabs)]:
        for fn, tab in tabledict.items():
            if append_to is not None:
                #todo this should deal with duplicates... maybe a validate new dataset function.
                # discard older data
                append_data = pd.read_csv(os.path.join(append_to, fn), index_col=0, )
                tab = pd.concat([append_data, tab], axis=axis)
            dupes = tab.index.duplicated()
            if dupes.any():
                LOG.warning()
            tab.to_csv(os.path.join(app_data_dir, fn), encoding='utf-8-sig')
    if not app_data_dir:
        return single_tables, comparisons_metadata, experiment_metadata


if __name__ == "__main__":

    print("Welcome to pipeline_to_viewer.py, version " + __version__)
    # Take commandline arguments
    parser = argparse.ArgumentParser(
        description="Convert the outputs of crispr_pipeline.py into tables of data formated for "
                    "crispr_screen_viewer package.")
    parser.add_argument('experiment_xlsx', nargs='+',
                        help="Paths to experiment .xlsx files")
    parser.add_argument('-o', '--outdir', metavar='PATH', required=True,
                        help="Directory to which all tables will be written. Created if not exant.")
    parser.add_argument('-r', '--results-dir', required=True,
                        help='Path containing the pipeline results used as input.')
    parser.add_argument('-a', '--append-to', metavar='PATH',
                        help="Directory containing existing screen viewer data files.",)
    parser.add_argument('-v', '--analysis-version', metavar='PREFIX',
                        help="--analysis-version used when writing pipeline output files. By default"
                             " the value in the experiment xlsx is used; this over-rides that.")
    parser.add_argument('-p', '--prefix', metavar='PREFIX',
                        help="File name prefix used when writing pipeline output files. By default"
                             " the value in the experiment xlsx is used; this over-rides that.")

    clargs = parser.parse_args()

    write_screens_data_for_viewer(
        analysis_xlsx=clargs.experiment_xlsx,
        results_root_dir=clargs.results_dir,
        app_data_dir=clargs.outdir,
        append_to=clargs.append_to,
        results_analysis_version=clargs.analysis_version,
        results_file_prefix=clargs.prefix,
    )