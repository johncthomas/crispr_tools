#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pathlib
import typing
from typing import Union, Any, Tuple, Dict

import pandas as pd
import numpy as np
import re
import os

from crispr_tools.data_classes import (
    analysis_short_names,
    analysis_score_names,
    get_treatment_str,
    AnalysisWorkbook
)

try:
    import gene_name_updater
except:
    print('gene_name_updater not installed')



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

#todo only one function, a high one, should require or expect the data to be stored in the
# /expid/analysis_ver/tables/ format. All other functions should be structure agnostic.

def get_statistic_table(analysis_path:str,
                        analysis_type:str, prepend_columns:str= '', outdir:str=None,
                        symbol_map:pd.Series=None, append_to=None,
                        chosen_statistics=('score', 'fdr', 'neg_p', 'pos_p')) -> typing.Dict[str, pd.DataFrame]:
    """For a specific ran analysis, get gene/comparison stats files, formated
    for use with CRISPR Screen Viewer. One table per chosen_statistics
    from a cirspr_pipeline analysis. Optionally join to existing tables.

    Args:
        analysis_path: Path to an analysis table.
        analysis_type: Type of infile, or infer from analysis_path.
        prepend_columns: ID of the experiment, must be alphanumeric ("_" is acceptable).
            Prefixed to output column headers.
        outdir: Directory to write score and FDR tables to. Return DataFrame if None.
        append_to: Directory of existing scores and lfc tables.
        symbol_map: Mapping object (dict, Series, func) giving gene names
            to be written to file.
        chosen_statistics: Names of statistics to get tables of. These should be level[1]
            column indicies of input table."""

    if type(symbol_map) is dict:
        symbol_map = pd.Series(symbol_map)

    # Check whether inputs are valid
    if not os.path.exists(analysis_path):
        raise Exception('File not found: ' + analysis_path) # Fail if infile not found

    # experiment ID can contain dash... might fail if it contains a dot?
    #todo test what is a valid expid, can it have a dot

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
        # if this changes then the lambda x: x in write_screens... needs to change.
        new_names = data.index.map(symbol_map)
        if new_names.isna().any():
            LOG.warning(f"NaN found in gene names, from original names: {data.index[new_names.isna()]}")
        data.index = new_names
        gene_dupes = data.index.duplicated()
        if gene_dupes.any():
            LOG.warning("Duplicate gene names found after mapping gene symbols: \n"
                        f"{', '.join(set(data.index[gene_dupes].dropna()))}")

    if prepend_columns:
        prepend_columns = prepend_columns + '.'

    tables = {}
    #for statcol, outstat in zip((score_name, 'fdr'), ('score', 'fdr')):
    for statcol in chosen_statistics:
        if statcol == 'score':
            outfile = f"{short_name}_score.csv"
            statcol = score_name
        else:
            outfile = f"{short_name}_{statcol}.csv"
        stats = data.xs(statcol, axis=1, level=1)

        stats.columns = [f"{prepend_columns}{x}" for x in stats.columns]

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
        symbol_mapping='symbol_mapping',
        drop_genes:typing.List[str]=['Non-targeting', 'Non-targeting control'],
) -> typing.Dict[str, pd.DataFrame]:
    """Pull all the data tables for all ran methods specified by an
    AnalysisWorkbook.

    symbol_mapping should be a path to symbol maps home"""

    LOG.info(f"Dropping {drop_genes}")

    all_results = {}

    # we are just assuming that group does not form part of the file names
    expid = analysis.expd['experiment_id']
    if analysis_version is None:
        analysis_version = analysis.expd['analysis_version']
    if results_prefix is None:
        results_prefix = analysis.expd['file_prefix']

    if symbol_mapping not in (False, None):
        lib = analysis.wb['Experiment details']['Library']

        symb_map = pd.read_csv(f'{symbol_mapping}/{lib}.csv',
                           index_col=0).Symbol
    else:
        # symb_map gets used with a df.index.map, so this just keeps the value.
        symb_map = lambda x: x

    methods = set()
    for _, row in analysis.wb['Analyses'].iterrows():
        # you can end up with blank rows in a workbook...
        if (not row['Control group']) and (not row['Method']):
            continue
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
            if drop_genes:

                tab.drop(drop_genes, inplace=True, errors='ignore')
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
        exp_id = analysis.expd['experiment_id']

        comp_row['Experiment ID'] = exp_id
        comp_row['Library'] = analysis.wb['Experiment details']['Library']
        comp_row['Comparison ID'] = f"{exp_id}.{ctrl}-{treat}"
        if pd.isna(comp_row['KO']) or (comp_row['KO'] == ''):
            comp_row['KO'] = 'WT'

        # Control group might include more information
        comp_row['Timepoint'] = group_name.split('_')[0]
        comp_row['Control group'] = group_name
        #comp_row['Source'] = source_dict[exp_id]
        comp_row['Available analyses'] = groups_methods[group_name]

        comparisons_metadata.append(comp_row)

    return pd.DataFrame(comparisons_metadata)

def gene_map_to_symbol_mapping_table(gmap:pd.Series):
    """Return a table with IDs etc, from a series mapping oldname->newname"""

    cols = ['Original_symbol', 'Symbol', 'Ensembl_gene_ID', 'NCBI_gene_ID', 'HGNC_ID']
    symbol_ids = gene_name_updater.symbol_ids_table.set_index('Approved_symbol').reindex(gmap.values)
    symbol_ids = symbol_ids.loc[:, cols[2:]]
    symbol_ids.index = gmap.index
    symbol_ids.insert(0, 'Symbol', gmap)
    symbol_ids.index.name = 'Original_symbol'
    symbol_ids = symbol_ids.loc[~symbol_ids.index.duplicated()]
    symbol_ids.columns = cols[1:]
    symbol_ids.fillna('X', inplace=True)
    return symbol_ids

def update_previously_mapped_names(symbols:pd.Series):
    """Pass series with originalSymbol->newishSymbol, newishSymbol
    will be updated, and a new symbol-IDs table generated."""
    # Test the mapped gene names against newer HGNC table
    res = gene_name_updater.update_gene_symbols(symbols.values, search_NCBI=False)
    updates = res['genes']

    for old, new in updates.loc[updates.index != updates.values].items():
        symbols[symbols == old] = new
    symbtab = gene_map_to_symbol_mapping_table(symbols)
    return symbtab


def previous_id_from_hgnc_complete(hgnc_fn, gns):
    hgnc = pd.read_csv(hgnc_fn, sep='\t', dtype=str, index_col='symbol')
    return hgnc.loc[:, ['hgnc_id', 'prev_symbol', 'entrez_id', 'ensembl_gene_id']].reindex(gns)

def write_screens_data_for_viewer(
        analysis_xlsx:ListOfPaths, results_root_dir:str, app_data_dir, append_to=None,
        results_analysis_version=None, results_file_prefix=None,
        hgnc_table=None, symbol_mapping='symbol_mapping',
        consistent_cell_lines=True,
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
    :param hgnc_table: path to hgnc_complete from https://www.genenames.org/download/archive/
        Used to include HGNC IDs and previous symbols in gene searches.
    :param symbol_mapping: if False, just use current names
    :param consistent_cell_lines: Change some cell line names for consistency.

    """
    # command line options...
    try:
        s = symbol_mapping.lower()
        if s in ('false', 'none'):
            symbol_mapping = False
    except: # maybe .lower() would crash if a Path is passed, but that's fine
        pass

    data_tables = {}
    comparisons_metadata = []
    experiment_details_list = []

    os.makedirs(app_data_dir, exist_ok=True)

    # Open xlsx on Mac create temporary files that end with .xlsx but start with ~$
    #   They are not XLSX files and will crash the pipeline.
    analysis_xlsx = [f for f in analysis_xlsx if not pathlib.Path(f).stem.startswith('~$')]

    filestr = "\n\t".join(analysis_xlsx)
    LOG.info(f'Parsing results for input experiment books: {filestr}')

    for fn in analysis_xlsx:
        LOG.info('\n\n')
        analysis = AnalysisWorkbook(fn)

        # data tables from the drugz/etc analyses
        tables = get_data_tables(
            analysis, results_root_dir, results_analysis_version,
            results_file_prefix, symbol_mapping=symbol_mapping
        )
        for k, tab in tables.items():
            try:
                data_tables[k].append(tab)
            except KeyError:
                data_tables[k] = [tab]

        comparisons_metadata.append(get_comparisons_metadata(analysis))
        experiment_details_list.append(analysis.wb['Experiment details'])

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

        # using ID in the analysis... should probably update everything to "Experiment name"
        "Experiment name":"Experiment ID",
        # Old "Citation" was actually the full reference.
        # As new citation is
        #    generated after the columns have been renamed this should be fine.
        "Citation":"Reference",
        'Analysis name':'Experiment ID',
        'Date screen completed (yyyy-mm-dd)':'Date published',
        'Date screen completed':'Date published',
        'analysis_version':'Analysis version',
        'file_prefix':'File prefix',
    }

    for expmet in experiment_details_list:
        todrop = []
        colmap = {k: k for k in expmet.index}
        for old, new in column_renamer.items():
            if old in colmap.keys():
                colmap[old] = new
                todrop.append(old)

        expmet.index = expmet.index.map(colmap)

    experiment_metadata = pd.DataFrame(experiment_details_list)
    experiment_metadata.set_index('Experiment ID', inplace=True, )

    # Short citation is generated from the full reference, automatic detection of clashing cites
    def find_date(reference):
        """Pull the date from a refernece in the MLA style."""
        pos_dates = set()
        for d in reference.split(' (')[1:]:
            d = d.split(')')[0]
            try:
                d = int(d)
                if 1900 < d < 2100:
                    pos_dates.add(d)
            except:
                pass

        pos_dates = list(pos_dates)
        if len(pos_dates) == 1:
            return pos_dates[0]
        LOG.warning(f"Multiple possible dates for reference '{reference}'. No citation created.")
        return '????'

    def short_cite_str(cite):
        year = find_date(cite)
        auth = cite.split(',')[0]
        return f"{auth} ({year})"

    # if there's no external data, and the details excel were put together with older template
    #   it's possible that we won't have a reference column in the whole table
    if "Reference" not in experiment_metadata.columns:
        experiment_metadata.loc[:, 'Reference'] = np.nan

    # should only be internal
    isinternal = lambda s: 'internal' in str(s).lower()
    noref = experiment_metadata.Reference.isna() | experiment_metadata.Reference.apply(isinternal)
    # fill with exp ID which is the index
    experiment_metadata.loc[noref, 'Reference'] = experiment_metadata.loc[noref].index
    experiment_metadata.loc[noref, 'Citation'] = experiment_metadata.loc[noref].index

    experiment_metadata.loc[~noref, 'Citation'] = experiment_metadata.loc[~noref, 'Reference'].apply(short_cite_str)

    single_tables = {}
    all_genes = set()
    for tabfn, tables in data_tables.items():
    # At the moment it's possible for their to be duplicate genes
        tab = pd.concat([t.loc[~t.index.duplicated()] for t in tables], axis=1)
        all_genes.update(tab.index)
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
                LOG.warning(f'Dupes in the index of {fn}!! {tab.index[dupes]}')
            tab.to_csv(os.path.join(app_data_dir, fn), encoding='utf-8-sig')

    if hgnc_table and app_data_dir:
        LOG.info('generating previous_and_id.csv')
        previd = previous_id_from_hgnc_complete(hgnc_table, list(all_genes))

        previd.to_csv(os.path.join(app_data_dir, 'previous_and_id.csv'))
    else:
        LOG.info('No HGNC table, no previous_and_id.csv being written')

    if consistent_cell_lines:
        cellfixes = {'HEK 293': 'HEK293', 'HEK 293A':'HEK293A', 'DLD1': 'DLD-1', 'RPE1': 'RPE-1'}
        LOG.info(f"Changing (for consistency) cell line names {', '.join(cellfixes.keys())}"
                 f" to {', '.join(cellfixes.values())}")
        for old, new in cellfixes.items():
            m = comparisons_metadata.loc[:, 'Cell line'] == old
            comparisons_metadata.loc[m, 'Cell line'] = new

    if not app_data_dir:
        return single_tables, comparisons_metadata, experiment_metadata

#todo if "File prefix" and "Analysis version" are set in details.xlsx does this work without those options.
if __name__ == "__main__":

    print("Welcome to pipeline_to_viewer.py")
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
                        help="--analysis-version of input results. Currently expects results files to be"
                             "found in '{results_dir}/{expid (from individual .xlsx)}/{analysis_version}'. ")
    parser.add_argument('--symbol-maps', default='symbol_mapping',
                        help='"NONE" or the directory containing files that map library gene names to those to'
                             'be displayed in app. If NONE the current names will be used. Default: ')
    parser.add_argument('-p', '--prefix', metavar='PREFIX',
                        help="File name prefix used when writing pipeline output files. By default"
                             " the value in the experiment xlsx is used; this over-rides that.")
    parser.add_argument('-t', '--hgnc-table',
                        help='Path to hgnc_complete .txt[.gz] (see https://www.genenames.org/download/archive/) '
                             'used to add HGNC ID and previous symbols to the gene search boxes.')

    LOG.info(f"Running from {os.getcwd()}")
    clargs = parser.parse_args()
    LOG.info(f"Parsed command line args: {clargs}")

    write_screens_data_for_viewer(
        analysis_xlsx=clargs.experiment_xlsx,
        results_root_dir=clargs.results_dir,
        app_data_dir=clargs.outdir,
        append_to=clargs.append_to,
        results_analysis_version=clargs.analysis_version,
        results_file_prefix=clargs.prefix,
        hgnc_table=clargs.hgnc_table,
    )

