#!/usr/bin/env python
"""Generate tables for use with crispr_screen_viewer.

Experiments, comparisons, per-analysis data tables,
gene mapping"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xlsxwriter

from crispr_tools.tools import tabulate_drugz, tabulate_mageck
from typing import Union, List, Dict, Any, Collection, Tuple
from crispr_tools.tools import is_ntctrl

import logging
LOG = logging.getLogger(__name__, )
logging.basicConfig()
LOG.setLevel(logging.INFO) # need to set the logger to the lowest level...

analysis_short_names = {'drugz':'drz', 'mageck':'mag'}
analysis_score_names = {'drugz':'normZ', 'mageck':'lfc'}



# statistc_colo = [scorek, 'fdr', 'neg_p', 'pos_p', 'fdr_log10']
# stat_renamed_cols = [scorek, 'FDR', 'Dropout p', 'Enrich p', 'Log10(FDR)']

ARROW = '→'

class AnalysisWorkbook:
    """Parse an Excel workbook describing how an experiment should be analysed.
    Generates the dict for JSON or calling .run_analyses

    Default method kwargs:
        mageck: --remove-zero both --remove-zero-threshold 2*pseudocount

    additional_options overwrites anything in the Excel file"""

    def __init__(self, xlfn, parse_workbook=True, counts_dir='.'):
        """
        Args:
            counts_dir: directory containing counts files specified in the
                Analyses sheet. This will be added to """

        self.safesplit = lambda x: x.replace(' ', '').split(',')
        self.fn = xlfn
        self.wb = self.load_analysis_workbook(xlfn)
        self.counts_dir = counts_dir
        if parse_workbook:
            self.expd = self.workbook_to_dict(counts_dir=counts_dir)

    def verify_wb_integrity(self, wb=None):
        if wb is None:
            wb = self.wb

        errors = []
        samp_deets = wb['Sample details']
        controls = wb['Control groups']
        analyses = wb['Analyses']

        if samp_deets.index.duplicated().any():
            errors.append('Duplicated replicate names in Sample details')

        for s in 'Control', 'Test':
            smps = controls[f'{s} sample']
            found = smps.isin(samp_deets.Sample)
            if (~found).any():
                errors.append(
                    f"Non-existent samples specified in {s} column of Control Groups sheet:\n   "
                    f"{', '.join(list(set(smps[~found])))}"
                )

        # check group names in Analyses
        missing_grps = set()
        analyses_grps = set()
        for c in analyses['Control group'].unique():
            for g in self.safesplit(c):
                analyses_grps.add(g)
                if g not in controls['Group'].unique():
                    missing_grps.add(g)

        if len(missing_grps):
            errors.append(f"Groups specified in Analyses sheet don't exist in Control groups: "
                          f"{', '.join(missing_grps)}")

        # check group names in Controls
        missing_grps = set()

        for c in controls['Group'].unique():
            if c not in analyses_grps:
                missing_grps.add(c)

        if len(missing_grps):
            errors.append(f"Groups defined in Controls sheet don't exist in Analyses: "
                          f"{', '.join(missing_grps)}")

        if errors:
            LOG.warning(f"Errors found in workbook {self.fn}.")
            LOG.warning('\n'.join(errors))
        else:
            LOG.info(f'No errors found in workbook {self.fn}.')

    def load_analysis_workbook(self, fn):
        # nothing needs to be numeric
        wb = pd.read_excel(fn, sheet_name=None, dtype=object, )
        wb = {k: tab.fillna('') for k, tab in wb.items()}

        wb['Replicate details'] = wb['Sample details'].set_index('Replicate', drop=False)
        wb['Sample details'] = wb['Replicate details'].drop_duplicates('Sample').set_index('Sample', drop=False)

        wb['Experiment details'] = pd.read_excel(
            fn, sheet_name='Experiment details',
            header=None, index_col=0,
        )[1].fillna('')

        return wb

    def workbook_to_dict(self, counts_dir=''):
        from crispr_tools.crispr_pipeline import pipeLOG, clargs_to_dict
        counts_dir = Path(counts_dir)
        pipeLOG.info(f'Parsing xlsx: {self.fn}')

        safesplit = self.safesplit

        repdeets = self.wb['Replicate details']
        ctrlgrps = self.wb['Control groups']

        # all "Experiment details" k:v are top level and don't need processingaa
        experiment_dictionary = {}

        # Some expd options will be in the workbook, store them with the appropriate
        #  internal key
        deets = self.wb['Experiment details'].to_dict()
        for col_ind, dict_key in {
            'Experiment name': 'experiment_id',
            'Analysis version': 'analysis_version',
            'Notes': 'notes',
            'File prefix': 'file_prefix',
            # some workbooks have them written in the internal style...
            'analysis_version': 'analysis_version',
            'file_prefix': 'file_prefix',
            # legacy labels...
            'Analysis name':'experiment_id',
            'Library':'library',
        }.items():
            try:
                experiment_dictionary[dict_key] = deets[col_ind]
            except KeyError:
                # Most won't exist and any other errors
                #   will get picked up later
                pass

        # control groups. produces:
        # {grpA:{ctrl_samp:[test_samp1, test_samp2, ...], ...}, ...}
        ctrl_dict = {}
        for grpn in ctrlgrps.Group.unique():
            cgrp = ctrlgrps[ctrlgrps.Group == grpn]
            # groupby gives the indexes, we want the Test sample values
            # gotta list() the reps cus json.dump can't handle an ndarray
            g = {c: list(cgrp.loc[i, 'Test sample'].values)
                 for c, i in cgrp.groupby('Control sample').groups.items()}
            ctrl_dict[grpn] = g
        experiment_dictionary['control_groups'] = ctrl_dict

        # sample_reps, just {sampA:[repA1, repA2], ...}
        experiment_dictionary['sample_reps'] = {
            k: list(v.values) for k, v in repdeets.groupby('Sample').groups.items()
        }

        # Analyses. Producing a list of dicts with "method" and "groups":List
        #   required keys. Should practically always be a fn_counts, but ignore
        #   if blank. kwargs & pseudocount added if not present
        analyses = []
        # pipeLOG.info('workbook_to_dict: Processing analyses')
        for _, row in self.wb['Analyses'].iterrows():
            pipeLOG.debug(str(row))

            for method in safesplit(row['Method']):
                andict = {'method': method,
                          'groups': safesplit(row['Control group'])}

                if 'kwargs' not in andict:
                    andict['kwargs'] = {}

                if row['Add pseudocount']:
                    psdc = int(row['Add pseudocount'])
                else:
                    psdc = 1
                andict['pseudocount'] = psdc

                if row['Counts file']:
                    andict['counts_file'] = os.path.join(
                        counts_dir,
                        row['Counts file']
                    )

                # add default args if any
                if not row['Arguments']:
                    if method == 'mageck':
                        pc = andict['pseudocount']
                        if pc > 1:
                            cutoff = 2 * pc
                        else:
                            cutoff = 0
                        andict['kwargs'] = {'remove-zero': 'both',
                                            'remove-zero-threshold': cutoff}
                else:
                    rowargs = row['Arguments']
                    try:

                        kwargs = eval(rowargs)
                        if not type(kwargs) is dict:
                            raise RuntimeError()
                    except:
                        pipeLOG.info(f"Parsing analysis Arguments, {rowargs}, as command line style")
                        # try parsing as command line options
                        kwargs = clargs_to_dict(rowargs)
                    andict['kwargs'] = kwargs
                    pipeLOG.info(f"Arguments: {rowargs} parsed to {kwargs}")

                # deal with paired, which is handled different in the programs
                if not pd.isna(row['Paired']):
                    if row['Paired'] and (method == 'mageck'):
                        andict['kwargs']['paired'] = ''
                    elif not row['Paired'] and (method == 'drugz'):
                        andict['kwargs']['unpaired'] = True

                # Name is an optional col
                if 'Name' in row and row['Name']:
                    andict['name'] = row['Name']
                else:
                    andict['name'] = '-'

                analyses.append(andict)
        experiment_dictionary['analyses'] = analyses

        # default filename prefix
        if 'file_prefix' not in experiment_dictionary:
            experiment_dictionary['file_prefix'] = experiment_dictionary['experiment_id']

        return experiment_dictionary

    def iter_comps(self, join_comps: str = False) -> Tuple[str, Union[str, Tuple]]:
        """Pass a string to join_comps to join the treatment and test sample
        names, otherwise a tuple of the sample names is returned."""
        for _, row in self.wb['Control groups'].iterrows():
            comp = row['Control sample'], row['Test sample']
            if join_comps:
                comp = join_comps.join(comp)
            yield (row.Group, comp)

    def check_analyses_ran(self, results_dir, analysis_version=None,
                           file_prefix=None, include_group_names=False):
        # tests for this should include groups that only run in some methods as well as missing comps
        self.safesplit = lambda x: x.replace(' ', '').split(',')
        safesplit = self.safesplit

        xd = self.wb['Experiment details']
        if file_prefix is None:
            file_prefix = xd['File prefix']
        if analysis_version is None:
            analysis_version = xd['Analysis version']

        full_dir = os.path.join(results_dir, xd['Experiment name'], analysis_version, )

        missing_comps = []
        for _, row in self.wb['Analyses'].iterrows():
            for grp in safesplit(row['Control group']):
                if include_group_names:
                    table_str = os.path.join(full_dir, 'tables', f"{grp}.{file_prefix}")
                else:
                    table_str = os.path.join(full_dir, 'tables', file_prefix)
                for meth in safesplit(row['Method']):
                    tab = pd.read_csv(table_str + f'.{meth}_table.csv', header=[0, 1], index_col=0)

                    for g, comp in self.iter_comps('-'):
                        if g == grp:
                            if comp not in tab.columns.levels[0]:
                                missing_comps.append(comp)
        if missing_comps:
            print(f'There are {len(missing_comps)} missing comparisons in {xd["Experiment name"]}')
        else:
            print(f'{xd["Experiment name"]} appears to have all comparisons written.')


def write_analysis_output_excel(
        results_table:Union[pd.DataFrame, str],
        control_groups:pd.DataFrame,
        out_prefix,
        analysis_type,
        one_file_per_ctrlgrp=True):
    """Write excel workbooks for each control group in an analysis"""

    if type(results_table) != pd.DataFrame:
        results_table = pd.read_csv(results_table, header=[0,1], index_col=0)

    def write_number_sheet(sheet, df):
        """Writes columns and row names formated as string,
        and numbers formated with 3 dp"""
        # write the index header
        sheet.write(0, 0, 'Gene')

        # write the stats headers
        for coli, c in enumerate(df.columns):
            # print(c)
            sheet.write(0, coli + 1, c, header_format)

        # row labels
        for row_i, gene in enumerate(df.index):
            sheet.write(1 + row_i, 0, gene, index_format)

        # data
        for col_i, col in enumerate(df):
            for row_i, val in enumerate(df[col]):
                try:
                    sheet.write(row_i + 1, col_i + 1, val, num_format)
                # catch nans
                except TypeError:
                    if pd.isna(val):
                        pass
                    else:
                        raise

    def color_scale_score(sheet, first_row, first_col, last_row, last_col):
        sheet.conditional_format(
            first_row, first_col, last_row, last_col,
            {
                'type': '3_color_scale',
                'min_type': 'num',
                'mid_type': 'num',
                'max_type': 'num',
                'min_value': -3,
                'mid_value': 0,
                'max_value': 3,
                'min_color': '#ff5030',
                'mid_color': '#ffffff',
                'max_color': '#6699ff'
            }
        )

    def color_scale_fdr(sheet, first_row, first_col, last_row, last_col):
        sheet.conditional_format(
            first_row, first_col, last_row, last_col,
            {
                'type': '2_color_scale',
                'min_type': 'num',
                'max_type': 'num',
                'min_value': 0.001,
                'max_value': 0.25,
                'min_color': '#3cdd3c',
                'max_color': '#ffffff'
            }
        )

    scorek = analysis_score_names[analysis_type] ###
    stat_cols = [scorek, 'fdr', 'neg_p', 'pos_p', 'fdr_log10']
    stat_renamed_cols = [scorek, 'FDR', 'Dropout p', 'Enrich p', 'Log10(FDR)']
    control_groups.loc[:, 'Timepoint'] = control_groups['Group'].map(lambda x: x.split('_')[0])

    if one_file_per_ctrlgrp:
        tpgrps = control_groups.groupby('Timepoint').groups
    else:
        tpgrps = {'all_results':control_groups.index}

    for grp, indx in tpgrps.items():
        filename = f'{out_prefix}.{grp}.xlsx'

        workbook = xlsxwriter.Workbook(filename, )
        index_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'right'})
        header_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'center'})
        num_format = workbook.add_format({'num_format': '0.000'})

        grpdict = {}
        T = control_groups.loc[indx]
        comps = T['Control sample'] + '-' + T['Test sample']

        # Write sheet of just scores
        score = results_table.xs(scorek, 1, 1).loc[:, comps]
        scoresheet = workbook.add_worksheet(name=scorek)
        write_number_sheet(
            scoresheet,
            score
        )
        color_scale_score(
            scoresheet,
            1, 1,
            score.shape[0]+1,
            score.shape[1]+1
        )

        # write all the per comparison sheets
        for comp in comps:
            grpdict[comp] = res = results_table.loc[:, (comp, stat_cols)]
            res.columns = stat_renamed_cols
            res.sort_values(scorek, inplace=True)
            sheet = workbook.add_worksheet(name=comp.replace('-', ARROW))

            write_number_sheet(sheet, res)
            last_row = res.shape[0] + 1

            color_scale_score(sheet, 1, 1, last_row, 1)
            color_scale_fdr(sheet, 1, 2, last_row, 2)

        workbook.close()


class CrisprCounts:
    """For holding and working with counts data from CRISPR screen or whatever.

    Args:
        fn_or_df: path to counts file or counts dataFrame.
        guide_hdr: Column header for guides
        gene_hdr: Column header for genes

    Props:
        counts: df holding numerical counts data
        guide_gene: mapping of counts.index to gene names
        gene_guides: mapping of gene name to guides (one to many)

    Methods (see method doc strings):
        apply_to_genes
        select_counts_by_genes
        counts_with_genes
        write_counts

    Notes: Should deal with dropping of rows; guide_gene etc are properties
    regenned using self.counts.index self._orig_guide_name every time
    they're referenced.
    Won't deal with addition of rows, create new instance with
    concatted count DFs or overwrite self._orig_guide_gene
    """
    def __init__(self, fn_or_df,
                 guide_hdr='guide', gene_hdr='gene',
                 sep='\t'):

        # prep, validate the dataframe
        # prep, validate the dataframe
        self.guidek = guide_hdr
        self.genek = gene_hdr
        if type(fn_or_df) != type(pd.DataFrame({})):
            df = pd.read_csv(fn_or_df, index_col=guide_hdr, sep=sep)
            if df.shape[1] == 0:
                raise RuntimeWarning(
                    "No counts in the loaded DF, check this_obj.counts, "
                    "maybe the wrong separator was used?"
                )
        else:
            df = fn_or_df
        # if df passed without indexing by guides, fix that.
        if self.guidek in df.columns:
            df.set_index(self.guidek, inplace=True)
        assert gene_hdr in df.columns

        g = df[self.genek]
        self._orig_guide_gene = g.copy()
        self._genegroups = g.groupby(g.values).groups
        df.drop(self.genek, axis=1, inplace=True)
        self.counts = df

        # Used for checking alterations to the df when
        #   selecting by genes
        self._indexid = id(df.index)
        # don't think we'll need this, but just in case
        self._columnid = id(df.columns)

    @property
    def guide_gene(self):
        curr_guides = self.counts.index
        return self._orig_guide_gene.reindex(curr_guides)

    @property
    def genes(self):
        return self.guide_gene.unique()

    @property
    def gene_guides(self):
        """Mapping of gene to guide ID"""
        # check for alterations to the index, recalc groups if required
        idx_id = self.counts.index
        if id(idx_id) != self._indexid:
            g = self.guide_gene
            grps = g.groupby(g.values).groups
            self._genegroups = grps
            self._indexid = idx_id

        return self._genegroups

    @property
    def apply_to_genes(self):
        return self.counts.groupby(self.guide_gene).apply

    def select_counts_by_genes(self, genes:Union[str, Collection[str]]):
        """Return rows from self.counts for specific gene(s). Pass a single name,
         or iterable of names"""
        if type(genes) == str:
            inds = self.gene_guides[genes]
        else:
            inds = []
            for g in genes:
                inds.extend(self.gene_guides[g])
        return self.counts.loc[inds]

    def get_counts_with_genes(self):
        """Get copy of counts dataframe including gene column"""
        cnt = self.counts.copy()
        cnt.insert(0, self.genek, self.guide_gene)
        return cnt

    def write_counts(self, filename, sep='\t'):
        """self.guide_gene is reinserted into self.counts and TSV file
        is written."""
        cnt = self.get_counts_with_genes()
        cnt.to_csv(filename, sep=sep)

    def rename_genes(self, oldname_newname:dict):
        og = self._orig_guide_gene
        for old, new in oldname_newname.items():
            og[og == old] = new
        self._orig_guide_gene = og


    # todo plotting methods


def load_symbol_map(fn) -> Dict[str,str]:
    """Return dict of library symbol to official symbol. Handle both old style .txt
    and more recent .csv"""
    fn = str(fn)
    smap = {}
    if fn.endswith('.txt'):
        with open(fn) as f:
            for line in f:
                a, b = line.strip().split('\t')
                smap[a] = b
                return smap
    elif fn.endswith('.csv'):
        tab = pd.read_csv(fn, index_col=0)
        return tab['Symbol'].to_dict()
    elif fn.endswith('.tsv'):
        tab = pd.read_csv(fn, sep='\t', index_col=0)
        return tab['Symbol'].to_dict()
    else:
        raise RuntimeError(f'Symbol map file must be headerless .txt; or .csv, or .tsv, file name provided: {fn}')



def load_yaml(paff):
    import yaml
    with open(paff) as f:
        obj = yaml.safe_load(f)
    return obj


def load_json(paff):
    """Load JSON from file path, ignoring lines that begin with #."""
    lines = []
    with open(paff) as f:
        for line in f:
            if line[0] != '#':
                lines.append(line)

        obj = json.loads('\n'.join(lines))
    return obj


def add_pseudocount_and_save(fn, pseudocount=4):
    fn = str(fn)
    alldub = pd.read_csv(fn, '\t', index_col=0)
    alldub.iloc[:, 1:] += pseudocount
    out_fn = fn.replace('.tsv', f'.psd{pseudocount}.tsv')
    print(out_fn)
    alldub.to_csv(out_fn, '\t')


# def map_and_ambig(gset:np.ndarray) -> (pd.Series, np.ndarray, np.ndarray):
    # now is gene_name_updater.update_gene_symbols


def print_fasta_style(genes, library, n_guides=3):
    """Useful for outputing sequences to BLAT or BLAST"""
    for g in genes:
        gds = library.loc[library.gene == g, 'seq'].values[:n_guides]
        for i, gd in enumerate(gds):
            print(f">{g}_{i}")
            print(gd)

# class Metadata:
#     """Class for holding the 3 metadata objects: analysis options, sample
#     details and experiment details.
#
#     Assumes that within source_dir there are expjson and
#     details directories with relevant files called {expid}.json/xlsx.
#
#     Use either source_xlsx to specify a single XLSX with all relevant info
#     or
#
#     source_dir expects that {source_dir}/expjson/{expid}.json
#     and {source_dir}/details/{expid}.xlsx both exist."""
#     def __init__(self, expid='', source_dir='', source_xlsx=''):
#         if not source_dir and not source_xlsx:
#             raise RuntimeError('One of source_dir or source_xlsx must be specified')
#
#         # parse metadata from multiple files in a directory
#         if source_dir:
#             d = {}
#             d['opt']  = load_json(
#                 os.path.join(
#                     source_dir,
#                     "expjson",
#                     f"{expid}.json"
#                 )
#             )
#
#             fn = os.path.join(
#                 source_dir,
#                 'details',
#                 f'{expid}.xlsx'
#             )
#
#             d['exp'] = pd.read_excel(fn, index_col=0, header=None).iloc[:, 0]
#             d['smp'] = pd.read_excel(fn, sheet_name='Sample details').drop_duplicates('Sample').set_index('Sample')
#         # parse metadata from single Excel workbook.
#         else:
#             d = {}
#             wb = AnalysisWorkbook(source_xlsx)
#             d['opt'] = wb.expd
#             d['smp']  = wb.wb['Sample details'].drop_duplicates('Sample').set_index('Sample')
#             d['exp']  = wb.wb['Experiment details']
#         self.d = d
#         self.opt = d['opt']
#         self.options = d['opt']
#         self.smp = d['smp']
#         self.sample_details = d['smp']
#         self.exp = d['exp']
#         self.experiment_details = d['exp']
#         self.exp_name = self.exp['Experiment name']
#         self.exp_id = self.exp_name

def is_nt(s):
    return pd.isna(s) or (s == 'None') or (s == 'DMSO') or (s=='') or (s=='WT')

def get_treatment_str(samp_deets:pd.DataFrame, ctrl:str, treat:str):
    """Return a string describing the treatment performed between the
    control and test samples.

    Args:
        samp_deets: dataframe, rows indexed by samp IDs, giving sample
            metadata.
        ctrl, treat: Sample IDs of a comparison.


    Divides possible treatments into chemical or KO;
    outputs a string for a chemical treatment, a gene-KO
    treatment, a chemical treatment in a KO background,
    or a KO in a chemical background. Depending o """

    fatarrow = '➤'

    t_deets, c_deets = samp_deets.loc[treat], samp_deets.loc[ctrl]

    # get chem treat str
    if not is_nt(t_deets.Treatment):
        if is_nt(c_deets.Treatment) or (c_deets.Treatment == t_deets.Treatment):
            chem_str = t_deets.Treatment
        # if both ctrl/treat samps are different treatments
        else:
            chem_str = f"{c_deets.Treatment}{fatarrow}{t_deets.Treatment}"
    else:
        chem_str = ''

    # Get the KO str.
    #  We assume that we'll never have a KO control and WT treat samp
    if is_nt(t_deets.KO):
        ko_str = ''
    else:
        if is_nt(c_deets.KO) or (c_deets.KO == t_deets.KO):
            ko_str = t_deets.KO+'-KO'
            # if it's just in a TP53 background we don't care that much
            if ko_str == 'TP53-KO':
                ko_str = ''
        else:
            ko_str = f"{c_deets.KO}-KO{fatarrow}{t_deets.KO}-KO"

    if chem_str and not ko_str:
        treatment_str = chem_str
    elif ko_str and not chem_str:
        treatment_str = ko_str
    elif not ko_str and not chem_str:
        treatment_str = 'No treatment'
    else:
        if t_deets.Treatment == c_deets.Treatment:
            # it's a KO in a chemical background
            treatment_str = f"{ko_str} (with {chem_str})"
        else:
            treatment_str = f"{chem_str} (in {ko_str} cells)"

    return treatment_str


# # pretty sure this is all replaced by code in pipeline_to_viewer
# def parse_analysis_metadata(analysis_workbook:str, suppress_expid=False):
#     """Write comparisons_metadata.csv and exeriments_metadata.csv.
#
#     Args:
#         suppress_expid: Don't prefix expid to comparison IDs"""
#
#     # parse out the tables from the workbook
#     deets = pd.read_excel(analysis_workbook, sheet_name=None)
#
#     rep_deets = deets['Replicate details']
#     samp_deets = rep_deets.drop_duplicates('Sample').set_index('Sample')
#
#     # Experiment details parsed to a dictionary
#     exp_deets = deets['Experiment details']
#     exp_deets = exp_deets.set_index(exp_deets.columns[0]).iloc[:, 0].to_dict()
#     expid = deets['Experiment details'].columns[-1]
#     exp_deets['ExpID'] = expid
#     exp_deets['Experiment name'] = expid
#     comp_table = deets['Control groups']
#
#     analyses = deets['Analyses']
#
#     # Get the analysis methods used for each group (which will
#     #   be output as 'Available analyses' for each comp)
#     safesplit = lambda x: x.replace(' ', '').split(',')
#     grp_method = {}
#     for _, row in analyses.iterrows():
#         for grp in safesplit(row['Control group']):
#             for meth in safesplit(row['Method']):
#                 try:
#                     grp_method[grp].append(meth)
#                 except KeyError:
#                     grp_method[grp] = [meth]
#
#     # construct the comparisons table; first as List[Dict]
#     comparisons_metadata = []
#     for _, row in comp_table.iterrows():
#         ctrl = row['Control sample']
#         test = row['Test sample']
#         ctrl_group = row['Group']
#
#         treat_row = samp_deets.loc[test].fillna('')
#
#         compid = f"{ctrl}-{test}"
#         if not suppress_expid:
#             compid = f"{expid}.{compid}"
#
#         out_row = {
#             'Comparison ID':compid,
#             'Experiment ID':expid,
#             # timepoint should just be 'endpoints', 'otherprior', or 'fromstart'.
#             'Timepoint':ctrl_group.split('_')[0],
#             'Treatment':get_treatment_str(samp_deets, ctrl, test),
#         }
#         # these values copied verbatim
#         for k in ('Dose', 'Growth inhibition %', 'Days grown', 'Cell line', 'KO'):
#             out_row[k] = treat_row[k]
#
#         out_row['Library'] = exp_deets['Library']
#
#         out_row['Available analyses'] = '|'.join(grp_method[ctrl_group])
#         out_row['Control group'] = ctrl_group
#         comparisons_metadata.append(out_row)
#     comparisons_metadata = pd.DataFrame(comparisons_metadata, ).set_index('Comparison ID')
#
#     comparisons_metadata.loc[comparisons_metadata.KO == '',  'KO'] = 'WT'
#
#     # Construct the experiments table, this is just copied from the
#     #   experiments sheet.
#     exp_deets['Treatments'] = '|'.join(samp_deets.Treatment.dropna().unique())
#     experiments_metadata = pd.DataFrame([exp_deets]).set_index('ExpID')
#
#     return experiments_metadata, comparisons_metadata
#
# def consolidate_metadata(workbooks:List[str], outdir=None):
#     """Consolidate metadata contained in provided list of analysis workbooks."""
#
#     # Consolidate
#     exp_metadatas = []
#     comp_metadatas = []
#     for fn in workbooks:
#         exp, comp = parse_analysis_metadata(fn)
#         exp_metadatas.append(exp)
#         comp_metadatas.append(comp)
#     comp_metadatas = pd.concat(comp_metadatas, axis=0)
#     exp_metadatas = pd.concat(exp_metadatas, axis=0)
#
#     if outdir:
#         comp_metadatas.to_csv(
#             os.path.join(outdir, 'comparisons_metadata.csv'),
#             encoding='utf-8-sig'
#         )
#         exp_metadatas.to_csv(
#             os.path.join(outdir, 'experiments_metadata.csv'),
#             encoding='utf-8-sig'
#         )
#
#     return exp_metadatas, comp_metadatas




if __name__ == '__main__':
    pass


