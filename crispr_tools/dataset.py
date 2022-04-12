import json
import os

import numpy as np
import pandas as pd

from crispr_tools.crispr_pipeline import AnalysisWorkbook
from crispr_tools.tools import tabulate_drugz, tabulate_mageck

"""Generate tables for use with crispr_screen_viewer.

Experiments, comparisons, per-analysis data tables,
gene mapping"""

from typing import Union, List, Dict, Any, Collection
from crispr_tools.tools import is_ntctrl

analysis_short_names = {'drugz':'drz', 'mageck':'mag'}

#todo no timepoint in comparisons??

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

class Metadata:
    """Class for holding the 3 metadata objects: analysis options, sample
    details and experiment details.

    Assumes that within source_dir there are expjson and
    details directories with relevant files called {expid}.json/xlsx.

    Use either source_xlsx to specify a single XLSX with all relevant info
    or

    source_dir expects that {source_dir}/expjson/{expid}.json
    and {source_dir}/details/{expid}.xlsx both exist."""
    def __init__(self, expid='', source_dir='', source_xlsx=''):
        if not source_dir and not source_xlsx:
            raise RuntimeError('One of source_dir or source_xlsx must be specified')

        # parse metadata from multiple files in a directory
        if source_dir:
            d = {}
            d['opt']  = load_json(
                os.path.join(
                    source_dir,
                    "expjson",
                    f"{expid}.json"
                )
            )

            fn = os.path.join(
                source_dir,
                'details',
                f'{expid}.xlsx'
            )

            d['exp'] = pd.read_excel(fn, index_col=0, header=None).iloc[:, 0]
            d['smp'] = pd.read_excel(fn, sheet_name='Sample details').drop_duplicates('Sample').set_index('Sample')
        # parse metadata from single Excel workbook.
        else:
            d = {}
            wb = AnalysisWorkbook(source_xlsx)
            d['opt'] = wb.expd
            d['smp']  = wb.wb['Sample details'].drop_duplicates('Sample').set_index('Sample')
            d['exp']  = wb.wb['Experiment details']
        self.d = d
        self.opt = d['opt']
        self.options = d['opt']
        self.smp = d['smp']
        self.sample_details = d['smp']
        self.exp = d['exp']
        self.experiment_details = d['exp']
        self.exp_name = self.exp['Experiment name']
        self.exp_id = self.exp_name


# libname_libmapfile = {
#     'Transomics Kinase': 'transomics_kinase_v2.mapping.1.txt',
#     'E3 ligases': 'hsLib1608RING_NoDupes.mapping.3.txt',
#     'Kinases V1': 'kinase_library_1.mapping.1.txt',
#     'Kinases V2 (Transomics)':'transomics_kinase_v2.mapping.2.txt',
#     'DUB/E2 V1': 'DUB_guides_no_dupes.mapping.1.txt',
#     'Brunello': 'broadgpp_brunello.mapping.2.txt',
#     'DDR Transomics': 'transomics_ddr_v2-1.mapping.1.txt',
#     'DUB/E2 V2': 'DUB_library2.with_qc.2_ntgenes.mapping.2.txt',
#     'Gattinara': 'broadgpp_gattinara.mapping.1.txt',
#     'Kosuke (human)': 'kosuke_yusa_library_V1.mapping.2.txt',
#     'Yusa (mouse) V2': 'Yusa mouse V2.txt',
#     'lncRNA+coding library K562': 'lncRNA_wg_mixed_no_dupes.2.mapping.2.txt',
#     'lncRNA+coding K562': 'lncRNA_wg_mixed_no_dupes.2.mapping.2.txt',
#     'TKOv2':'Durocher_TKOv2.txt',
#     'TKOv3':'Durocher_TKOv3.txt',
# }

def construct_data_table(name_expath:Union[Dict[str,str],List[str]],
                         metadata:List[Metadata],
                         gene_mappers:Dict[str,Dict]=None,
                         outdir=None):
    """Write score and FDR CSV for analyses in name_expath.

    Expected dir structures: source_dir/experiment_id/analysis_version/etc.

    Args:
        name_expath: name of experiments and directory containing the results.
                If they are the same just pass a list.
        metadata: List of Metadata, matched to name_expath by Metadata.exp_id
        gene_mappers: Dict of dict, for each expid inner dict has map of gene
                symbols.
        file_prefix: prefix with which files were written. ExpID by default.
        outdir: Where consolidated data files will be written.
        """

    if type(name_expath) is list:
        name_expath = {n:n for n in name_expath}

    results = {}

    for mdat in metadata:
        expid = mdat.exp_id

        d = mdat.d
        opt = d['opt']
        ran_analyses = set([a['method'] for a in opt['analyses']])

        prepaff = os.path.join(
            name_expath[expid],
            opt['analysis_version'],
        )

        def tab_to_score_sig(tab, scorekey) -> (pd.DataFrame, pd.DataFrame):
            if gene_mappers:
                symbol_map = gene_mappers[expid]
                idx = tab.index.map(symbol_map)
                if idx.isna().any():
                    raise RuntimeError(f'Gene symbol mapping failed for experiment {expid}' )
                tab.index=idx
            tab = tab.loc[~tab.index.duplicated()]
            scoretab = tab.xs(scorekey, 1,1)
            sigtab = tab.xs('fdr', 1,1)
            cols =  scoretab.columns.map(lambda x: f"{expid}.{x}")
            scoretab.columns = cols
            sigtab.columns = cols
            return (scoretab, sigtab)

        for analysis in ran_analyses:
            short_name = analysis_short_names[analysis]
            if short_name+'_score' not in results:
                results[short_name+'_score'] = []
                results[short_name+'_fdr'] = []

            paff = os.path.join(
                prepaff,
                'tables',
                f"{opt['file_prefix']}.{analysis}_table.csv"
            )

            score, sig = tab_to_score_sig(
                pd.read_csv(paff, header=[0,1], index_col=0),
                {'mageck':'lfc', 'drugz':'normZ'}[analysis]
            )
            results[f"{short_name}_score"].append(score)
            results[f"{short_name}_fdr"].append(sig)

    results = {k:pd.concat(tabs, axis=1) for k, tabs in results.items()}
    if outdir is not None:
        for k, tab in results.items():
            tab.to_csv(os.path.join(outdir, f"{k}.csv"))

    return results

def consolidate_metadata(metadata:List[Metadata],
                         outdir=None):
    """Write comparisons_metadata.csv and exeriments_metadata.csv. List of
    Metadata objects generated from details Excel workbooks."""

    comparisons_table = []
    experiment_metadata = []

    def refmt_ko(s):
        if (s == 'WT'):
            return 'WT'
        #return ' + '.join([gene+'-KO' for gene in s.split('|')])
        return '+'.join(s.split('|'))+'-KO'

    for mdat in metadata:
        deets = mdat.d
        opt = deets['opt']
        smp = deets['smp']
        exp = deets['exp']
        exp_id = deets['opt']['experiment_id']
        smp.loc[:, 'KO'] = smp.KO.fillna('WT')
        smp.loc[smp.KO == '', 'KO'] = 'WT'

        print(exp_id)
        # Shorter name for this frequently used key...
        exp['ExpID'] = exp['Experiment name']
        # But let's keep the old one so that all the fields in the
        #   Excel are present in the dict
        #del exp['Experiment name']
        experiment_metadata.append(
            exp.to_dict()
        )

        # iterate through control groups, comparisons
        for grp_name, ctrlgrp in opt['control_groups'].items():

            # get the analyses ran on this control group, we use the short name here
            ran_analyses = sorted(list(set(
                [analysis_short_names[a['method']] for a in opt['analyses'] if grp_name in a['groups']]
            )))

            for ctrl, treats in ctrlgrp.items():
                for treat in treats:
                    comp_row = {}
                    comp_row['Ctrl samp'] = ctrl
                    comp_row['Treat samp'] = treat
                    comp_row['Timepoint'] = timepoint = grp_name.split('_')[0]

                    # get treatments string, multiple treats joined by "&".
                    #   Treatments are defined by different values appearing in the
                    #   ctrl vs treat sample KO or Treatment columns
                    treatments = [] # list of observed treatments
                    for col in 'KO', 'Treatment':
                        tc, cc = smp.loc[treat, col], smp.loc[ctrl, col]
                        if col == 'KO':
                            tc, cc = refmt_ko(tc), refmt_ko(cc)

                        # If there's not a treatment continue (assuming that absense of a chem/KO
                        #   will never be a treatment
                        if pd.isna(tc) or tc == 'None' or tc == 'DMSO' or tc == 'WT' or tc=='':
                            continue
                        if (tc != cc) or (timepoint != 'endpoints'):
                            # if the background is neutral just list treatment
                            if pd.isna(cc) or (cc == 'None') or (cc == 'DMSO') or \
                                    (cc == 'WT') or (cc == '') or (cc == tc):

                                treatments.append(tc)
                            else:
                                treatments.append(f"{cc} ➤ {tc}")

                    comp_row['Treatment'] = ' & '.join(treatments)

                    for k in ('Dose', 'Growth inhibition %', 'Days grown',
                              'Cell line',  'KO'):
                        comp_row[k] = smp.loc[treat, k]

                    comp_row['Experiment ID'] = exp_id
                    comp_row['Library'] = exp['Library']
                    comp_row['Comparison ID'] = f"{exp_id}.{ctrl}-{treat}"
                    comp_row['Available analyses'] = '|'.join(ran_analyses)

                    comp_row['Control group'] = grp_name

                    comparisons_table.append(comp_row)

    # this is the table that will be displayed when comparisons are selected
    # An experiment focused table with, e.g. author and notes will perhaps come later.
    columns = ['Comparison ID',
               'Experiment ID',
               'Ctrl samp',
               'Treat samp',
               'Timepoint',
               'Treatment',
               'Dose',
               'Growth inhibition %',
               'Days grown',
               'Cell line',
               'KO',
               'Library',
               'Available analyses',
               'Control group',]

    # incase I forget to update columns if adding new
    assert  all([c in columns for c in comp_row.keys()]) and all([c in comp_row.keys() for c in columns])

    comparisons_table = pd.DataFrame(comparisons_table, columns=columns).set_index('Comparison ID')
    experiment_metadata = pd.DataFrame(experiment_metadata).set_index('ExpID')

    try:
        experiment_metadata = experiment_metadata.drop(
            [np.nan], axis=1)
    except:
        pass

    if outdir:
        comparisons_table.to_csv(
            os.path.join(outdir, 'comparisons_metadata.csv'),
            encoding='utf-8-sig'
        )
        experiment_metadata.to_csv(
            os.path.join(outdir, 'experiments_metadata.csv'),
            encoding='utf-8-sig'
        )

    return comparisons_table



# def get_treatment_str(comp_row):
#
#     def is_nt(s):
#         return pd.isna(s) or (s == 'None') or (s == 'DMSO') or (s=='')
#
#     t_deets, c_deets = comp_row.loc[treat], comp_row.loc[ctrl]
#
#     # get chem treat str
#     if not is_nt(t_deets.Treatment):
#         if is_nt(c_deets.Treatment):
#             chem_str = t_deets.Treatment
#         else:
#             chem_str = f"{c_deets.Treatment}{ARROW}{t_deets.Treatment}"
#     # the chem_str is only blank when there was non-genetic perturbation
#     else:
#         chem_str = ''
#
#     if not t_deets.KO:
#         ko_str = ''
#     else:
#         if not c_deets.KO:
#             ko_str = t_deets.KO+'-KO'
#         else:
#             ko_str = f"{c_deets.KO}-KO{ARROW}{t_deets.KO}-KO"
#
#     if chem_str and not ko_str:
#         treatment_str = chem_str
#     elif ko_str and not chem_str:
#         treatment_str = ko_str
#     else:
#         if t_deets.Treatment == c_deets.Treatment:
#             # it's a KO in a chemical background
#             treatment_str = f"{ko_str} (with {chem_str})"
#         else:
#             treatment_str = f"{chem_str} (in {ko_str} cells)"
#     return treatment_str



if __name__ == '__main__':

    os.chdir('/Users/johnc.thomas/Dropbox/crispr/almu_review')

    expid = 'Alvarez-QuilonDurocher2020'

    metadats = []
    exppaths = {}
    gene_mappers = {}
    mdat = Metadata(source_xlsx=f'details/{expid}.xlsx')
    exppaths[expid] = f'ran_analyses/{expid}'

    # Analysis version _shouldn't_ be fixed in any metadata, but it needs
    #    restating all the time. This feels messy.
    mdat.d['opt']['analysis_version'] = 'db1'
    mdat.d['opt']['file_prefix'] = 'x'
    metadats.append(mdat)

    mapper = f"symbol_mapping/{mdat.d['exp'].Library}.csv"
    gene_mappers[expid] = load_symbol_map(mapper)


    metadata = consolidate_metadata(
        metadats,
        outdir=None,
    )
    print(metadata.Treatment)
