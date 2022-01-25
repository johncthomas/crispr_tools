#!/usr/bin/env python3
import json
import logging, os, pathlib, datetime, inspect, argparse, sys

from typing import List, Dict
from copy import copy
from subprocess import call, check_output
from pathlib import Path
from itertools import combinations


from attrdict import AttrDict

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from crispr_tools import qc, tools, jacks_tools
import pprint

try:
    from jacks.jacks_io import runJACKS
    from jacks.infer import LOG as jacksLOG
    jacksLOG.setLevel(logging.WARNING)
except ImportError:
    # print('To run jacks you need to install JACKS,\n', 'https://github.com/felicityallen/JACKS/tree/master/jacks\n'
    #       "You can still run Mageck though, if it's installed.")
    def runJACKS(*a, **k):
        raise ModuleNotFoundError('JACKS not installed!!!')

from crispr_tools.drugz import drugZ_analysis
from crispr_tools.tools import list_not_str

# with open(pathlib.Path(__file__).parent/'version.txt') as f:
#     __version__ = f.readline().replace('\n', '')

class PipelineOptionsError(Exception):
    """Errors in the configuration file that would prevent the pipeline from running"""
    pass

pipeLOG = logging.getLogger('pipeline')
pipeLOG.setLevel(logging.INFO)
#logging.getLogger('matplotlib').setLevel(logging.WARNING)

#from crispr_tools.count_reads import count_reads, count_batch, map_counts
from crispr_tools.tools import plot_read_violins, plot_ROC, plot_volcano, tabulate_mageck, tabulate_drugz
from crispr_tools.jacks_tools import tabulate_score, scores_scatterplot, mahal_nocov
from crispr_tools.version import __version__


#todo: refactor the analysis so it's modulor and not just a massive function
#todo check validity of every file and option before starting.
#  eg:
#   do all reps exsist in the count file
#   are samples in controls defined in the sample reps
#todo handle different analyses better, making it easy to add new ones and turn of during analysis
#todo Either run_{method} or call_{method}
    # in particular there should be functions for each that take the exact same arguments.





def clargs_to_dict(argstring) -> Dict[str, str]:
    """process command line args of format:
    --arg1 val --flag --arg2 val2

    returns {'arg1':'val', 'flag':'', 'arg2':'val2'}"""
    kwargs = {}
    while '  ' in argstring:
        argstring = argstring.replace('  ', ' ')
    splitargs = argstring.split()
    for i in range(len(splitargs)):
        thing = splitargs[i]
        if thing.startswith('-'):
            key = thing.replace('-', '')
            if i+1==len(splitargs) or splitargs[i+1].startswith('-'):
                kwargs[key] = ''
            else:
                kwargs[key] = splitargs[i+1]
    return kwargs

#todo test excel kwargs

class AnalysisWorkbook:
    """Parse the excel sheet describing how an experiment should be analysed.
    Generates the dict for JSON or calling .run_analyses

    Default method kwargs:
        mageck: --remove-zero both --remove-zero-threshold 2*pseudocount

    additional_options overwrites anything in the excel file"""

    def __init__(self, xlfn, parse_workbook=True, **additional_options):
        #self.xpk = 'Experiment details'
        self.wb = self.load_analysis_workbook(xlfn)
        if parse_workbook:
            self.expd = self.workbook_to_dict(**additional_options)

    def varefy_wb_integrity(self):
        assert not self.wb['Sample details'].index.duplicated().any()

    def load_analysis_workbook(self, fn):
        # nothing needs to be numeric
        wb = pd.read_excel(fn, sheet_name=None, dtype=object, )
        wb = {k:tab.fillna('') for k,tab in wb.items()}

        wb['Sample details'] = wb['Sample details'].set_index('Replicate')

        wb['Experiment details'] = pd.read_excel(
            fn, sheet_name='Experiment details',
            header=None, index_col=0,
        )[1].fillna('')

        #wb['Sample details'] = wb['Sample details'].set_index('Replicate')

        return wb

    def workbook_to_dict(self, **additional_options):
        # couple of options are cells that can contain comma split values
        # pretty easy to inadvertently slip a space in there...
        safesplit = lambda x: x.replace(' ', '').split(',')

        repdeets = self.wb['Sample details']
        ctrlgrps = self.wb['Control groups']

        # all "Experiment details" k:v are top level and don't need processing
        experiment_dictionary = {}

        deets = self.wb['Experiment details'].to_dict()
        for wbk, jsonk in {'Experiment name': 'experiment_id',
                       'Analysis version': 'analysis_version',
                       'Notes': 'notes',
                       'File prefix': 'file_prefix'}.items():
            try:
                experiment_dictionary[jsonk] = deets[wbk]
            except KeyError:
                # Most won't exist and any other errors
                #   will get picked up later
                pass

        # control groups. Results in:
        # {grpA:{ctrl_samp:[test_samp1, test_samp2, ...], ...}, ...}
        ctrl_dict = {}
        for grpn in ctrlgrps.Group.unique():
            cgrp = ctrlgrps[ctrlgrps.Group == grpn]
            # groupby gives the indexes, we want the Test sample values
            # gotta list() the reps cus json.dump can't handle an ndarray
            g = {c:list(cgrp.loc[i, 'Test sample'].values)
                 for c, i in cgrp.groupby('Control sample').groups.items()}
            ctrl_dict[grpn] = g
        experiment_dictionary['control_groups'] = ctrl_dict

        # sample_reps, just {sampA:[repA1, repA2], ...}
        experiment_dictionary['sample_reps'] = {
            k:list(v.values) for k, v in repdeets.groupby('Sample').groups.items()
        }

        # Analyses. Producing a list of dicts with "method" and "groups":List
        #   required keys. Should practically always be a fn_counts, but ignore
        #   if blank. kwargs & pseudocount added if not present
        analyses = []
        for _, row in self.wb['Analyses'].iterrows():
            for method in safesplit(row['Method']):
                andict = {'method':method,
                          'groups':safesplit(row['Control group'])}

                if 'kwargs' not in andict:
                    andict['kwargs'] = {}

                if row['Add pseudocount']:
                    andict['pseudocount'] = int(row['Add pseudocount'])
                else:
                    andict['pseudocount'] = 1
                if row['Counts file']:
                    andict['counts_file'] = row['Counts file']

                # add default args if any
                if not row['Arguments']:
                    if method == 'mageck':
                        pc = andict['pseudocount']
                        if pc > 1:
                            cutoff = 2*pc
                        else:
                            cutoff = 0
                        andict['kwargs'] = {'remove-zero':'both',
                                            'remove-zero-threshold':cutoff}
                else:
                    rowargs = row['Arguments']
                    try:
                        kwargs = eval(rowargs)
                        if not type(kwargs) is  dict:
                            raise RuntimeError()
                    except:
                        # try parsing as command line options
                        kwargs = clargs_to_dict(rowargs)
                    andict['kwargs'] = kwargs

                # deal with paired, which is handled different in the programs
                if row['Paired'] and (method == 'mageck'):
                    andict['kwargs']['paired'] = ''
                elif not row['Paired'] and (method == 'drugz'):
                    andict['kwargs']['unpaired'] = True

                analyses.append(andict)
        experiment_dictionary['analyses'] = analyses

        for k, v in additional_options.items():
            experiment_dictionary[k] = v

        if 'file_prefix' not in experiment_dictionary:
            experiment_dictionary['file_prefix'] = experiment_dictionary['experiment_id']

        return experiment_dictionary



class StreamToLogger(object):
   """
   File-like stream object that redirects writes to a logger instance.
   """
   def __init__(self, logger, log_level=logging.INFO):
      self.logger = logger
      self.log_level = log_level
      self.linebuf = ''

   def write(self, buf):
      for line in buf.rstrip().splitlines():
         self.logger.log(self.log_level, line.rstrip())



# make the jacks repmap: rep, samp, ctrl
def call_jacks(sample_reps:Dict[str,List[str]], ctrlmap, count_fn, prefix, kwargs:Dict=None):
    if kwargs is None:
        kwargs={}
    repmap_fn = prefix+'.repmap.tsv'
    write_repmap(sample_reps, ctrlmap, repmap_fn)
    #pipeLOG.info(f"Running JACKS, mode and median normalised with {repmap_fn}")

    jacks_args = (count_fn, repmap_fn, count_fn, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',)

    runJACKS(*jacks_args, outprefix=prefix, **kwargs)


def set_logger(log_fn):

    hndlr = logging.FileHandler(log_fn, 'w')
    # hndlr.setLevel(logging.INFO)
    pipeLOG.setLevel(logging.INFO)
    pipeLOG.addHandler(hndlr)
    try:
        jacksLOG.addHandler(hndlr)
    except:
        pass
    #mageckLOG.addHandler(hndlr)


def write_repmap(sample_reps:Dict[str, list], ctrlmap:Dict[str, list], repmap_fn:str):
    """Write a JACKS format replicate file.

    sample_reps link the replicate names in the count file to given sample
    names. ctrlmap specifies which samples are controls for which other samples."""
    reptab = [['rep', 'samp', 'ctrl']]
    for ctrl_samp, treat_samps in ctrlmap.items():
        for ctrl_rep in sample_reps[ctrl_samp]:
            reptab.append([ctrl_rep, ctrl_samp, ctrl_samp])
        for treat in treat_samps:
            for rep in sample_reps[treat]:
                reptab.append([rep, treat, ctrl_samp])
                # print([rep, samp, ctrl])

    with open(repmap_fn, 'w') as f:
        for line in reptab:
            f.write('\t'.join(line) + '\n')


def call_drugZ_batch(sample_reps:Dict[str, list],
                     control_map:Dict[str, list],
                     counts_file:str,
                     prefix:str,
                     kwargs:dict=None,
                     pseudocount=1,):
    """output files written to {file_prefix}.{ctrl}-{treat}.tsv

    One file per comparison, in the default drugz format. Use tabulate_drugz to
    get a standardised table"""
    kwargs = kwargs
    if kwargs is None:
        kwargs = {}

    dzargs = AttrDict()
    dzargs.infile = counts_file
    # defaults
    dzargs.index_column = 0
    dzargs.minobs = 1
    dzargs.half_window_size = 500
    dzargs.quiet = False
    dzargs.pseudocount = pseudocount
    dzargs.fc_outfile=''
    dzargs.unpaired = False

    # not implimented in drugZ at time of writing
    dzargs.remove_genes = None
    dzargs.update(kwargs)

    for ctrl_samp, treat_samples in control_map.items():

        dzargs.control_samples = ','.join(
            list_not_str( sample_reps[ctrl_samp] )
        )

        treat_samples = list_not_str(treat_samples)
        for treat_samp in treat_samples:

            dzargs.drug_samples = ','.join(
                list_not_str(sample_reps[treat_samp])
            )

            dzargs.drugz_output_file = f"{prefix}.{ctrl_samp}-{treat_samp}.tsv".replace('//', '/')

            drugZ_analysis(dzargs)

    pipeLOG.info('Finished drugZ')


def call_mageck(control_samp:str, treat_samp:str, sample_reps:Dict[str, List[str]],
                counts_file:str, prefix:str, kwargs:List[Dict], logger=None, dryrun=False):
    """Call mageck from the command line:
        "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}{ctrlnm}-{sampnm}"

    Args:
        control_samp: Name of control sample as written in sample_reps.
        treat_samp:   as above, but for treatment.
        sample_reps: Dict giving mapping of sample to replicates, replicate names match
            column headers in counts file.
        prefix: Prefix added to all created files, should include desired output dir.
        kwargs: Keyword arguments to be passed on to mageck. If kwarg is flag it should
            be in the format {'flag-kwarg':''}

    """
    mageck_str = "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}.{ctrlnm}-{sampnm}"
    # get the replicate strings
    ctrls = ','.join(sample_reps[control_samp])
    treats = ','.join(sample_reps[treat_samp])
    s = mageck_str.format(
        counts=counts_file,
        treat=treats,
        ctrl=ctrls,
        outprefix=prefix,
        ctrlnm=control_samp,
        sampnm=treat_samp
    )
    mag_additional_args = []
    if kwargs:

        mag_additional_args = [f"--{_k} {_v}" for _k, _v in kwargs.items()]

    mag_args = s.split() + mag_additional_args
    # for some reason mageck fails to understand mag_args if you don't use shell=True
    command = ' '.join(mag_args)
    if logger:
        logger.info(command)
    print(command)
    if not dryrun:
        call(command, shell=True)


# CALL SIGNATURE MUST MATCH CALL JACKS and any others
#todo one function/class as wrapper for all analysis methods
def call_mageck_batch(sample_reps:Dict[str, list],
                      control_map:Dict[str, list],
                      counts_file:str,
                      prefix:str,
                      kwargs:dict=None,
                      pseudocount=1,):

    """Run mageck analyses using comparisons specified in control_map."""

    mageck_pairs_done = [] # used to avoid repeating comps in EXTRA
    call("which mageck".split())
    pipeLOG.info('Running MAGeCK version ' + check_output(["mageck", "-v"]).decode())

    # Deal with pseudocount. Mageck doesn't have this facility, so we need
    #   to write a temporary counts file.
    # It adds 1 before logging anyway.
    tmp_fn = None
    if pseudocount > 1:
        counts = pd.read_csv(counts_file, sep='\t', index_col=0)
        num_cols = counts.dtypes != object
        counts.loc[:, num_cols] += pseudocount
        tmp_fn = prefix+'temp_counts.tsv'
        counts.to_csv(tmp_fn, sep='\t')
        counts_file = tmp_fn

    for ctrl_samp, treat_samples in control_map.items():

        if type(treat_samples) is str:
            treat_samples = [treat_samples]

        for treat in treat_samples:
            if treat == ctrl_samp:
                continue

            call_mageck(ctrl_samp, treat, sample_reps, counts_file, prefix, kwargs)
            mageck_pairs_done.append((ctrl_samp, treat))

    # delete the file with additional pseudocount, if there is one.
    if tmp_fn is not None:
        os.remove(tmp_fn)

# When adding new functions, need to deal with pseudocount option in the main func
#   and pairedness in the AnalysisWorkbook
analysis_functions = {'mageck':call_mageck_batch, 'drugz':call_drugZ_batch}
analysis_tabulate = {'mageck':tabulate_mageck, 'drugz':tabulate_drugz}
available_analyses = analysis_functions.keys()


# def iter_comps(comparisons: List[dict], tab: pd.DataFrame=None):
#     for comparison in comparisons:
#         for ctrl_samp, comp_samps in comparison.items():
#             for comp in comp_samps:
#                 if tab is not None:
#                     if comp not in tab.columns or ctrl_samp not in tab.columns:
#                         continue
#                 yield ctrl_samp, comp


def process_control_map(controls, samples):
    """Get mappings ctrl_samp->[test_samps,] for each control group,
    dealing with keywords ALL & EXCEPT.

    Returns a copy of the input controls dict."""

    controls = copy(controls)
    for ctrl_grp, ctrlmap in controls.items():
        for ctrl, samps in ctrlmap.items():
            if type(samps) == str:
                samps = [samps]
            # if you add keywords you'll need to add them to the database add_data module
            if samps == ['ALL']:
                samps = [s for s in samples if s!=ctrl]
            elif type(samps) == dict:
                # maybe could have other keywords in the future
                samps = [s for s in samples if s not in samps['EXCEPT']]
            ctrlmap[ctrl] = samps
    return controls


def check_multicontrol_to_treat(control_map):
    """Check if a single treatment sample is paired to multiple controls
    (which is not allowed by JACKS at least).

    Returns:
        {control_group:{treat:[control1, control2]}}
        for any conflicts are found. An empty dict if there's no problems."""
    # store conflicts by treat:[ctrl1, ctrl2, ...]
    conflicts = {grp:{} for grp in control_map.keys()}
    conflict_found = False
    for grp, controls in control_map.items():
        for a, b in combinations(controls.keys(), 2):

            for treatSamp in controls[a]:
                if treatSamp in controls[b]:
                    conflict_found = True
                    try:
                        conflicts[grp][treatSamp].append(a)
                    except KeyError:
                        conflicts[grp][treatSamp] = [a]
                    conflicts[grp][treatSamp].append(b)
    if conflict_found:
        return {grp:cnflct for grp, cnflct in conflicts.items() if cnflct}

def validate_required_arguments(arguments:dict):
    """Check for required options, the basic format of `analyses` options,
    and that sample names match in control_groups and sample_reps.

    If validation fails, raise a RuntimeError"""



    # check the required top level options are all present
    required_options = ['sample_reps', 'experiment_id', 'analysis_version',
                        'file_prefix', 'control_groups', 'analyses']
    option_is_missing = lambda op: (op not in arguments.keys()) or (not arguments[op])
    missing_options = [op for op in required_options if option_is_missing(op)]
    if missing_options:
        raise RuntimeError(
            f'Required options missing (or valueless) from config file: {",".join(missing_options)}'
        )


    # check the required analyses options present
    required_analysis_opts = ['method']
    for k in required_analysis_opts:
        for ans in arguments['analyses']:
            if k not in ans:
                raise RuntimeError(
                    f'Required analysis {k} option missing from analysis {ans}.'
                )

    # and the analyses dics themselves are formatted correctly.
    analysis_opts_types = {'method':str, 'kwargs':dict, 'groups':list, 'counts_file':str,
                           'label':str}
    for analysis_dict in arguments['analyses']:
        for opt, tipe in analysis_opts_types.items():
            try:
                if type(analysis_dict[opt]) != tipe:
                    raise RuntimeError(
                        f"Analysis option types not as expected: {arguments['analyses']}"
                    )
            # there are optional options, ignore them if they're missing
            except KeyError:
                pass

    samples = arguments['sample_reps'].keys()
    missing_samples = set()
    for _, ctrl_samps in arguments['control_groups'].items():
        for ctrl, samps in ctrl_samps.items():
            if ctrl not in samples:
                missing_samples.add(ctrl)
            for smp in samps:
                if smp not in samples:
                    missing_samples.add(smp)
    missing_samples = list(missing_samples)
    if missing_samples:
        raise RuntimeError(f"Samples named in control_groups missing in sample_reps:\n{', '.join(missing_samples)}")

def process_arguments(arguments:dict) -> dict:
    """Deal with special keywords from the experiment dictionary.
    Varefy integrity of options."""
    # (currently optional arguments are handled in the pipeline, but
    #  alternatively I could have them set here
    #  so changes to how the arguments
    #  work would require only changingn things here.
    #  This is probably fine for now.)

    # raise exception for invalid arguments
    validate_required_arguments(arguments)

    samples = arguments['sample_reps'].keys()
    controls = process_control_map(arguments['control_groups'], samples)
    arguments['control_groups'] = controls

    # generate output_dir from exp_id and analysis_version
    arguments["output_dir"] = os.path.join(arguments['output_dir'],
                                           arguments["experiment_id"],
                                           arguments["analysis_version"])

    pipeLOG.info(f'output_dir = {arguments["output_dir"]}')

    for k, v in arguments['sample_reps'].items():
        if type(v) == str:
            arguments['sample_reps'][k] = [v]

    reps_list = []
    for r in arguments['sample_reps'].values():
        reps_list.extend(r)

    # get list of all counts files used
    counts_files = set()
    counts_not_top_level = False
    try:
        # this top level not required
        counts_files.add(arguments['counts_file'])
    except KeyError:
        counts_not_top_level = True

    # add analysis specific counts files, also check method exists
    for ans in arguments['analyses']:
        if 'counts_file' in ans:
            counts_files.add(ans['counts_file'])
        else:
            if counts_not_top_level:
                raise PipelineOptionsError('No counts file set for analysis')
        method = ans['method']
        if method not in available_analyses:
            raise PipelineOptionsError(f'Unknown method: "{method}" specified.')


    # Check a counts file has been specified
    counts_files = list(counts_files)
    if not counts_files:
        raise PipelineOptionsError('No counts file set for analysis')

    for cfn in counts_files:
        with open(cfn) as f:
            line = next(f)
            if not '\t' in line:
                raise PipelineOptionsError(
                    f"No tabs in first line of count file {cfn}"+
                    f", is it comma seperated?\n\t{line}"
                )
            # check for missing/miss-spelled replicates
            cnt_reps = line.strip().split('\t')
            missing_reps = [r for r in reps_list if r not in cnt_reps]
            if missing_reps:
                raise ValueError(
                    f"These replicates not found in count file: {missing_reps}"
                    f"\nCount columns: {', '.join(cnt_reps)}"
                )

    # Control groups to be included, defaults to all if None
    if 'run_groups' not in arguments:
        arguments['run_groups'] = None
    run_groups = arguments['run_groups']
    if run_groups is not None:
        if type(run_groups) is str:
            arguments['run_groups'] = run_groups.split(',')



    return arguments


def run_analyses(output_dir, file_prefix,
                 sample_reps:Dict[str, list],
                 control_groups:Dict[str, Dict[str, list]],
                 analyses:List[dict],
                 methods_kwargs:Dict=None,
                 dont_log=False, #todo replace dont_log with some kind of verbosity thing
                 compjoiner='-', #tools.ARROW,
                 use_group_as_label=False,
                 notes='',
                 skip_analyses=None,
                 run_groups:List[str]=None, counts_file=None):

    """Run batches of CRISPR analyses."""

    #pipeLOG.warning(f'Unrecognised arguments passed to run_analyses:\n   {unrecognised_kwargs}')
    pipeLOG.info(f"notes: {notes}")

    if skip_analyses is None:
        skip_analyses = []

    if run_groups is None:
        group_included = lambda x: True
    else:
        group_included = lambda x: x in run_groups
        pipeLOG.info(f'Running only control groups: {run_groups}')

    # Create the root experiment directory
    # Can't make directory and subdir at the same time, so iterate through the tree
    p = str(output_dir)
    for i in range(len(p)):
        #print(p)
        d = p.split('/')[:i+1]
        try:
            os.mkdir('/'.join(d))
        except FileExistsError:
            pass

    if not dont_log:
        t = '{}-{}-{}_{}h{}m{}s.{}'.format(*datetime.datetime.now().timetuple()[:-1])
        set_logger(str(Path(output_dir, file_prefix + f'log_{t}.txt')))

    # pull the args from the function
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    del frame # help garbage collection
    pipeLOG.info('\n\t'.join([f"{x} = {values[x]}" for x in args]))
    pipeLOG.info(f"Full output directory = {os.path.realpath(output_dir)}")

    # create the directories
    methods_used = list(set([ans['method'] for ans in analyses]))

    for analysis_method in methods_used:
        call(['mkdir', '-p', str(Path(output_dir, analysis_method))])
        #call(['mkdir', str(Path(output_dir, analysis_method, 'tables'))])
        call(['mkdir', str(Path(output_dir, analysis_method, 'files'))])
    table_dir = str(Path(output_dir, 'tables'))
    call(['mkdir', table_dir])

    ######################
    ## Run the analyses ##
    ######################
    ran_analyses = set()
    for analysis_dict in analyses:
        analysis_method = analysis_dict['method']
        if analysis_method in skip_analyses:
            pipeLOG.info(f"Skipping analysis method {analysis_method}")
            continue

        # The structure here (which may be pointlessly complicated) is that
        #   default arguments for an analysis method can be set at the top
        #   level of the JSON, but individual analysis_dicts may contain
        #   arguments that override these.
        try:
            default_method_kwargs = methods_kwargs[analysis_method]
        except:
            default_method_kwargs = {}

        # if we do not specify groups, use all control groups.
        if 'groups' not in analysis_dict:
            groups = {g:{} for g in control_groups.keys()}
        else:
            groups = analysis_dict['groups']

        try:
            if analysis_dict['kwargs']:
                curr_kwargs = analysis_dict['kwargs']
            else:
                curr_kwargs = default_method_kwargs
        except KeyError:
            curr_kwargs = default_method_kwargs

        # go through the selected groups for this analysis and run the thing
        for grp in groups:
            if not group_included(grp):
                pipeLOG.info(f"Skipping group {grp}")
                continue
            if use_group_as_label:
                labstr = grp
            if 'label' in analysis_dict:
                labstr = '.'+analysis_dict['label']
            else:
                labstr = ''

            ctrl_map = control_groups[grp]
            try:
                curr_counts = analysis_dict['counts_file']
            except KeyError:
                curr_counts = counts_file

            pipeLOG.info(
                f"Running batch {analysis_method}\n\tgroup: {grp}\n\twith options: {analysis_dict}\n"
                f"\tWith kwargs: {curr_kwargs}"
            )

            get_prefix = lambda out_type: str(Path(output_dir, analysis_method, out_type, f"{file_prefix}{labstr}"))
            # all these things must have the same signature
            out_prefix = get_prefix('files')

            # info for tabulating results after the fact.
            ran_analyses.add(
                (analysis_method, out_prefix, f"{file_prefix}{labstr}")
            )
            analysis_func = analysis_functions[analysis_method]
            analysis_func(sample_reps, ctrl_map, curr_counts, out_prefix, curr_kwargs, )

        # tabulate the analyses
        for analysis_method, results_prefix, table_file_prefix in list(ran_analyses):
            tab = analysis_tabulate[analysis_method](results_prefix, compjoiner)
            tabfn = str(os.path.join(output_dir, 'tables', f'{table_file_prefix}.{analysis_method}_table.csv'))
            pipeLOG.info(f'writing table: {tabfn}')
            tab.to_csv(tabfn, encoding='utf-8-sig')


if __name__ == '__main__':

    print(__version__)
    parser = argparse.ArgumentParser(
        description="Run mageck and jacks analyses using a JSON file.\n "
                    "Use arguments below to override JSON options"
    )

    parser.add_argument(
        'config_file', metavar='EXCEL_OR_JSON',
        help = ("path to .xlsx or .json file specifying arguments. JSON must must contain `sample_reps`"
                ", `analyses` & `control_groups` keywords and values. If passing an excel, --analysis-version"
                " must be set. "
                "Other options may be overridden by command line arguments.")
        )
    parser.add_argument('--counts', metavar='COUNTS', help='Path to counts file',
                        default=None, dest='counts_file')
    parser.add_argument('--output-dir', metavar='PATH', default='.',
                        help=('Path to where results directory will be created if not current'
                              ' directory. Experiment_id and analysis_version will determine results dir.') )
    parser.add_argument('--file-prefix', default=None,
                        help="String to form identifying prefix for all files generated.")
    parser.add_argument('--skip-analyses', metavar='list,of,progs', default=None,
                        help='Filter analyses to be run.')
    parser.add_argument('--run-groups', metavar='list,of,groups', default=None,
                        help='Specify control groups (as defined in exp dict) to be included.')
    parser.add_argument('--dont-log', action='store_true', dest='dont_log', default=None,
                        help="Don't write a log file.")
    parser.add_argument('--analysis-version', default=None,
                        help='Output files will be stored in a directory of the above name, within the experiment dir.')

    # get arguments from the command line and the YAML
    clargs = parser.parse_args() # need to assign this before calling vars() for some reason
    cmd_line_args = vars(clargs)
    file_prefix = cmd_line_args['file_prefix']
    if file_prefix is None:
         file_prefix = ''

    # load the configuration file, stripping out the "comment" lines
    config_file = cmd_line_args['config_file']
    if config_file.endswith('xlsx'):
        config_file_args = AnalysisWorkbook(config_file).workbook_to_dict()
    elif config_file.endswith('json'):
        json_str = '\n'.join([l for l in open(config_file) if l.replace(' ', '')[0] != '#'])

        def dict_raise_on_duplicates(ordered_pairs):
            """Reject duplicate keys."""
            d = {}
            for k, v in ordered_pairs:
                if k in d:
                    raise ValueError("Duplicate key in JSON: %r" % (k,))
                else:
                    d[k] = v
            return d
        config_file_args = json.loads(json_str, object_pairs_hook=dict_raise_on_duplicates)
    else:
        raise RuntimeError(f"Unrecognised file type: {config_file}. Only accepts .json or .xlsx")

    # over write yml_args with any specified in the command line
    for k, v in cmd_line_args.items():
        if v is not None:
            config_file_args[k] = v

    # validate and process arguments
    args = process_arguments(config_file_args)

    # These don't need to be passed to the pipeline.
    for k in ['config_file', 'notes', 'experiment_id', 'analysis_version',
              'labels']:
        if k in args:
            del args[k]

    run_analyses(**args)

