#!/usr/bin/env python3
import json
import logging, os, pathlib, datetime, inspect, argparse, sys

from typing import List, Dict
from copy import copy
from subprocess import call, check_output
from pathlib import Path
from itertools import combinations

from argparse import Namespace
from attrdict import AttrDict, AttrMap

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

class ConfigurationError(Exception):
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


"""Go from FastQ files to completed JACKS/MAGeCK analysis. 
fq->mapped_counts&violin plots are one command (count_reads.count_batch() ) 
counts->charts&tables another."""



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
                     kwargs:dict=None, ):
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
    dzargs.pseudocount = 1
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
                      skip_extra_mageck=True):

    """Run mageck analyses using comparisons specified in control_map."""

    mageck_pairs_done = [] # used to avoid repeating comps in EXTRA
    call("which mageck".split())
    pipeLOG.info('Running MAGeCK version ' + check_output(["mageck", "-v"]).decode())

    for ctrl_samp, treat_samples in control_map.items():

        if type(treat_samples) is str:
            treat_samples = [treat_samples]

        for treat in treat_samples:
            if treat == ctrl_samp:
                continue

            call_mageck(ctrl_samp, treat, sample_reps, counts_file, prefix, kwargs)
            mageck_pairs_done.append((ctrl_samp, treat))
        # EXTRA
        if not skip_extra_mageck:
            # run all combinations of samples that share at least one control sample
            # replace the control group name in the string
            pref_split = prefix.split('.')
            pref_split[-2] = 'EXTRA'
            prefix_extra = '.'.join(pref_split)
            for c, t in combinations(treat_samples, 2):
                if (c, t) not in mageck_pairs_done and (t, c) not in mageck_pairs_done:
                    mageck_pairs_done.append((c, t))
                    # maybe don't tablulate these at the moment.
                    call_mageck(c, t, sample_reps, counts_file, prefix_extra, kwargs)

def validate_expd(expd:dict):
    pass


analysis_functions = {'mageck':call_mageck_batch, 'drugz':call_drugZ_batch}
analysis_tabulate = {'mageck':tabulate_mageck, 'drugz':tabulate_drugz}
available_analyses = analysis_functions.keys()

def iter_comps(comparisons: List[dict], tab: pd.DataFrame=None):
    for comparison in comparisons:
        for ctrl_samp, comp_samps in comparison.items():
            for comp in comp_samps:
                if tab is not None:
                    if comp not in tab.columns or ctrl_samp not in tab.columns:
                        continue
                yield ctrl_samp, comp


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
    required_options = ['sample_reps', 'experiment_id', 'analysis_version', 'counts_file',
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

def process_arguments(arguments:dict):
    """deal with special keywords from the experiment yaml, and allow some
    ambiguous syntax in the yaml. Also do some checking of validity"""
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

    with open(arguments['counts_file']) as f:
        line = next(f)
        if not '\t' in line:
            raise ConfigurationError(
                f"No tabs in first line of count file {arguments['counts_file']}"+
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

    # turn the csv analyses into an actual list
    for k in 'skip_groups', 'skip_analyses':
        if (k in arguments) and arguments[k]:
            arguments[k] = arguments[k].split(',')
    return arguments


def run_analyses(counts_file, output_dir, file_prefix,
                 sample_reps:Dict[str, list],
                 control_groups:Dict[str, Dict[str, list]],
                 analyses:List[dict],
                 methods_kwargs:Dict=None,
                 dont_log=False,
                 compjoiner='-', #tools.ARROW,
                 use_group_as_label=False,
                 notes='',
                 skip_analyses=None,
                 skip_groups=None):

    """Run batches of CRISPR analyses."""

    #pipeLOG.warning(f'Unrecognised arguments passed to run_analyses:\n   {unrecognised_kwargs}')
    pipeLOG.info(f"notes: {notes}")

    if skip_analyses is None:
        skip_analyses = []
    if skip_groups is None:
        skip_groups = []

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
        try:
            default_method_kwargs = methods_kwargs[analysis_method]
        except:
            default_method_kwargs = {}

        # if we do not specify groups, use all control groups.
        if 'groups' not in analysis_dict:
            groups = {g:{} for g in control_groups.keys()}
        else:
            groups = analysis_dict['groups']

        # These kwargs set at the analysis level, can be overridden by args at
        #   the control group level. New, analysis level, args overwrite experiment
        #   level args
        # try:
        #     curr_kwargs = {**default_method_kwargs, **analysis_dict['kwargs']}
        # except KeyError:
        #     curr_kwargs = copy(default_method_kwargs)
        try:
            if analysis_dict['kwargs']:
                curr_kwargs = analysis_dict['kwargs']
            else:
                curr_kwargs = default_method_kwargs
        except KeyError:
            curr_kwargs = default_method_kwargs


        # go through the selected groups for this analysis and run the thing
        for grp in groups:
            if grp in skip_groups:
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
            ran_analyses.add((analysis_method, out_prefix, f"{file_prefix}{labstr}"))
            analysis_func = analysis_functions[analysis_method]
            analysis_func(sample_reps, ctrl_map, curr_counts, out_prefix, curr_kwargs)

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
        'config_file', metavar='JSON_PATH',
        help = "path to .json file specifying arguments. At a minimum must contain `sample_reps`"
                  ", `analyses` & `control_groups` keywords and values. "
                  "Other options may be overridden by command line arguments."
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
    parser.add_argument('--skip-groups', metavar='list,of,groups', default=None,
                        help='Filter control groups (as defined in exp dict) to be included.')
    parser.add_argument('--dont-log', action='store_true', dest='dont_log', default=None,
                        help="Don't write a log file.")

    # get arguments from the command line and the YAML
    clargs = parser.parse_args() # need to assign this before calling vars() for some reason
    cmd_line_args = vars(clargs)

    # load the configuration file, stripping out the "comment" lines
    config_file = cmd_line_args['config_file']
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

    # over write yml_args with any specified in the command line
    for k, v in cmd_line_args.items():
        if v is not None:
            config_file_args[k] = v

    # validate and process arguments
    args = process_arguments(config_file_args)

    # These don't need to be passed to the pipeline.
    for k in ['config_file', 'notes', 'experiment_id', 'analysis_version',
              'labels']:
        try:
            del args[k]
        # if it doesn't exist...
        except KeyError:
            # we don't need to delete it
            pass

    run_analyses(**args)

