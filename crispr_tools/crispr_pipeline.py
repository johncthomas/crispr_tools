#!/usr/bin/env python3
import json
import logging, os, datetime, inspect, argparse
import gzip

from typing import List, Dict
from copy import copy
from subprocess import call, check_output
from pathlib import Path
from itertools import combinations


from attrdictionary import AttrDict

import pandas as pd
#import matplotlib.pyplot as plt
from crispr_tools import data_classes

ARROW = '→'
# todo don't write drugz tables until everything is run
# have default prefix
# outdir + prefix issue if both not local (tests should include weird combinations)
# Validation should include sample names with dots and dashes!


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

def file_exists_or_zero(fn):
    """Determine whether the output drugz/mageck file already exists and contains
    data. Returns True if the file exists and is not zero bytes. Returns False
    if the file does not exist or it does exist but is zero bytes. Used for
    checkpointing purposes. Analysis steps will skip if the output drugz/mageck file
    exists and is not size zero bytes (function returns True)."""
    if os.path.exists(fn) and os.path.getsize(fn) != 0:
        return True
    else:
        return False

class PipelineOptionsError(Exception):
    """Errors in the configuration file that would prevent the pipeline from running"""
    pass

logging.basicConfig() # required for setLevel to work if you don't set a specific handler
pipeLOG = logging.getLogger('pipeline')
pipeLOG.setLevel(logging.INFO)

#from crispr_tools.count_reads import count_reads, count_batch, map_counts
#from crispr_tools.tools import plot_read_violins, plot_ROC, plot_volcano, tabulate_mageck, tabulate_drugz
from crispr_tools.tools import tabulate_mageck, tabulate_drugz
#from crispr_tools.jacks_tools import tabulate_score, scores_scatterplot, mahal_nocov
from crispr_tools.version import __version__

#todo analysis programs should be defined as classes that do the following:
    # handle kwargs, have a .call method, have a .tabulate_results method,
    # write tables for the screen viewer
    # have phenotype and significance scores available as same attributes,
    # inheret from class that has comparison selection built in
    #   i.e. where you pass two samples and get a result
    # a short name attribute
    # This would make it more easy to add new analyses and mean that
    #   I stop scattering around the parts of what should be a unified concept.

#todo test excel kwargs
#todo clear up how args are processed when going through a workbook is unclear
    # process_workbook (counts_file from cmd line passed separately - what else needs to be passed?)
    # the dictionary from above passed to process args
    # other command line args applied

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



def get_treatment_str(samp_deets, treat, ctrl):

    def is_nt(s):
        return pd.isna(s) or (s == 'None') or (s == 'DMSO') or (s=='')

    t_deets, c_deets = samp_deets.loc[treat], samp_deets.loc[ctrl]

    # get chem treat str
    if not is_nt(t_deets.Treatment):
        if is_nt(c_deets.Treatment):
            chem_str = t_deets.Treatment
        else:
            chem_str = f"{c_deets.Treatment}{ARROW}{t_deets.Treatment}"
    # the chem_str is only blank when there was non-genetic perturbation
    else:
        chem_str = ''

    if not t_deets.KO:
        ko_str = ''
    else:
        if not c_deets.KO:
            ko_str = t_deets.KO+'-KO'
        else:
            ko_str = f"{c_deets.KO}-KO{ARROW}{t_deets.KO}-KO"

    if chem_str and not ko_str:
        treatment_str = chem_str
    elif ko_str and not chem_str:
        treatment_str = ko_str
    else:
        if t_deets.Treatment == c_deets.Treatment:
            # it's a KO in a chemical background
            treatment_str = f"{ko_str} (with {chem_str})"
        else:
            treatment_str = f"{chem_str} (in {ko_str} cells)"
    return treatment_str


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


def set_logger(log_fn, level=logging.INFO):

    hndlr = logging.FileHandler(log_fn, 'w')
    hndlr.setLevel(level)

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
                     pseudocount=1,
                     drop_guide_less_than=0,
                     drop_guide_type=None):
    """output files written to {file_prefix}.{ctrl}-{treat}.tsv

    One file per comparison, in the default drugz format. Use tabulate_drugz to
    get a standardised table

    use drop_guide_* args remove guides that appear below given abundance.
    Valid drop_guide_type is 'any', 'all'. These can also be set via kwargs."""

    if kwargs:
        try:
            drop_guide_less_than = int(kwargs['drop_guide_less_than'])
            drop_guide_type = kwargs['drop_guide_type']
        except KeyError:
            pass

    if drop_guide_less_than:
        assert drop_guide_type in ('any', 'all', 'both')

    kwargs = kwargs
    if kwargs is None:
        kwargs = {}
    #todo remove attrdict usage
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

    # To drop guides we will be writing new temp count files using
    #   this DF.
    if drop_guide_less_than:
        tmpfn = f'{prefix}.tmp_hsfdok.tsv'
        dzargs.infile = tmpfn
        cnt = pd.read_csv(counts_file, index_col=0, sep='\t')


    # go through control and treatments to be analysed
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
            
            if file_exists_or_zero(dzargs.drugz_output_file):
                os.remove(dzargs.drugz_output_file)

            if drop_guide_less_than:
                reps = list_not_str(sample_reps[ctrl_samp]) \
                       + list_not_str(sample_reps[treat_samp])

                cnt_lessthan = (cnt.loc[:, reps] < drop_guide_less_than)
                # if add more types here, add new keywords to assertion at top of func
                if drop_guide_type == 'any':
                    bad_guides = cnt_lessthan.any(1)
                #elif drop_guide_type in ('both', 'all'):
                else:
                     bad_guides = cnt_lessthan.all(1)

                cnt.loc[~bad_guides, ['gene']+reps].to_csv(tmpfn, sep='\t')
                pipeLOG.info(
                    f"call_drugZ_batch: {sum(bad_guides)} guides removed "
                    f"from {ctrl_samp}-{treat_samp}, using less than {drop_guide_less_than} "
                    f"with '{drop_guide_type}' method."
                )

            drugZ_analysis(dzargs)

    if drop_guide_less_than:
        os.remove(tmpfn)
    pipeLOG.info('Finished drugZ')


def call_mageck(control_samp:str, treat_samp:str, sample_reps:Dict[str, List[str]],
                counts_file:str, prefix:str, kwargs:Dict,
                logger=None, dryrun=False, ):
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

    # Using a temporary file to deal with .gz and add pseudocount if set
    tmp_fn = prefix+'temp_counts.tsv'
    counts = pd.read_csv(counts_file, sep='\t', index_col=0)

    # Deal with pseudocount. Mageck doesn't have this facility, so we need
    #   to write a temporary counts file.
    # It adds 1 before logging anyway.
    if (pseudocount > 1):
        num_cols = counts.dtypes != object
        counts.loc[:, num_cols] += pseudocount

    counts.to_csv(tmp_fn, sep='\t')
    counts_file = tmp_fn

    for ctrl_samp, treat_samples in control_map.items():
        
        if type(treat_samples) is str:
            treat_samples = [treat_samples]

        for treat in treat_samples:
            if treat == ctrl_samp:
                continue
            
            outfn = '{outprefix}.{ctrlnm}-{sampnm}.gene_summary.txt'.format(outprefix=prefix, ctrlnm=ctrl_samp, sampnm=treat)
            if not file_exists_or_zero(outfn):

                call_mageck(ctrl_samp, treat, sample_reps, counts_file, prefix, kwargs)
                mageck_pairs_done.append((ctrl_samp, treat))
                
            else:
                pipeLOG.info('Output MAGeCK analysis exists: ' + outfn + ". Not running MAGeCK.")


    # delete the file with additional pseudocount, if there is one.
    if tmp_fn is not None:
        os.remove(tmp_fn)

# When adding new functions, need to deal with pseudocount option in the main func
#   and pairedness in the AnalysisWorkbook
def dry_function(*args, **kwargs):
    pass

analysis_functions = {'mageck':call_mageck_batch, 'drugz':call_drugZ_batch, }
analysis_tabulate = {'mageck':tabulate_mageck, 'drugz':tabulate_drugz, }
available_analyses = analysis_functions.keys()
analysis_functions['dry'] = dry_function
analysis_tabulate['dry'] = dry_function

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

#todo do we need to check for multiple controls anymore without jacks?
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
                         'analyses', 'file_prefix',]

    analyses = arguments['analyses']
    control_groups = arguments['control_groups']

    option_is_missing = lambda op: (op not in arguments.keys()) or (not arguments[op])
    missing_options = [op for op in required_options if option_is_missing(op)]
    if missing_options:
        raise RuntimeError(
            f'Required options missing (or valueless) from config file: {",".join(missing_options)}'
        )

    # check the required analyses options present
    required_analysis_opts = ['method']
    for k in required_analysis_opts:
        for ans in analyses:
            if k not in ans:
                raise RuntimeError(
                    f'Required analysis {k} option missing from analysis {ans}.'
                )

    # and the analyses dics themselves are formatted correctly.
    analysis_opts_types = {'method':str, 'kwargs':dict, 'groups':list, 'counts_file':str,
                           'label':str}
    for analysis_dict in analyses:
        for opt, tipe in analysis_opts_types.items():
            try:
                if type(analysis_dict[opt]) != tipe:
                    raise RuntimeError(
                        f"Analysis option types not as expected: {analyses}"
                    )
            # there are optional options, ignore them if they're missing
            except KeyError:
                pass

    samples = arguments['sample_reps'].keys()
    if any(['-' in s for s in samples]):
        raise RuntimeError(f'One or more sample names contain "-", which is not allowed.')
    # Find samples specified in comparisons but not found in sample_reps
    missing_samples = set()
    for _, ctrl_samps in control_groups.items():
        for ctrl, samps in ctrl_samps.items():
            if ctrl not in samples:
                missing_samples.add(ctrl)
            for smp in samps:
                if smp not in samples:
                    missing_samples.add(smp)
    missing_samples = list(missing_samples)
    if missing_samples:
        raise RuntimeError(f"Samples named in control_groups missing in sample_reps:\n{', '.join(missing_samples)}")

    missing_groups = set()
    for ans in analyses:
        for grp in ans['groups']:
            if grp not in control_groups.keys():
                missing_groups.add(grp)
        if missing_groups:
            raise RuntimeError(
                f'Group(s) in Analyses not found in Control Groups: {missing_groups}.\n'
                f'Control Groups groups {list(control_groups.keys())}'
            )



def process_arguments(arguments:dict, delete_unrequired_args=True) -> dict:
    """Deal with special keywords from the experiment dictionary.
    Varefy integrity of options.

    If delete_unrequired_args, arguments that come from AnalysisWorkbook
    but not required by run_analyses are removed."""
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
        if cfn.endswith('.gz'):
            open_count_file = lambda fn: gzip.open(fn, 'rt')
        else:
            open_count_file = open
        with open_count_file(cfn) as f:
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

    # Remove control groups not specified in run_groups, if run_groups set.
    if 'run_groups' in arguments and arguments['run_groups'] is not None:
        run_groups = arguments['run_groups']
        if type(run_groups) is str:
            run_groups = run_groups.split(',')
        remove_groups = set()
        for grp in arguments['control_groups']:
            if grp not in run_groups:
                remove_groups.add(grp)
        pipeLOG.info(f'Removing control groups: {", ".join(sorted(remove_groups))}')
        for grp in remove_groups:
            del arguments['control_groups'][grp]
        if len(arguments["control_groups"]) == 0:
            raise PipelineOptionsError(f'No control groups left with run_groups option: {run_groups}')
        pipeLOG.info(f'Running pipeline with remaining groups: {",".join(arguments["control_groups"].keys())}')
    # remove before calling run_analyses
    if 'run_groups' in arguments:
        del arguments['run_groups']

    # remove analyses not specified in run_analysis if any set
    if 'run_analyses' in arguments and arguments['run_analyses'] is not None:
        analyses_to_run = arguments['run_analyses']
        if type(analyses_to_run) is str:
            analyses_to_run = analyses_to_run.split(',')
        filtered_analyses = []
        removed_analyses = set()
        for analysis in arguments['analyses']:
            if analysis['name'] in analyses_to_run:
                filtered_analyses.append(analysis)
            else:
                removed_analyses.add(analysis['name'])
        arguments['analyses'] = filtered_analyses
        pipeLOG.info(f"Removed analyses named: {', '.join(sorted(removed_analyses))}")
        if len(arguments['analyses']) == 0:
            raise PipelineOptionsError(f'No analyses left with run_analyses option: {analyses_to_run}')
    if 'run_analyses' in arguments:
        del arguments['run_analyses']


    # These don't need to be passed to the pipeline.
    if delete_unrequired_args:
        for k in ['config_file', 'notes', 'experiment_id', 'analysis_version',
                  'labels']:
            if k in arguments:
                del arguments[k]

    return arguments


def run_analyses(output_dir, file_prefix,
                 sample_reps:Dict[str, list],
                 control_groups:Dict[str, Dict[str, list]],
                 analyses:List[dict],
                 counts_file=None,
                 methods_kwargs:Dict=None,
                 dont_log=False,  #todo replace dont_log with some kind of verbosity thing
                 compjoiner='-',  #tools.ARROW,
                 notes='',
                 skip_method=None,
                 dry_run=False,
                 **kwarghole,
                 ):

    """Run batches of CRISPR analyses.

    Args:
        output_dir: Directory within which all files will be written.
        file_prefix: Prefix for all written files (should not contain directories)
        sample_reps: Dictionary defining sample_name to list of replicate names.
        control_groups: Groups of maps; which samples will be used in comparisons.
            {ctrl_grp_name: {ctrl_sample: [test_samp1, ...]}}
        analyses: List of analyses to be performed. Each dict contains all required
            options.
        counts_file: counts files specified per analysis take precident, this value
            used only when one is not specified in an analysis.
        methods_kwargs: Options to be passed to analysis methods. {analysis_name: {}}
        dont_log: Set to True to supress writing logs. Will be preplaced by some
            kind of verbosity setting.
        compjoiner: In written files, character(s) used to join ctrl_samp & treat_samp.
    """

    if kwarghole:
        pipeLOG.warning(f"Unexpected keyword arguments: {kwarghole}")

    #pipeLOG.warning(f'Unrecognised arguments passed to run_analyses:\n   {unrecognised_kwargs}')
    pipeLOG.info(f"notes: {notes}")

    if skip_method is None:
        skip_method = []

    # Create the root experiment directory
    p = str(output_dir)
    os.makedirs(p, exist_ok=True)

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

    def mkdir(p):
        """os.mkdir but it ignores FileExistsError"""
        try:
            os.mkdir(p)
        except FileExistsError:
            pass
    print(output_dir)
    mkdir(output_dir)
    mkdir(str(Path(output_dir, 'tables')))
    for analysis_method in methods_used:

        mkdir(str(Path(output_dir, analysis_method)))
        mkdir(str(Path(output_dir, analysis_method, 'files')))

    ######################
    ## Run the analyses ##
    ######################
    ran_analyses = set()

    for analysis_dict in analyses:

        if dry_run:
            analysis_method = 'dry'
        else:
            analysis_method = analysis_dict['method']

        if analysis_method in skip_method:
            pipeLOG.info(f"Skipping analysis method {analysis_method}")
            continue

        # The structure here (which may be pointlessly complicated) is that
        #   default arguments for an analysis method can be set at the top
        #   level of the JSON, but individual analysis_dicts may contain
        #   arguments that override these.
        try:
            default_method_kwargs = methods_kwargs[analysis_method]
            if not default_method_kwargs:
                default_method_kwargs = {}
        except:
            default_method_kwargs = {}

        # if we do not specify groups, use all control groups.
        if 'groups' not in analysis_dict:
            groups = {g:{} for g in control_groups.keys()}
        else:
            groups = list_not_str(analysis_dict['groups'])

        try:
            if analysis_dict['kwargs']:
                curr_kwargs = analysis_dict['kwargs']
            else:
                curr_kwargs = default_method_kwargs
        except KeyError:
            curr_kwargs = default_method_kwargs

        # go through the selected groups for this analysis and run the thing
        for grp in groups:

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
            if not dry_run:
                tab.to_csv(tabfn, encoding='utf-8-sig')


def load_configuration_file(config_filename, counts_dir='.') -> dict:
    """
    Args:
        config_filename: An .xlsx or .json containing information to
            run analyses.
        counts_dir: A directory containing the counts files, if not
            the current working dir.
    """

    # Load the configuration file
    config_file = config_filename

    if config_file.endswith('xlsx'):
            config_file_args = data_classes.AnalysisWorkbook(config_file, counts_dir=counts_dir).expd

    # Stripping out the "comment" lines
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

    return config_file_args

if __name__ == '__main__':

    pipeLOG.info(str(__version__))
    parser = argparse.ArgumentParser(
        description="Run mageck and jacks analyses using a JSON file.\n "
                    "Use arguments below to override JSON options"
    )

    parser.add_argument(
        'config_file', metavar='FILE',
        help = ("path to .xlsx or .json file specifying arguments. JSON must must contain `sample_reps`"
                ", `analyses` & `control_groups` keywords and values. "
                "Other options may be overridden by command line arguments.")
        )
    parser.add_argument('--counts', metavar='FILE/DIR',
                        help='Path to counts file, or directory containing counts file(s)'
                             ' specified in config_file',
                        default=None, dest='counts_file')
    parser.add_argument('--output-dir', metavar='PATH', default='.',
                        help=('Path to where results directory will be created if not current'
                              ' directory. Experiment_id and analysis_version will determine results dir.') )
    parser.add_argument('--file-prefix', default='result',
                        help="String to form identifying prefix for all files generated.")
    parser.add_argument('--skip-method', metavar='list,of,progs', default=None,
                        help='Filter analyses methods to be run.')
    parser.add_argument('--run-groups', metavar='list,of,groups', default=None,
                        help='Specify control groups (as defined in exp dict) to be included. All included by default.')
    parser.add_argument('--run-analyses', metavar='list,of,names', default=None,
                        help='Specify analyses by name to run. Names defined on Analyses sheet or {"name":name} in'
                             ' analyses dictionaries. All included by default.')
    parser.add_argument('--dont-log', action='store_true', dest='dont_log', default=None,
                        help="Don't write a log file.")
    parser.add_argument('--analysis-version', default=None,
                        help='Output files will be stored in a directory of the above name, within the experiment dir.')
    parser.add_argument('--dry-run', action='store_true', default=False,
                        help='Tries to validate all the options specified without actually running anything. ' 
                             'Not guaranteed to spot all issues.')

    ## any new args will be passed to run_analyses by default, so delete em before
    # parser.add_argument('--debug', action='store_true', default=False,
    #                     help="Set pipeLOG level to debug")

    # get arguments from the command line and the YAML
    clargs = parser.parse_args() # need to assign this before calling vars() for some reason
    cmd_line_args = vars(clargs)

    # if cmd_line_args['debug']:
    #     pipeLOG.setLevel(logging.DEBUG)
    #     pipeLOG.debug('Debugging level is set')
    # del cmd_line_args['debug']

    # if the command line arg gives a directory for counts, parse that
    counts_dir = ''
    if 'counts_file' in cmd_line_args and cmd_line_args['counts_file']:
        cntpath = cmd_line_args['counts_file']
        if os.path.isdir(cntpath):
            counts_dir = cntpath
            del cmd_line_args['counts_file']
        elif not os.path.isfile(cntpath):
            del cmd_line_args['counts_file']

    config_file_args = load_configuration_file(cmd_line_args['config_file'], counts_dir)
    #print(config_file_args)
    # over write config file args with any specified in the command line
    for k, v in cmd_line_args.items():
        if v is not None:
            config_file_args[k] = v

    # validate and process arguments
    args = process_arguments(config_file_args)

    run_analyses(**args)

