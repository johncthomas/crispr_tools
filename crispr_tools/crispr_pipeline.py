#!/usr/bin/env python3

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

try:
    from jacks.jacks_io import runJACKS
except ImportError:
    print('To run jacks you need to install JACKS,\n', 'https://github.com/felicityallen/JACKS/tree/master/jacks\n' 
          "You can still run Mageck though, if it's installed.")
    def runJACKS(*a, **k):
        raise ModuleNotFoundError('JACKS not installed!!!')

from crispr_tools.drugz import drugZ_analysis

from crispr_tools.tools import list_not_str

# with open(pathlib.Path(__file__).parent/'version.txt') as f:
#     __version__ = f.readline().replace('\n', '')

pipeLOG = logging.getLogger('pipeline')
pipeLOG.setLevel(logging.INFO)

from jacks.infer import LOG as jacksLOG
jacksLOG.setLevel(logging.WARNING)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

#from crispr_tools.count_reads import count_reads, count_batch, map_counts
from crispr_tools.tools import plot_read_violins, plot_ROC, plot_volcano, tabulate_mageck, plot_volcano_from_mageck
from crispr_tools.jacks_tools import tabulate_score, scores_scatterplot, mahal_nocov
from crispr_tools.version import __version__


#todo: refactor the analysis so it's modulor and not just a massive function
#todo check validity of every file and option before starting.
#  eg:
#   do all reps exsist in the count file
#   are samples in controls defined in the sample reps
#todo QC by default!!
#todo save yaml in output dir
#todo handle different analyses better, making it easy to add new ones and turn of during analysis


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


def call_JACKS(fn_counts, fn_repmap, outprefix, norm='mode'):
    # it'll probably be easiest to just write a temp repmap tsv from the xlsx
    xl = pd.ExcelFile(fn_repmap)
    for ctrlgroup in xl.sheet_names:
        repmap = xl.parse(ctrlgroup)
        repmap.to_csv('tmp.repmap.tsv', '\t')
        runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                 norm_type=norm, outprefix=outprefix+'.'+ctrlgroup+'.')

def set_logger(log_fn):

    hndlr = logging.FileHandler(log_fn, 'w')
    # hndlr.setLevel(logging.INFO)
    pipeLOG.setLevel(logging.INFO)
    pipeLOG.addHandler(hndlr)
    jacksLOG.addHandler(hndlr)
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


# class PseudoArgs:
#     """for holding attributes"""
#     # https://stackoverflow.com/questions/4984647/accessing-dict-keys-like-an-attribute
#     def __init__(self, *args, **kwargs):
#         super(PseudoArgs, self).__init__(*args, **kwargs)
#         self.__dict__ = self





def run_drugZ(fn_counts, outdir, file_prefix,
                 sample_reps:Dict[str, list],
                 controls:Dict[str, Dict[str, list]], drugz_args:dict=None):

    dzargs = AttrDict()
    dzargs.infile = fn_counts
    # defaults
    dzargs.index_column = 0
    dzargs.minobs = 1
    dzargs.half_window_size = 500
    dzargs.quiet = False
    dzargs.pseudocount = 1
    dzargs.fc_outfile=''

    # not implimented in drugZ at time of writing
    dzargs.remove_genes = None

    dzargs.update(drugz_args)

    call(['mkdir', '-p', str(Path(outdir, 'drugz', 'files'))])

    pipeLOG.info('Running drugZ')

    for control_group, control_map in controls.items():
        for ctrl_samp, treat_samples in control_map.items():

            dzargs.control_samples = ','.join(
                list_not_str( sample_reps[ctrl_samp] )
            )

            treat_samples = list_not_str(treat_samples)
            for treat_samp in treat_samples:

                dzargs.drug_samples = ','.join(
                    list_not_str(sample_reps[treat_samp])
                )

                dzargs.drugz_output_file = f"{outdir}/drugz/files/{file_prefix}.{ctrl_samp}-{treat_samp}.tsv".replace('//', '/')

                drugZ_analysis(dzargs)



    pipeLOG.info('Finished drugZ')


def call_mageck(control_samp:str, treat_samp:str, sample_reps:Dict[str, List[str]],
                fn_counts:str, prefix:str, kwargs:dict):
    mageck_str = "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}{ctrlnm}-{sampnm}"
    # get the replicate strings
    ctrls = ','.join(sample_reps[control_samp])
    treats = ','.join(sample_reps[treat_samp])
    s = mageck_str.format(
        counts=fn_counts,
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
    pipeLOG.info(' '.join(mag_args))
    call(' '.join(mag_args), shell=True)

def run_mageck_batch(sample_reps:Dict[str, list],
                     control_map:Dict[str, list],
                     count_fn:str,
                     prefix:str,
                     mag_kwargs:dict=None,
                     skip_extra_mageck=True):

    """Run mageck analyses using comparisons specified in control_map."""

    mageck_pairs_done = [] #todo remove this, or do something with it
    call("which mageck".split())
    pipeLOG.info('Running MAGeCK version ' + check_output(["mageck", "-v"]).decode())

    for ctrl_samp, treat_samples in control_map.items():

        if type(treat_samples) is str:
            treat_samples = [treat_samples]

        for treat in treat_samples:
            if treat == ctrl_samp:
                continue

            call_mageck(ctrl_samp, treat, sample_reps, count_fn, prefix, mag_kwargs)
            mageck_pairs_done.append((ctrl_samp, treat))
        if not skip_extra_mageck:
            # run all combinations of samples that share at least one control sample
            # replace the control group name in the string
            pref_split = prefix.split('.')
            pref_split[-2] = 'EXTRA'
            prefix = '.'.join(pref_split)
            for c, t in combinations(treat_samples, 2):
                if (c, t) not in mageck_pairs_done and (t, c) not in mageck_pairs_done:
                    mageck_pairs_done.append((c, t))
                    # maybe don't tablulate these at the moment.
                    call_mageck(c, t, sample_reps, mag_kwargs)


def run_analysis(fn_counts, outdir, file_prefix,
                 sample_reps:Dict[str, list],
                 controls:Dict[str, Dict[str, list]],
                 comparisons: List[dict],
                 volc_labels=15, scatter_labels=10,
                 charts_only = False, jacks_eff_fn=None,
                 skip_mageck = False, skip_jacks = False, skip_drugz=True,
                 skip_charts=False,
                 dont_log=False, exp_name='', analysis_name='',
                 skip_extra_mageck = False, jacks_kwargs:Dict=None,
                 mageck_kwargs:dict=None,
                 ctrl_genes:str=None, notes='', **unused_args):


    #call(['mkdir', outdir])

        #haaaacky, but the fact it cant make a dir and subdir at the same time is annoying
    p = str(outdir)
    for i in range(len(p)):
        print(p)
        d = p.split('/')[:i+1]
        try:
            os.mkdir('/'.join(d))
        except FileExistsError:
            pass

    if jacks_kwargs is None:
        jacks_kwargs = {}

    if ctrl_genes:
        jacks_kwargs['ctrl_genes'] = ctrl_genes
        #todo write a temp file for mageck, and clean up all temp files (put paths in a list)

    for analysis_str in ('jacks_mode', 'jacks_median', 'mageck'):
        call(['mkdir', '-p', str(Path(outdir, analysis_str))])
        call(['mkdir', str(Path(outdir, analysis_str, 'volcano'))])
        call(['mkdir', str(Path(outdir, analysis_str, 'scatter'))])
        call(['mkdir', str(Path(outdir, analysis_str, 'tables'))])
        call(['mkdir', str(Path(outdir, analysis_str, 'files'))])
        #call(['mkdir', str(Path(outdir, 'scattercharts'))])

    if not dont_log:
        t = '{}-{}-{}_{}h{}m{}s.{}'.format(*datetime.datetime.now().timetuple()[:-1])
        set_logger(str(Path(outdir, file_prefix + f'log_{t}.txt')))

        sys.stderr = StreamToLogger(pipeLOG, logging.ERROR)
    pipeLOG.info('exp_name = ' + exp_name)
    pipeLOG.info('analysis_name = ' + analysis_name)
    pipeLOG.info('Version = '+__version__)
    pipeLOG.info('Working dir = '+os.getcwd())
    pipeLOG.info('outdir = '+outdir)

    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    pipeLOG.info('\n\t'.join([f"{x} = {values[x]}" for x in args]))
    pipeLOG.info(f"Full outdir = {os.path.realpath(outdir)}")

    pipeLOG.info('Control groups being used: '+', '.join(controls.keys()))

    for ctrlgroup, ctrlmap in controls.items():
        # File name prefixes
        fnprefix = file_prefix + '.' + ctrlgroup + '.'
        pipeLOG.info(f'\nRunning {ctrlgroup}, writing to {fnprefix}')
        jmode_prefix = str(Path(outdir, 'jacks_mode', 'files', fnprefix))
        jmed_prefix = str(Path(outdir, 'jacks_median', 'files', fnprefix))
        mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))

        #################
        # Run JACKS
        if not charts_only and not skip_jacks:
            # make the jacks repmap: rep, samp, ctrl
            repmap_fn = outdir + ctrlgroup + '.repmap.tsv'
            write_repmap(sample_reps, ctrlmap, repmap_fn)
            pipeLOG.info(f"Running JACKS, mode and median normalised with {repmap_fn}")
            if jacks_eff_fn:
                jacks_eff = jacks_eff_fn.format(ctrlgroup)
                if os.path.isfile(jacks_eff):
                    pipeLOG.info(f'\tUsing {jacks_eff} for efficacy')
                else:
                    pipeLOG.info(f"\tcould not find {jacks_eff}")
                    raise FileNotFoundError
            else:
                jacks_eff = None

            jacks_args = (fn_counts, repmap_fn, fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',)

            runJACKS(*jacks_args, norm_type='mode', outprefix=jmode_prefix, reffile=jacks_eff, **jacks_kwargs)
            runJACKS(*jacks_args, norm_type='median', outprefix=jmed_prefix, reffile=jacks_eff, **jacks_kwargs)

        #########
        # *run MAGECK
        if not charts_only and not skip_mageck:
            run_mageck_batch(sample_reps, controls[ctrlgroup], fn_counts, mf_prefix, mageck_kwargs, skip_extra_mageck)
            # mageck_pairs_done = []
            # def _run_mageck(_ctrl_samp, _treat, prefix = mf_prefix):
            #     mageck_str = "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}{ctrlnm}-{sampnm}"
            #     # get the replicate strings
            #     ctrls = ','.join(sample_reps[_ctrl_samp])
            #     treats = ','.join(sample_reps[_treat])
            #     s = mageck_str.format(
            #         counts=fn_counts,
            #         treat=treats,
            #         ctrl=ctrls,
            #         outprefix=prefix,
            #         ctrlnm=_ctrl_samp,
            #         sampnm=_treat
            #     )
            #     mag_additional_args = []
            #     if mageck_kwargs:
            #         mag_additional_args = [f"--{_k} {_v}" for _k, _v in mageck_kwargs.items()]
            #     mag_args = s.split()+mag_additional_args
            #     # for some reason mageck fails to understand mag_args if you don't use shell=True
            #     pipeLOG.info(' '.join(mag_args))
            #     call(' '.join(mag_args), shell=True)
            #
            # call("which mageck".split())
            # pipeLOG.info('Running MAGeCK version '+check_output(["mageck","-v"]).decode())
            #
            # for ctrl_samp, treat_samps in ctrlmap.items():
            #     if type(treat_samps) is str:
            #         treat_samps = [treat_samps]
            #     for treat in treat_samps:
            #         if treat == ctrl_samp:
            #             continue
            #         # ctrls = ','.join(sample_reps[_ctrl_samp])
            #         # treats = ','.join(sample_reps[_treat])
            #         _run_mageck(ctrl_samp, treat)
            #         mageck_pairs_done.append((ctrl_samp, treat))
            #     if not skip_extra_mageck:
            #         # run all combinations of samples that share at least one control sample
            #         for c, t in combinations(treat_samps, 2):
            #             if (c,t) not in mageck_pairs_done and (t,c) not in mageck_pairs_done:
            #                 mageck_pairs_done.append((c,t))
            #                 # maybe don't tablulate these at the moment.
            #                 _run_mageck(c, t, prefix=mf_prefix.replace(ctrlgroup, 'EXTRA'))


        if not charts_only and not skip_drugz:
            pipeLOG.info('Running DrugZ.')
            run_drugZ(fn_counts, outdir, file_prefix, sample_reps, controls, {})

        ############
        #* Get tables
        #todo probably have a dict of tables that is populated as we go. And probably just use snakemake
        analyses_used = []
        if not skip_jacks:
            scoresmode = tabulate_score(jmode_prefix)
            scoresmed = tabulate_score(jmed_prefix)
            analyses_used.extend([('jacks_mode', 'jacks_score',  scoresmode),
                                  ('jacks_median', 'jacks_score', scoresmed)])
            if not charts_only:
                pipeLOG.info('Writing JACKS tables')
                modfn = str(Path(outdir, 'jacks_mode', 'tables', fnprefix+'jacks_scores.csv'))
                medfn = str(Path(outdir, 'jacks_median', 'tables', fnprefix + 'jacks_scores.csv'))
                pipeLOG.info('Writing: '+modfn+'\n\t& '+medfn)
                scoresmode.to_csv(modfn)
                scoresmed.to_csv(medfn)

        if not skip_mageck:
            scores_mageck = tabulate_mageck(mf_prefix)
            analyses_used.append(('mageck', 'lfc', scores_mageck))
            if not charts_only:
                magfn = str(Path(outdir, 'mageck', 'tables', fnprefix + 'mageck_table.csv'))
                pipeLOG.info('Writing MAGeCK tables '+magfn)
                scores_mageck.to_csv(magfn)

        #################
        #* volcano charts of all
        # use the universal one, have label options as part of args, maybe as a dict.
        if not skip_charts:
            for analysis_str, xkey, analysis_tab in analyses_used:
                pipeLOG.info('Volcano plot, '+analysis_str)
                #analysis = analysis
                for exp in analysis_tab.columns.levels[0]:
                    if analysis_str.startswith('jacks'):
                        chart_title = '{} from {} ({} normalisation)'.format(exp, ctrlgroup, analysis_str)
                    else:
                        # split by fnpref just in case it has a dash in it...
                        #print(exp)
                        ctrlnm, treatnm = exp.split('-')
                        chart_title = '{} from {} ({} analysis)'.format(treatnm, ctrlnm, analysis_str)
                    plot_volcano(
                        xkey, 'fdr_log10', analysis_tab[exp], chart_title,
                        volc_labels,volc_labels,p_thresh=0.1,
                        outfn=str(Path(outdir, analysis_str, 'volcano', fnprefix+exp+'.'+analysis_str+'_volcano.png'))
                    )
                    plt.close()

            # call(['tar', '-zcf',
            #       str(Path(outdir, analysis_str+'volcano_charts.tar.gz')), #file to be created
            #       # change to the target dir so it doesn't include the path in the tar
            #       '-C', str(Path(outdir, analysis_str, 'volcano')), '.'])


        analyses_tups = []
        if not skip_jacks:
            analyses_tups.extend([('jacks_mode', scoresmode), ('jacks_median', scoresmed)])
        # if not skip_mageck:
        #     analyses_tups.append(('mageck', scores_mageck))

        # do jacks scatter and mahal table. tups should never be empty
        for analysis_str, analysis_tab in analyses_tups:
            pipeLOG.info(f'Doing comparisons of {analysis_str}')
            mahals = pd.DataFrame(index=analysis_tab.index)
            comparisons_done = False
            for comp_from, comp_to in iter_comps(comparisons, analysis_tab):
                comparisons_done = True
                pipeLOG.info(f"\t{comp_from} vs {comp_to}")
                #scatter

                if not skip_charts:
                    scores_scatterplot(comp_from, comp_to, analysis_tab,
                                       label_pos=scatter_labels, label_neg=scatter_labels,
                                       dist_name='Mahalanobis distance')
                    plt.title(f"{comp_from} vs {comp_to} ({ctrlgroup}, JACKS)")
                    fn = str(
                        Path(outdir, analysis_str, 'scatter', fnprefix+"{}_vs_{}.scatter.png".format(comp_from, comp_to))
                    )
                    pipeLOG.info('Writing: '+fn)
                    plt.savefig(fn, dpi=150)
                    plt.close()

                # pop mahal table
                A, B = analysis_tab[comp_from], analysis_tab[comp_to]
                _, _, mahal = mahal_nocov(A['jacks_score'], A['stdev'], B['jacks_score'], B['stdev'])
                mahals.loc[:, f"{comp_from} vs {comp_to}"] = mahal

            if comparisons_done:
                mahalfn = Path(outdir, analysis_str, 'tables', fnprefix+f"mahalanobis_distances.csv")
                pipeLOG.info('Writing: '+str(mahalfn))
                mahals.to_csv(mahalfn)

            # # wipes over the one produced above if comptab exists.
            # call(['tar', '-zcf', str(Path(outdir, analysis_str + 'scatter_charts.tar.gz')),
            #       '-C', str(Path(outdir, analysis_str, 'scatter')), '.'])
            #         #str(Path(outdir, analysis_str, 'scatter'))])
            # call(['tar', '-zcf', str(Path(outdir, analysis_str + 'tables.tar.gz')),
            #       '-C', str(Path(outdir, analysis_str, 'tables')), '.'])
            # #str(Path(outdir, analysis_str, 'tables'))])

    #* END OF control group loop
    ###################

    # Do all the mageck comparisons at the end as we need the NT->treat for significance
    if not skip_mageck:
        pipeLOG.info('Doing comparisons of mageck')
        mag_tables = {}
        # Get the tables, including the intersample combinations
        ctrlgroups = list(controls.keys())
        if not skip_extra_mageck:
            ctrlgroups += ['EXTRA']
        for ctrlgroup in ctrlgroups:
            fnprefix = file_prefix + '.' + ctrlgroup + '.'
            mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))
            scores_mageck = tabulate_mageck(mf_prefix)
            mag_tables[ctrlgroup] = scores_mageck
        # go through each comparison pair and check for the existence of the right measures.
        for comp_from, comp_to in iter_comps(comparisons):
            lfc_tabs = []
            fdr_tabs = []
            for ctrlgroup, tab in mag_tables.items():
                # column headers are "ctrl-samp"
                exp = pd.Series(tab.columns.levels[0], index=tab.columns.levels[0])
                samps = set([c.split('-')[1] for c in exp])
                # so, if a comp is specified with NT->T then we want to use NT->T for fdr and
                # ctrl->[NT/T] LFCs for the scatter, then name by ctrl.NT_T
                if comp_from in samps and comp_to in samps:
                    treat_header = exp[exp.str.contains(comp_to)][0]
                    ctrl_header = exp[exp.str.contains(comp_from)][0]
                    # some weird groups will exist as "ctrlgroup", so check the ctrls actually match
                    if ctrl_header.split('-')[0] == treat_header.split('-')[0]:
                        lfc_tabs.append((ctrlgroup, ctrl_header,  treat_header,))
                # if f'{comp_to}-{comp_from}' in exp:
                #     fdr_tabs.append((ctrlgroup, f'{comp_to}-{comp_from}'))
                if f'{comp_from}-{comp_to}' in exp:
                    fdr_tabs.append((ctrlgroup, f'{comp_from}-{comp_to}'))
                else: # no source of sig
                    pipeLOG.info(f'{comp_from}->{comp_to} not possible, missing comparison.')
                    continue
            # print(lfc_tabs)
            # print('**', fdr_tabs)
            if not skip_charts:
                for fdrtab, fdrexp in fdr_tabs:
                    for ctrlgroup, ctrl_header,  treat_header in lfc_tabs:
                        pipeLOG.info(f'MAgeck results comp using {ctrl_header} {treat_header}, {fdrexp}')

                        scores_scatterplot(ctrl_header, treat_header, mag_tables[ctrlgroup], True, scatter_labels, scatter_labels,
                                           distance=mag_tables[fdrtab].loc[:, (fdrexp, 'fdr_log10')],
                                           min_label_dist=0.3,
                                           dist_name='log10(FDR)')
                        plt.title(f"{ctrl_header} vs {treat_header} (MAGeCK)")
                        fn = str(
                            Path(outdir, 'mageck', 'scatter', file_prefix+ '.' + ctrlgroup + ".{}_vs_{}.scatter.png".format(comp_from, comp_to))
                        )
                        pipeLOG.info('Writing: '+fn)
                        plt.savefig(fn, dpi=150)
                        plt.close()

    import shutil
    shutil.copy(fn_counts, Path(outdir, Path(fn_counts).name))

    pipeLOG.info("Done. "+str(datetime.datetime.now()))


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
        #todo: if ctrlmap is a list assume a single control group and don't crash
        #print('ctrlmap',  ctrlmap)
        for ctrl, samps in ctrlmap.items():
            if type(samps) == str:
                samps = [samps]
            # if you add keywords you'll need to add them to the database add_data module
            if samps == ['ALL']:
                samps = copy(list(samples))
            elif type(samps) == dict:
                # maybe could have other keywords in the future
                samps = [s for s in samples if s not in samps['EXCEPT']]
            ctrlmap[ctrl] = samps
    return controls

def process_arguments(arguments:dict):
    """deal with special keywords from the experiment yaml, and allow some
    ambiguous syntax in the yaml. Also do some checking of validity"""
    samples = arguments['sample_reps'].keys()
    controls = process_control_map(arguments['controls'], samples)
    arguments['controls'] = controls

    # in case strings are passed instead of lists
    if 'comparisons' in arguments:
        for comp in arguments['comparisons']:
            print(arguments['comparisons'])
            for k, v in comp.items():
                if type(v) == str:
                    comp[k] = [v]

    if 'sample_reps' in arguments:
        for k, v in arguments['sample_reps'].items():
            if type(v) == str:
                arguments['sample_reps'][k] = [v]

    with open(arguments['fn_counts']) as f:
        line = next(f)
        if not '\t' in line:
            raise ValueError('No tabs in file, is it comma seperated?\n\t'+line)

    return arguments


def do_qcs(count, expd, jobs = ('pca', 'clustermap', 'violinplots', 'regressions')):
    """Perform PCAs, clustermap, regressions between clones using both
    normalised counts and lfcs."""

    lncounts = tools.size_factor_normalise(count)
    samp_reps = expd['sample_reps']
    ctrl_map  = expd['controls']

    plot_read_violins(lncounts, size_norm=False, log=False)


if __name__ == '__main__':
    print(__version__)
    parser = argparse.ArgumentParser(
        description="Run mageck and jacks analyses using a YAML file.\n "
                    "Use arguments below to override yaml options"
    )
    parser.add_argument('fn_yaml', metavar="YAML",
                        help="path to .yml file specifying arguments. At a minimum must contain `sample_reps`"
                             " & `controls` keywords and values. "
                             "Other options are overridden by command line arguments.")
    parser.add_argument('--fn_counts', metavar='COUNTS', help='Path to counts file', default=None)
    parser.add_argument('--outdir', metavar='OUTDIR', default=None,
                        help='Path to where results files will be stored, a '
                        "directory structure will be created.")
    parser.add_argument('--file_prefix', metavar='PREFIX', default=None,
                        help="String to form identifying prefix for all files generated.")
    parser.add_argument('-v', '--volc-labels',dest='volc_labels', metavar='N', type=int, default=None,
                        help="Number of depleted/enriched genes to label in volcano charts.")
    parser.add_argument('-s', '--scatter-labels', metavar='N', type=int, default=None, dest='scatter_labels',
                        help="Number of genes to label on either side of The Line in scatter charts.")
    parser.add_argument('--charts-only', action='store_true', default=None,
                        help="Don't run MAGeCK or JACKS, just produce plots with existing files.")
    parser.add_argument('-j', '--jacks-efficacy', metavar='EFF_FILE', dest='jacks_eff_fn', default=None,
                        help='Path to efficacy file to be used in JACKS analysis. '
                             'If the provided string contains {}, the repmap ctrl group will be substituted'
                             'into the filename. Optional.')
    parser.add_argument('--skip-jacks', action='store_true', dest='skip_jacks', default=None,
                        help="don't run JACKS analysis or try to plot from JACKS analyses")

    parser.add_argument('--skip-mageck', action='store_true', dest='skip_mageck', default=None,
                        help="don't run MAGeCK analysis or try to plot from MAGeCK analyses")

    parser.add_argument('--skip-charts', action='store_true', dest='skip_charts', default=None,
                        help="don't produce any charts")

    parser.add_argument('--dont-log', action='store_true', dest='dont_log', default=None,
                        help="Don't write a log file.")



    # get arguments from the command line and the YAML
    clargs = parser.parse_args() # need to assign this before calling vars() for some reason
    cmd_args = vars(clargs)
    yml_args = yaml.safe_load(open(cmd_args['fn_yaml']))
    del cmd_args['fn_yaml']

    # generate outdir from exp_name and analysis_name
    if "outdir" not in yml_args:

        yml_args["outdir"] = os.path.join(yml_args["exp_name"],  yml_args["analysis_name"])


    # over write yml_args with any specified in the command line
    for k, v in cmd_args.items():
        if v is not None:
            yml_args[k] = v
    print(yml_args)
    args = process_arguments(yml_args)
    del args['labels']
    #print(args)

    run_analysis(**args)



    # run_analysis(**vars(clargs))
    # os.chdir('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/matylda_HS715-20')
    # run_analysis(**{'fn_counts': '/Users/johnc.thomas/Dropbox/crispr/counts_all/matylda_HS715-20.counts.tsv',
    #  'fn_repmap': 'matylda_HS715-20.repmap.3.xlsx', 'outdir': 'take5', 'file_prefix': 'HS715-20', 'labeldep': 30,
    #  'labelenr': 10, 'charts_only': False})