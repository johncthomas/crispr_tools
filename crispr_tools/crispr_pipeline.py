#!python3
import logging, os, pathlib, datetime, inspect, argparse, sys

from typing import List, Dict
from copy import copy
from subprocess import call, check_output
from pathlib import Path
from itertools import combinations

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    from jacks.jacks_io import runJACKS
except ImportError:
    print('To run jacks you need to install JACKS,\n', 'https://github.com/felicityallen/JACKS/tree/master/jacks\n' 
          "You can still run Mageck though, if it's installed.")
    def runJACKS(*a, **k):
        raise ModuleNotFoundError('JACKS not installed!!!')

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

#todo add sample labels to sample_details for nice chart labels
#todo? put each chunk of analysis into its own function

"""Go from FastQ files to completed JACKS/MAGeCK analysis. 
fq->mapped_counts&violin plots are one command (count_reads.count_batch() ) 
counts->charts&tables another."""

class StreamToLogger(object):
   """
   Fake file-like stream object that redirects writes to a logger instance.
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


def run_analysis(fn_counts, outdir, file_prefix,
                 sample_reps:Dict[str, list],
                 controls:Dict[str, Dict[str, list]],
                 comparisons: List[dict],
                 volc_labels=15, scatter_labels=10,
                 charts_only = False, jacks_eff_fn=None,
                 skip_mageck = False, skip_jacks = False, dont_log=False,
                 skip_extra_mageck = False, notes=''):


    call(['mkdir', outdir])


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
            reptab = [['rep', 'samp', 'ctrl']]
            for ctrl_samp, treat_samps in ctrlmap.items():
                for ctrl_rep in sample_reps[ctrl_samp]:
                    reptab.append([ctrl_rep, ctrl_samp, ctrl_samp])
                for treat in treat_samps:
                    for rep in sample_reps[treat]:
                        reptab.append([rep, treat, ctrl_samp])
                        # print([rep, samp, ctrl])
                        del rep
                    del treat
            repmap_fn = outdir + ctrlgroup + '.repmap.tsv'
            with open(repmap_fn, 'w') as f:
                for line in reptab:
                    f.write('\t'.join(line) + '\n')


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
            #jacks_kwargs =  dict()
            runJACKS(*jacks_args, norm_type='mode', outprefix=jmode_prefix, reffile=jacks_eff)
            runJACKS(*jacks_args, norm_type='median', outprefix=jmed_prefix, reffile=jacks_eff)

        #########
        # *run MAGECK
        if not charts_only and not skip_mageck:
            mageck_pairs_done = []
            def _run_mageck(_ctrl_samp, _treat, prefix = mf_prefix):
                mageck_str = "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}{ctrlnm}-{sampnm}"
                # get the replicate strings
                ctrls = ','.join(sample_reps[_ctrl_samp])
                treats = ','.join(sample_reps[_treat])
                s = mageck_str.format(
                    counts=fn_counts,
                    treat=treats,
                    ctrl=ctrls,
                    outprefix=prefix,
                    ctrlnm=_ctrl_samp,
                    sampnm=_treat
                )
                pipeLOG.info(s)
                call(s.split())

            pipeLOG.info('Running MAGeCK version '+check_output(["mageck","-v"]).decode())

            for ctrl_samp, treat_samps in ctrlmap.items():
                for treat in treat_samps:
                    if treat == ctrl_samp:
                        continue
                    _run_mageck(ctrl_samp, treat)
                    mageck_pairs_done.append((ctrl_samp, treat))
                if not skip_extra_mageck:
                    # run all combinations of samples that share at least one control sample
                    for c, t in combinations(treat_samps, 2):
                        if (c,t) not in mageck_pairs_done and (t,c) not in mageck_pairs_done:
                            mageck_pairs_done.append((c,t))
                            # maybe don't tablulate these at the moment.
                            _run_mageck(c, t, prefix=mf_prefix.replace(ctrlgroup, 'EXTRA'))

        ############
        #* Get tables
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

def process_arguments(arguments):

    samples = arguments['sample_reps'].keys()
    controls = arguments['controls']
    for ctrl_grp, ctrlmap in controls.items():
        #todo: if ctrlmap is a list assume a single control group and don't crash
        #print('ctrlmap',  ctrlmap)
        for ctrl, samps in ctrlmap.items():
            if type(samps) == str:
                samps = [samps]
            if samps == ['ALL']:
                samps = copy(list(samples))
            elif type(samps) == dict:
                # maybe could have other keywords in the future
                samps = [s for s in samples if s not in samps['EXCEPT']]
            ctrlmap[ctrl] = samps
    arguments['controls'] = controls

    # in case strings are passed instead of lists
    for comp in arguments['comparisons']:
        print(arguments['comparisons'])
        for k, v in comp.items():
            if type(v) == str:
                comp[k] = [v]
    for k, v in arguments['sample_reps'].items():
        if type(v) == str:
            arguments['sample_reps'][k] = [v]
    return arguments

if __name__ == '__main__':

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
    parser.add_argument('--dont-log', action='store_true', dest='dont_log', default=None,
                        help="Don't write a log file.")

    # get arguments from the command line and the YAML
    clargs = parser.parse_args() # need to assign this before calling vars() for some reason
    cmd_args = vars(clargs)
    yml_args = yaml.safe_load(open(cmd_args['fn_yaml']))
    del cmd_args['fn_yaml']

    # over write yml_args with any specified in the command line
    for k, v in cmd_args.items():
        if v is not None:
            yml_args[k] = v

    args = process_arguments(yml_args)

    #print(args)

    run_analysis(**args)



    # run_analysis(**vars(clargs))
    # os.chdir('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/matylda_HS715-20')
    # run_analysis(**{'fn_counts': '/Users/johnc.thomas/Dropbox/crispr/counts_all/matylda_HS715-20.counts.tsv',
    #  'fn_repmap': 'matylda_HS715-20.repmap.3.xlsx', 'outdir': 'take5', 'file_prefix': 'HS715-20', 'labeldep': 30,
    #  'labelenr': 10, 'charts_only': False})