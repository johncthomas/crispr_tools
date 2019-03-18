#!python3
import logging, os, pathlib, datetime, inspect, argparse

from subprocess import call
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from jacks.jacks_io import runJACKS

with open(pathlib.Path(__file__).parent/'version.txt') as f:
    __version__ = f.readline().replace('\n', '')

pipeLOG = logging.getLogger('pipeline')
pipeLOG.setLevel(logging.INFO)

from jacks.infer import LOG as jacksLOG
jacksLOG.setLevel(logging.WARNING)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

from crispr_tools.count_reads import count_reads, count_batch, map_counts
from crispr_tools.tools import plot_read_violins, plot_ROC, plot_volcano, tabulate_mageck, plot_volcano_from_mageck
from crispr_tools.jacks_tools import tabulate_score, scores_scatterplot, mahal_nocov

#todo use YAML instead of Excel books (excel reading isn't default part of pandas...)
#todo don't ouput xlsx for comps
#todo add sample labels to sample_details for nice chart labels


"""Go from FastQ files to completed JACKS/MAGeCK analysis. 
fq->mapped_counts&violin plots are one command (count_reads.count_batch() ) 
counts->charts&tables another."""

mageck_str = "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}{ctrlnm}-{sampnm}"

def call_magecks(fn_counts, fn_repmap, outprefix):
    mageck_str = "mageck test -k {counts} -t {treat} -c {ctrl} -n {outprefix}{ctrlnm}-{sampnm}"
    xl = pd.ExcelFile(fn_repmap)
    for ctrlgroup in xl.sheet_names:
        repmap = xl.parse(ctrlgroup)

        samp_map = {}
        for samp, groupi in repmap.groupby('samp').groups.items():
            samp_map[samp] = [repmap.iloc[i, 0] for i in groupi]

        # comparison by sample names
        comps = list(set([(row.ctrl, row.samp) for _, row in repmap.iterrows()]))

        for ctrl, treat in comps:
            ctrls = ','.join(samp_map[ctrl])
            treats = ','.join(samp_map[treat])
            s = mageck_str.format(
                counts=fn_counts,
                treat=treats,
                ctrl=ctrls,
                outprefix=outprefix+'.'+ctrlgroup+'.',
                ctrlnm=ctrl,
                sampnm=treat
            )
            call(s, shell=True, )


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


def run_analysis(fn_counts, fn_design_sheet, outdir, file_prefix, volc_labels=15, scatter_labels=10,
                 charts_only = False, jacks_eff_fn=None,
                 skip_mageck = False, skip_jacks = False, dont_log=False):

    pipeLOG.info("Version = "+__version__)

    call(['mkdir', outdir])
    if not dont_log:
        t = str(datetime.datetime.now()).replace(' ', '_').replace('/', '|')
        set_logger(str(Path(outdir, file_prefix + f'log_{t}.txt')))

    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    pipeLOG.info('\n\t'.join([f"{x} = {values[x]}" for x in args]))
    pipeLOG.info(f"Full outdir = {os.path.realpath(outdir)}")

    for analysis_str in ('jacks_mode', 'jacks_median', 'mageck'):
        call(['mkdir', str(Path(outdir, analysis_str))])
        call(['mkdir', str(Path(outdir, analysis_str, 'volcano'))])
        call(['mkdir', str(Path(outdir, analysis_str, 'scatter'))])
        call(['mkdir', str(Path(outdir, analysis_str, 'tables'))])
        call(['mkdir', str(Path(outdir, analysis_str, 'files'))])
        #call(['mkdir', str(Path(outdir, 'scattercharts'))])

    xl = pd.ExcelFile(fn_design_sheet)
    if 'repmaps' in xl.sheet_names:
        ctrl_groups = xl.parse('repmaps', header=None).iloc[:, 0].values
    else:
        ctrl_groups = [s for s in xl.sheet_names if any([x in s for x in('D2', 'pre', 'NT', 'DMSO')])]

    pipeLOG.info('Control groups being used: '+', '.join(ctrl_groups))

    comptab = xl.parse('sample_details', index_col=0)
    if 'comps' not in comptab.columns:
        raise ValueError


    for ctrlgroup in ctrl_groups:
        fnprefix = file_prefix + '.' + ctrlgroup + '.'
        pipeLOG.info(f'\nRunning {ctrlgroup}, writing to {fnprefix}')
        jmode_prefix = str(Path(outdir, 'jacks_mode', 'files', fnprefix))
        jmed_prefix = str(Path(outdir, 'jacks_median', 'files', fnprefix))

        mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))

        repmap = xl.parse(ctrlgroup, index_col=0)
        repmap.to_csv('tmp.repmap.tsv', '\t')

        if not charts_only and not skip_jacks:
            pipeLOG.info("Running JACKS, mode and median normalised")
            if jacks_eff_fn:
                jacks_eff = jacks_eff_fn.format(ctrlgroup)
                if os.path.isfile(jacks_eff):
                    pipeLOG.info(f'\tUsing {jacks_eff} for efficacy')
                else:
                    pipeLOG.info(f"\tcould not find {jacks_eff}, calculating efficacy")
                    jacks_eff = None
            else:
                jacks_eff = None

            runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                     norm_type='mode', outprefix=jmode_prefix, reffile=jacks_eff)
            runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                     norm_type='median', outprefix=jmed_prefix, reffile=jacks_eff)

        # comparison by sample names
        mag_pairs = list(set([(row.ctrl, row.samp) for _, row in repmap.iterrows() if row.ctrl != row.samp]))
        # group samples from the repmap, call jacks with groups of samples
        samp_map = repmap.groupby('samp').groups
        if not charts_only and not skip_mageck:
            pipeLOG.info('Running MAGeCK')
            for ctrl, treat in mag_pairs:
                ctrls = ','.join(samp_map[ctrl])
                treats = ','.join(samp_map[treat])
                s = mageck_str.format(
                    counts=fn_counts,
                    treat=treats,
                    ctrl=ctrls,
                    outprefix=mf_prefix,
                    ctrlnm=ctrl,
                    sampnm=treat
                )
                call(s, shell=True, )

        analyses_used = []
        # get results tables of all
        if not skip_jacks:
            scoresmode = tabulate_score(jmode_prefix)
            scoresmed = tabulate_score(jmed_prefix)
            analyses_used.extend([('jacks_mode', 'jacks_score',  scoresmode),
                                  ('jacks_median', 'jacks_score', scoresmed)])
            if not charts_only:
                pipeLOG.info('Writing JACKS tables')
                scoresmode.to_csv(str(Path(outdir, 'jacks_mode', 'tables', fnprefix+'jacks_scores.csv')))
                scoresmed.to_csv(str(Path(outdir, 'jacks_median', 'tables', fnprefix + 'jacks_scores.csv')))

        if not skip_mageck:
            scores_mageck = tabulate_mageck(mf_prefix)
            analyses_used.append(('mageck', 'lfc', scores_mageck))
            if not charts_only:
                pipeLOG.info('Writing MAGeCK tables')
                scores_mageck.to_csv(str(Path(outdir, 'mageck', 'tables', fnprefix+'mageck_table.csv')))

        #volcano charts of all
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

            call(['tar', '-zcf',
                  str(Path(outdir, analysis_str+'volcano_charts.tar.gz')), #file to be created
                  # change to the target dir so it doesn't include the path in the tar
                  '-C', str(Path(outdir, analysis_str, 'volcano')), '.'])


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
            for comp_to, comp_from in iter_comps(comptab.comps, analysis_tab):
                comparisons_done = True
                pipeLOG.info(f"\t{comp_from} vs {comp_to}")
                #scatter


                scores_scatterplot(comp_from, comp_to, analysis_tab,
                                   label_pos=scatter_labels, label_neg=scatter_labels,
                                   dist_name='Mahalanobis distance')
                plt.title(f"{comp_from} vs {comp_to} ({ctrlgroup}, JACKS)")
                plt.savefig(str(
                    Path(outdir, analysis_str, 'scatter', fnprefix+"{}_vs_{}.scatter.png".format(comp_from, comp_to))
                ), dpi=150)
                plt.close()

                # pop mahal table
                A, B = analysis_tab[comp_from], analysis_tab[comp_to]
                _, _, mahal = mahal_nocov(A['jacks_score'], A['stdev'], B['jacks_score'], B['stdev'])
                mahals.loc[:, f"{comp_from} vs {comp_to}"] = mahal

            if comparisons_done:
                mahals.to_csv(Path(outdir, analysis_str, 'tables', fnprefix+f"mahalanobis_distances.csv"))

            # wipes over the one produced above if comptab exists.
            call(['tar', '-zcf', str(Path(outdir, analysis_str + 'scatter_charts.tar.gz')),
                  '-C', str(Path(outdir, analysis_str, 'scatter')), '.'])
                    #str(Path(outdir, analysis_str, 'scatter'))])
            call(['tar', '-zcf', str(Path(outdir, analysis_str + 'tables.tar.gz')),
                  '-C', str(Path(outdir, analysis_str, 'tables')), '.'])
            #str(Path(outdir, analysis_str, 'tables'))])

    # Do all the mageck comparisons at the end as we need the NT->treat for significance
    if not skip_mageck:
        pipeLOG.info('Doing comparisons of mageck')
        mag_tables = {}
        # Get the tables
        for ctrlgroup in ctrl_groups:
            fnprefix = file_prefix + '.' + ctrlgroup + '.'
            mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))
            scores_mageck = tabulate_mageck(mf_prefix)
            mag_tables[ctrlgroup] = scores_mageck

        # go through each comparison pair and check for the existence of the right measures.
        for comp_to, comp_from in iter_comps(comptab.comps):
            lfc_tabs = []
            fdr_tabs = []
            for ctrlgroup, tab in mag_tables.items():
                # second of column headers is the sample in mageck
                exp = pd.Series(tab.columns.levels[0], index=tab.columns.levels[0])
                samps = set([c.split('-')[1] for c in exp])
                # so, if a comp is specified with NT->T then we want to use NT->T for fdr and
                # ctrl->[NT/T] LFCs for the scatter, then name by ctrl.NT_T
                if comp_from in samps and comp_to in samps:
                    samphead = exp[exp.str.contains(comp_to)][0]
                    comphead = exp[exp.str.contains(comp_from)][0]
                    # some weird groups will exist as "ctrlgroup", so check the ctrls actually match
                    if comphead.split('-')[0] == samphead.split('-')[0]:
                        lfc_tabs.append((ctrlgroup, samphead, comphead))
                if f'{comp_to}-{comp_from}' in exp:
                    fdr_tabs.append((ctrlgroup, f'{comp_to}-{comp_from}'))
                elif f'{comp_from}-{comp_to}' in exp:
                    fdr_tabs.append((ctrlgroup, f'{comp_from}-{comp_to}'))
                else: # no source of sig
                    pipeLOG.info(f'{comp_from}->{comp_to} not possible, missing comparison.')
                    continue
            # print(lfc_tabs)
            # print('**', fdr_tabs)
            for fdrtab, fdrexp in fdr_tabs:
                for ctrlgroup, samphead, comphead in lfc_tabs:
                    pipeLOG.info(f'MAgeck results comp using {samphead} {comphead}, {fdrexp}')

                    scores_scatterplot(comphead, samphead, mag_tables[ctrlgroup], True, scatter_labels, scatter_labels,
                                       distance=mag_tables[fdrtab].loc[:, (fdrexp, 'fdr_log10')],
                                       min_label_dist=0.3,
                                       dist_name='log10(FDR)')
                    plt.title(f"{comphead} vs {samphead} (MAGeCK)")
                    plt.savefig(str(
                        Path(outdir, 'mageck', 'scatter', file_prefix+ '.' + ctrlgroup + ".{}_vs_{}.scatter.png".format(comp_from, comp_to))
                    ), dpi=150)
                    plt.close()

    pipeLOG.info("Done. "+str(datetime.datetime.now()))

def iter_comps(comp_series: pd.Series, tab: pd.DataFrame=None):
    for samp, comps in comp_series.iteritems():
        if comps is not np.NaN:
            comps = comps.replace(' ', '').split(',')
            # comps generally controls
            for comp in comps:
                if tab is not None:
                    if comp not in tab.columns or samp not in tab.columns:
                        continue
                yield samp, comp


if __name__ == '__main__':
    print(__version__)
    parser = argparse.ArgumentParser(
        description="Run mageck and jacks analyses using a design_sheet.xlsx and counts.tsv"
    )
    parser.add_argument('fn_counts', metavar='COUNTS', help='Path to counts file')
    parser.add_argument('fn_design_sheet', metavar='REPMAP',
                        help='Path to excel file containing repmaps etc.')
    parser.add_argument('outdir', metavar='OUTDIR', help='Path to where results files will be stored, a '
                        "directory structure will be created.")
    parser.add_argument('file_prefix', metavar='PREFIX', help="String to form identifying prefix for all files generated.")
    parser.add_argument('-v', '--volc-labels',dest='volc_labels', metavar='N', type=int, default=30,
                        help="Number of depleted/enriched genes to label in volcano charts.")
    parser.add_argument('-s', '--scatter-labels', metavar='N', type=int, default=10, dest='scatter_labels',
                        help="Number of genes to label on either side of The Line in scatter charts.")
    parser.add_argument('--charts-only', action='store_true', default=False,
                        help="Don't run MAGeCK or JACKS, just produce plots with existing files.")
    parser.add_argument('-j', '--jacks-efficacy', metavar='EFF_FILE', dest='jacks_eff_fn',
                        help='Path to efficacy file to be used in JACKS analysis. '
                             'If the provided string contains {}, the repmap ctrl group will be substituted'
                             'into the filename. Optional.')
    parser.add_argument('--skip-jacks', action='store_true', dest='skip_jacks',
                        help="don't run JACKS analysis or try to plot from JACKS analyses")
    parser.add_argument('--skip-mageck', action='store_true', dest='skip_mageck',
                        help="don't run MAGeCK analysis or try to plot from MAGeCK analyses")
    parser.add_argument('--dont-log', action='store_true', dest='dont_log',
                        help="Don't write a log file.")
    clargs = parser.parse_args()
    run_analysis(**vars(clargs))
    # os.chdir('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/matylda_HS715-20')
    # run_analysis(**{'fn_counts': '/Users/johnc.thomas/Dropbox/crispr/counts_all/matylda_HS715-20.counts.tsv',
    #  'fn_repmap': 'matylda_HS715-20.repmap.3.xlsx', 'outdir': 'take5', 'file_prefix': 'HS715-20', 'labeldep': 30,
    #  'labelenr': 10, 'charts_only': False})