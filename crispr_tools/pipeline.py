import logging, os, pathlib
import gzip
import argparse
from subprocess import call
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from jacks.jacks_io import runJACKS
from jacks.infer import LOG as jacksLOG
jacksLOG.setLevel(logging.WARNING)

logging.getLogger('matplotlib').setLevel(logging.WARNING)

from crispr_tools.count_reads import count_reads, count_batch, map_counts
from crispr_tools.tools import plot_read_violins, plot_ROC, plot_volcano, tabulate_mageck, plot_volcano_from_mageck
from crispr_tools.jacks_tools import tablulate_score, scores_scatterplot, mahal_nocov

__version__ = '1.7.5b1'


#todo include log file and identification?
#todo add sample labels to repmaps for nice chart labels
#todo toy data for test

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
    sheeti = 0

    xl = pd.ExcelFile(fn_repmap)
    for ctrlgroup in xl.sheet_names:
        repmap = xl.parse(ctrlgroup)
        repmap.to_csv('tmp.repmap.tsv', '\t')
        runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                 norm_type=norm, outprefix=outprefix+'.'+ctrlgroup+'.')

def run_analysis(fn_counts, fn_repmap, outdir, file_prefix, labeldep = 20, labelenr = 20,
                 charts_only = False, skip_mageck = False, skip_jacks = False):

    # comps end up on

    # if fileprefix[-1] != '.':
    #     fileprefix += '.'

    # prep
    # if os.path.isdir(outdir):
    #     input('Output directory '+outdir+' already exists.\nPress enter to overwrite any existing results...')


    call(['mkdir', outdir])

    for analysis in ('jacks_mode', 'jacks_median', 'mageck'):
        call(['mkdir', str(Path(outdir, analysis))])
        call(['mkdir', str(Path(outdir, analysis, 'volcano'))])
        call(['mkdir', str(Path(outdir, analysis, 'scatter'))])
        call(['mkdir', str(Path(outdir, analysis, 'tables'))])
        call(['mkdir', str(Path(outdir, analysis, 'files'))])
        #call(['mkdir', str(Path(outdir, 'scattercharts'))])

    xl = pd.ExcelFile(fn_repmap)
    if 'repmaps' in xl.sheet_names:
        ctrls = xl.parse('repmaps', header=None).iloc[:, 0].values
    else:
        ctrls = [s for s in xl.sheet_names if any([x in s for x in('D2', 'pre', 'NT', 'DMSO')])]
    print('control groups:' , ctrls)
    for ctrlgroup in ctrls:
        fnprefix = file_prefix + '.' + ctrlgroup + '.'
        jmode_prefix = str(Path(outdir, 'jacks_mode', 'files', fnprefix))
        jmed_prefix = str(Path(outdir, 'jacks_median', 'files', fnprefix))

        mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))

        repmap = xl.parse(ctrlgroup, index_col=0)
        repmap.to_csv('tmp.repmap.tsv', '\t')

        if  not charts_only and not skip_jacks:
            runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                     norm_type='mode', outprefix=jmode_prefix)
            runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                     norm_type='median', outprefix=jmed_prefix)

        samp_map = repmap.groupby('samp').groups

        # comparison by sample names
        mag_pairs = list(set([(row.ctrl, row.samp) for _, row in repmap.iterrows() if row.ctrl != row.samp]))

        if not charts_only and not skip_mageck:
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
            scoresmode = tablulate_score(jmode_prefix)
            scoresmed = tablulate_score(jmed_prefix)
            analyses_used.extend([('jacks_mode', 'jacks_score',  scoresmode),
                                  ('jacks_median', 'jacks_score', scoresmed)])
            if not charts_only:
                scoresmode.to_excel(str(Path(outdir, 'jacks_mode', 'tables', fnprefix+'jacks_scores.xlsx')))
                scoresmed.to_excel(str(Path(outdir, 'jacks_median', 'tables', fnprefix + 'jacks_scores.xlsx')))

        if not skip_mageck:
            magtab = tabulate_mageck(mf_prefix)
            analyses_used.append(('mageck', 'lfc', magtab))
            if not charts_only:
                magtab.to_excel(str(Path(outdir, 'mageck', 'tables', fnprefix+'mageck_table.xlsx')))

        #volcano charts of all
        # use the universal one, have label options as part of args, maybe as a dict.
        for analysis, xkey, tab in analyses_used:
            #analysis = analysis
            for exp in tab.columns.levels[0]:
                if analysis.startswith('jacks'):
                    chart_title = '{} from {} ({} normalisation)'.format(exp, ctrlgroup, analysis)
                else:
                    # split by fnpref just in case it has a dash in it...
                    #print(exp)
                    ctrlnm, treatnm = exp.split('-')
                    chart_title = '{} from {} ({} analysis)'.format(treatnm, ctrlnm, analysis)
                plot_volcano(
                    xkey, 'fdr_log10', tab[exp], chart_title,
                    labeldep,labelenr,p_thresh=0.1,
                    outfn=str(Path(outdir, analysis, 'volcano', fnprefix+exp+'.'+analysis+'_volcano.png'))
                )
                plt.close()

            call(['tar', '-zcf', str(Path(outdir, analysis+'volcano_charts.tar.gz')) ,str(Path(outdir, analysis, 'volcano'))])

        #scattercharts using comparisons.csv

        #compstr = [s for s in xl.sheet_names if s.lower().startswith('comp')][0]
        comptab = xl.parse('sample_deets', index_col=0)
        if 'comps' not in comptab.columns:
            raise ValueError
        # except:
        #     print('** No comparisons')
        #     continue

        if not skip_jacks:
            for analysis, tab in ('jacks_mode', scoresmode), ('jacks_median', scoresmed):
                mahals = pd.DataFrame(index=tab.index)
                for samp, comps in comptab.iterrows():
                    if comps.comps is not np.NaN:
                        comps = comps.comps.replace(' ', '').split(',')
                        # comps generally controls
                        for comp in comps:
                            if comp not in tab.columns or samp not in tab.columns:
                                continue
                            # plot
                            scores_scatterplot(comp, samp, tab, label_pos=labelenr, label_neg=labeldep)
                            plt.savefig(str(
                                Path(outdir, analysis, 'scatter', fnprefix+"{}_vs_{}.scatter.png".format(comp, samp))
                            ), dpi=150)
                            plt.close()

                            # pop mahal table
                            A, B = tab[comp], tab[samp]
                            _, _, mahal = mahal_nocov(A['jacks_score'], A['stdev'], B['jacks_score'], B['stdev'])
                            mahals.loc[:, f"{comp} vs {samp}"] = mahal
                mahals.to_csv(Path(outdir, analysis, 'tables', fnprefix+f"mahalanobis_distances.csv"))
                # wipes over the one produced above if comptab exists.
                call(['tar', '-zcf', str(Path(outdir, analysis + 'scatter_charts.tar.gz')),
                      str(Path(outdir, analysis, 'scatter'))])
                call(['tar', '-zcf', str(Path(outdir, analysis + 'tables.tar.gz')),
                  str(Path(outdir, analysis, 'tables'))])

def iter_comps(comp_series: pd.Series, tab: pd.DataFrame):
    for samp, comps in comp_series.iteritems():
        if comps is not np.NaN:
            comps = comps.replace(' ', '').split(',')
            # comps generally controls
            for comp in comps:
                if comp not in tab.columns or samp not in tab.columns:
                    continue
                yield samp, comp

def test_pipeline(**kwargs):
    p = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    testp = str(p.parent/'tests')
    run_analysis(testp+'/test_counts.tsv', testp+'/test_pipeline.xlsx', testp+'/ran_test', 'itsatest',
             labeldep=10, labelenr=10, **kwargs)

if __name__ == '__main__':



    print(__version__)
    parser = argparse.ArgumentParser(description='Run mageck and jacks analyses using a repmap.xlsx and'
                                     " counts.tsv")
    parser.add_argument('fn_counts', metavar='COUNTS', help='Path to counts file')
    parser.add_argument('fn_repmap', metavar='REPMAP', help='Path to repmap that excel file')
    parser.add_argument('outdir', metavar='OUTDIR', help='Path to where results files will be stored, a '
                        "directory structure will be created.")
    parser.add_argument('file_prefix', metavar='PREFIX', help="String to form identifying prefix for all files generated.")
    parser.add_argument('-d', '--labeldep', metavar='N', type=int, default=10, help="Number of depleted genes to label"
                        "in charts.")
    parser.add_argument('-e', '--labelenr', metavar='N', type=int, default=10, help="Number of enriched genes to label"
                        "in charts.")
    parser.add_argument('--charts-only', action='store_true', default=False,
                        help="Don't run MAGeCK or JACKS, just produce plots with existing files.")
    parser.add_argument('--skip-jacks', action='store_true', dest='skip_jacks',
                        help="don't run JACKS analysis or try to plot from JACKS analyses")
    parser.add_argument('--skip-mageck', action='store_true', dest='skip_mageck',
                        help="don't run MAGeCK analysis or try to plot from MAGeCK analyses")
    clargs = parser.parse_args()
    print(vars(clargs))
    run_analysis(**vars(clargs))
    # os.chdir('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/matylda_HS715-20')
    # run_analysis(**{'fn_counts': '/Users/johnc.thomas/Dropbox/crispr/counts_all/matylda_HS715-20.counts.tsv',
    #  'fn_repmap': 'matylda_HS715-20.repmap.3.xlsx', 'outdir': 'take5', 'file_prefix': 'HS715-20', 'labeldep': 30,
    #  'labelenr': 10, 'charts_only': False})