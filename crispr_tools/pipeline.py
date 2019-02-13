import logging, os
import gzip
import argparse
from subprocess import call
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from jacks.jacks_io import runJACKS
from jacks.infer import LOG as jacksLOG
jacksLOG.setLevel(logging.WARNING)

logging.getLogger('matplotlib').setLevel(logging.WARNING)

from crispr_tools.count_reads import count_reads, count_batch, map_counts
from crispr_tools.tools import plot_read_violins, plot_ROC, plot_volcano, tabulate_mageck, plot_volcano_from_mageck
from crispr_tools.jacks_tools import tablulate_score, tabulate_jacks, plot_volcano_from_scoretable, scores_scatterplot


#todo fix chart titles for mag volc, and x title for jacks volc

#todo put zips in the analysis root
#todo include log file and identification?
#todo add sample labels to repmaps for nice chart labels
#todo do jacks mode and median norm by default, or choose the appropriate one


"""Go from FastQ files to completed JACKS/MAGeCK analysis. 
fq->mapped_counts&violin plots are one command 
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

def run_analysis(fn_counts, fn_repmap, outdir, file_prefix, labeldep = 30, labelenr = 10,
                 charts_only = False):
    # if fileprefix[-1] != '.':
    #     fileprefix += '.'

    # prep
    # if os.path.isdir(outdir):
    #     input('Output directory '+outdir+' already exists.\nPress enter to overwrite any existing results...')
    call(['mkdir', outdir])

    for resd in ('jacks_mode', 'jacks_median', 'mageck'):
        call(['mkdir', str(Path(outdir, resd))])
        call(['mkdir', str(Path(outdir, resd, 'charts'))])
        call(['mkdir', str(Path(outdir, resd, 'tables'))])
        call(['mkdir', str(Path(outdir, resd, 'files'))])

    xl = pd.ExcelFile(fn_repmap)
    for ctrlgroup in xl.sheet_names:
        fnprefix = file_prefix + '.' + ctrlgroup + '.'
        jmode_prefix = str(Path(outdir, 'jacks_mode', 'files', fnprefix))
        jmed_prefix = str(Path(outdir, 'jacks_median', 'files', fnprefix))

        mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))

        repmap = xl.parse(ctrlgroup, index_col=0)
        repmap.to_csv('tmp.repmap.tsv', '\t')

        if not charts_only:
            runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                     norm_type='mode', outprefix=jmode_prefix)
            runJACKS(fn_counts, 'tmp.repmap.tsv', fn_counts, 'rep', 'samp', None, 'ctrl', 'guide', 'gene',
                     norm_type='median', outprefix=jmed_prefix)

        samp_map = repmap.groupby('samp').groups

        # comparison by sample names
        comps = list(set([(row.ctrl, row.samp) for _, row in repmap.iterrows() if row.ctrl != row.samp]))

        if not charts_only:
            for ctrl, treat in comps:
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

        # get results tables of all
        scoresmode = tablulate_score(jmode_prefix)
        scoresmed = tablulate_score(jmed_prefix)
        magtab = tabulate_mageck(mf_prefix)

        if not charts_only:
            # write tables
            scoresmode.to_excel(str(Path(outdir, 'jacks_mode', 'tables', fnprefix+'jacks_scores.xlsx')))
            scoresmed.to_excel(str(Path(outdir, 'jacks_median', 'tables', fnprefix + 'jacks_scores.xlsx')))
            magtab.to_excel(str(Path(outdir, 'mageck', 'tables', fnprefix+'mageck_table.xlsx')))

        #volcano charts of all
        # use the universal one, have label options as part of args, maybe as a dict.
        for annm, xkey, tab in ('jacks_mode', 'jacks_score',  scoresmode), \
                               ('jacks_median', 'jacks_score', scoresmed), \
                               ('mageck', 'lfc', magtab):
            resd = annm
            for exp in tab.columns.levels[0]:
                if annm.startswith('jacks'):
                    chart_title = '{} vs {} ({} normalisation)'.format(exp, ctrlgroup, annm)
                else:
                    ctrlnm, treatnm = exp.split('-')
                    chart_title = '{} vs {} ({} analysis)'.format(treatnm, ctrlnm, annm)

                plot_volcano(
                    xkey, 'fdr_log10', tab[exp], chart_title,
                    labeldep,labelenr,p_thresh=0.1,
                    outfn=str(Path(outdir, resd, 'charts', fnprefix+exp+'.'+annm+'_volcano.png'))
                )
                plt.close()

            call(['tar', '-zcf', str(Path(outdir, resd+'_charts.tar.gz')) ,str(Path(outdir, resd, 'charts'))])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run mageck and jacks analyses using a repmap.xlsx and'
                                     " counts.tsv")
    parser.add_argument('fn_counts', metavar='COUNTS', help='Path to counts file')
    parser.add_argument('fn_repmap', metavar='REPMAP', help='Path to repmap that excel file')
    parser.add_argument('outdir', metavar='OUTDIR', help='Path to where results files will be stored, a '
                        "directory structure will be created.")
    parser.add_argument('file_prefix', metavar='PREFIX', help="String to form identifying prefix for all files generated.")
    parser.add_argument('-d', '--labeldep', metavar='N', type=int, default=30, help="Number of depleted genes to label"
                        "in charts.")
    parser.add_argument('-e', '--labelenr', metavar='N', type=int, default=10, help="Number of enriched genes to label"
                        "in charts.")
    parser.add_argument('--charts-only', action='store_true', default=False,
                        help="Don't run MAGeCK or JACKS, just produce plots with existing files.")
    clargs = parser.parse_args()
    #print(vars(clargs))
    run_analysis(**vars(clargs))