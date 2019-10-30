print(__file__)

# import logging
# import sys
# import seaborn as sns
# from jacks_tools import tabulate_score
# prefix = '~/Dropbox/crispr/screens_analysis/matylda_742-3^miniscreen-dubv2/take1/jacks_median/files/mat_742-3.D2.'
# from jacks_tools import *
#
# genes = pd.read_table(prefix + '_gene_JACKS_results.txt', sep='\t', index_col=0)
# genes_index = sorted(genes.index)
# genes = genes.reindex(genes_index)
# genesstd = pd.read_table(prefix + '_gene_std_JACKS_results.txt', sep='\t', index_col=0)
# genesstd = genesstd.reindex(genes_index)
# ps = genes / genesstd
# ps = ps.apply(norm.cdf)
#
# # multiindex DF for each experiment giving results
# sig_cols = pd.MultiIndex.from_product((ps.columns, ['jacks_score', 'fdr_pos', 'fdr_neg', 'fdr_log10', 'stdev']),
#                                       names=('exp', 'stat'))
# sig_df = pd.DataFrame(index=genes_index, columns=sig_cols)
#
# for exp in ps.columns:
#    sig_df.loc[:, (exp, 'fdr_neg')] = multipletests(ps[exp], method='fdr_bh')[1]
#    sig_df.loc[:, (exp, 'fdr_pos')] = multipletests(1 - ps[exp], method='fdr_bh')[1]
#    sig_df.loc[:, (exp, 'jacks_score')] = genes[exp]
#    sig_df.loc[:, (exp, 'stdev')] = genesstd[exp]
#    score_table = sig_df[exp].copy()
#    # get one FDR
#    pos = score_table['jacks_score'] > 0
#    sns.distplot(score_table['jacks_score'], bins=50)
#    plt.show()
#    print(sum(score_table['jacks_score'] > 0))
#    break
#
#
#
# # class StreamToLogger(object):
# #    """
# #    Fake file-like stream object that redirects writes to a logger instance.
# #    """
# #    def __init__(self, logger, log_level=logging.INFO):
# #       self.logger = logger
# #       self.log_level = log_level
# #       self.linebuf = ''
# #
# #    def write(self, buf):
# #       for line in buf.rstrip().splitlines():
# #          self.logger.log(self.log_level, line.rstrip())
# #
# # # logging.basicConfig(
# # #    level=logging.DEBUG,
# # #    format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
# # #    filename="out.log",
# # #    filemode='a'
# # # )
# #
# # import os
# # print(os.getcwd())
# #
# # pipeLOG = logging.getLogger('pipeline')
# #
# # pipeLOG.setLevel(logging.INFO)
# # log_fn = 'test_logsterr.txt'
# # hndlr = logging.FileHandler(log_fn, 'w')
# # pipeLOG.setLevel(logging.INFO)
# # pipeLOG.addHandler(hndlr)
# #
# # sys.stderr = StreamToLogger(pipeLOG, logging.ERROR)
# # pipeLOG.info('fffffaaa')
# # raise Exception('Test to standard error')
#
#
# # import yaml, copy
# # def process_yaml(fn):
# #     # arguments = yaml.safe_load(open('/Users/johnc.thomas/Dropbox/crispr/pycharm/jts_crispr_processing/crispr_tools_yaml/tests/test_pipeline.yaml'))
# #     arguments = yaml.safe_load(open(fn))
# #     # print(arguments)
# #     # repmaps = arguments['repmap']
# #     samples = arguments['sample_reps'].keys()
# #     controls = arguments['controls']
# #     for ctrl_grp, ctrlmap in controls.items():
# #
# #         for ctrl, samps in ctrlmap.items():
# #             if samps == 'ALL':
# #                 samps = copy(list(samples))
# #             elif type(samps) == dict:
# #                 # maybe could have other keywords in the future
# #                 samps = [s for s in samples if s not in samps['EXCEPT']]
# #             else:
# #                 if type(samps) == str:
# #                     samps = [samps]
# #             ctrlmap[ctrl] = samps
# #     arguments['controls'] = controls
# #
# #     # in case strings are passed instead of lists
# #     for comp in arguments['comparisons']:
# #
# #         for k, v in comp.items():
# #             if type(v) == str:
# #                 comp[k] = [v]
# #     for k, v in arguments['sample_reps'].items():
# #         if type(v) == str:
# #             arguments['sample_reps'][k] = [v]
# #     return arguments
# #
# # process_yaml('paco-maryam_628-9_separate.yaml')
#
# # import os, sys
# # sys.path.insert(0, '.')
# #
# # CURRENT_DIR = os.path.dirname(__file__)
# #
# # # We can not import `xgboost.libpath` in setup.py directly since xgboost/__init__.py
# # # import `xgboost.core` and finally will import `numpy` and `scipy` which are setup
# # # `install_requires`. That's why we're using `exec` here.
# # libpath_py = os.path.join(CURRENT_DIR, 'crispr_tools/tests/test_counts.tsv')
# # libpath = {'__file__': libpath_py}
# # print(libpath)
#
# #exec(compile(open(libpath_py, "rb").read(), libpath_py, 'exec'), libpath, libpath)
#
# #LIB_PATH = libpath['find_lib_path']()
# #print("Install libxgboost from: %s" % LIB_PATH)
# #print([os.path.relpath(p) for p in LIB_PATH])
#
#
# # import numpy as np
# # #from crispr_tools.jacks_tools import scores_scatterplot
# # from crispr_tools import *
# # import matplotlib.pyplot as plt
# #
# # import logging, os, pathlib, datetime, inspect, argparse
# #
# # from subprocess import call
# # from pathlib import Path
# #
# # import numpy as np
# # import pandas as pd
# # import matplotlib.pyplot as plt
# #
# # os.chdir('/Users/johnc.thomas/Dropbox/crispr/screens_analysis/')
# #
# # #fn_design_sheet = "rebecca.601-2_608-9/rebecca.601-2_608-9.experiment_design.xlsx"
# # fn_design_sheet = 'paco-maryam.628-9/paco-maryam.628-9.experiment_design.xlsx'
# # outdir = "paco-maryam.628-9/take1"
# # file_prefix = "paco-maryam.628-9"
# # volc_labels = 30
# # scatter_labels = 10
# # charts_only = True
# # jacks_eff_fn = None
# # skip_mageck = False
# # skip_jacks = True
# # dont_log = False
# #
# # xl = pd.ExcelFile(fn_design_sheet)
# # if 'repmaps' in xl.sheet_names:
# #     ctrl_groups = xl.parse('repmaps', header=None).iloc[:, 0].values
# # else:
# #     ctrl_groups = [s for s in xl.sheet_names if any([x in s for x in ('D2', 'pre', 'NT', 'DMSO')])]
# #
# # #pipeLOG.info('Control groups being used: ' + ', '.join(ctrl_groups))
# #
# # comptab = xl.parse('sample_details', index_col=0)
# # if 'comps' not in comptab.columns:
# #     raise ValueError
# #
# #
# # def iter_comps(comp_series: pd.Series, tab: pd.DataFrame=None):
# #     for samp, comps in comp_series.iteritems():
# #         if comps is not np.NaN:
# #             comps = comps.replace(' ', '').split(',')
# #             # comps generally controls
# #             for comp in comps:
# #                 if tab is not None:
# #                     if comp not in tab.columns or samp not in tab.columns:
# #                         continue
# #                 yield samp, comp
# #
# # if not skip_mageck:
# #     #pipeLOG.info('Doing comparisons of mageck')
# #     mag_tables = {}
# #     # Get the tables
# #     for ctrlgroup in ctrl_groups:
# #         fnprefix = file_prefix + '.' + ctrlgroup + '.'
# #         mf_prefix = str(Path(outdir, 'mageck', 'files', fnprefix))
# #         scores_mageck = tabulate_mageck(mf_prefix)
# #         mag_tables[ctrlgroup] = scores_mageck
# #
# #     # go through each comparison pair and check for the existence of the right measures.
# #     # The right measures are A->B for FDR, and some_ctrl->A & some_ctrl->B for LFC
# #     for comp_to, comp_from in iter_comps(comptab.comps):
# #         lfc_tabs = []
# #         fdr_tabs = []
# #         for ctrlgroup, tab in mag_tables.items():
# #             # second of column headers is the sample in mageck
# #             exp = pd.Series(tab.columns.levels[0], index=tab.columns.levels[0])
# #             samps = set([c.split('-')[1] for c in exp])
# #             # so, if a comp is specified with NT->T then we want to use NT->T for fdr and
# #             # ctrl->[NT/T] LFCs for the scatter, then name by ctrl.NT_T
# #             if comp_from in samps and comp_to in samps:
# #                 samphead = exp[exp.str.contains(comp_to)][0]
# #                 comphead = exp[exp.str.contains(comp_from)][0]
# #                 # some weird groups will exist as "ctrlgroup", so check the ctrls actually match
# #                 if comphead.split('-')[0] == samphead.split('-')[0]:
# #                     lfc_tabs.append((ctrlgroup, samphead, comphead))
# #             if f'{comp_to}-{comp_from}' in exp:
# #                 fdr_tabs.append((ctrlgroup, f'{comp_to}-{comp_from}'))
# #             elif f'{comp_from}-{comp_to}' in exp:
# #                 fdr_tabs.append((ctrlgroup, f'{comp_from}-{comp_to}'))
# #             else:  # no source of sig
# #                 #pipeLOG.info(f'{comp_from}->{comp_to} not possible, missing comparison.')
# #                 #print('noooop')
# #                 continue
# #         # print('comp', comp_from, comp_to)
# #         # print('lfc', lfc_tabs)
# #         # print('fdr', fdr_tabs)
# #         # print(lfc_tabs)
# #         # print('**', fdr_tabs)
# #         for fdrtab, fdrexp in fdr_tabs:
# #             for ctrlgroup, samphead, comphead in lfc_tabs:
# #                 fn = file_prefix + '.' + ctrlgroup + ".{}_vs_{}.scatter.png".format(comp_from, comp_to)
# #                 call(f"cp {outdir}/mageck/scatter_/{fn} {outdir}/mageck/scatter/{fn} ", shell=True)
# #                 #      shell=True)
# #                 pass
# #                 #pipeLOG.info(f'MAgeck results comp using {samphead} {comphead}, {fdrexp}')
# #                 #
# #                 # scores_scatterplot(comphead, samphead, mag_tables[ctrlgroup], True, 1,1,
# #                 #                    distance=mag_tables[fdrtab].loc[:, (fdrexp, 'fdr_log10')],
# #                 #                    min_label_dist=0.3)
# #                 # plt.title(f"{comphead} vs {samphead} (MAGeCK) jjjjjjjjj")
# #                 # plt.savefig(str(
# #                 #     Path(outdir, 'mageck', 'scatter',
# #                 #          file_prefix + '.' + ctrlgroup + ".{}_vs_{}.scatter.png".format(comp_from, comp_to))
# #                 # ), dpi=150)
# #                 # plt.close()
# #
# #
# #
# # # p = '/Users/johnc.thomas/Dropbox/crispr/screens_analysis/rebecca.601-2_608-9/take1/mageck/files/rebecca.601-9.{}.'
# # #
# # # mag_tables = {}
# # # for grp in 'D2', 'pretreat', 'NT':
# # #     mag_tables[grp] = tabulate_mageck(p.format(grp))
# # #
# # # print(mag_tables['D2'].columns.levels[0])
# # #
# # # ctrlgroup='D2'
# # # comphead = 'A_48h-A_DMSO'
# # # samphead = 'A_48h-A_IC10'
# # # scores_scatterplot(comphead, samphead, mag_tables[ctrlgroup], True, 10, 10,
# # #                    distance=mag_tables['NT'].loc[:, ('A_DMSO-A_IC10', 'fdr_log10')],
# # #                    min_label_dist=0.3)
# # # plt.show()