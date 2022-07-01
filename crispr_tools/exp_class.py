import pandas as pd

np = pd.np
from typing import Union, List, Dict
import yaml
from crispr_tools import *
from attrdict import AttrDict
import os
from crispr_tools.crispr_pipeline import process_control_map
from crispr_tools.tools import clonal_lfcs
from copy import copy
import pathlib
class CrisprExperiment:
    def __init__(self,

                 count: Union[str, pd.DataFrame],
                 expd: Union[str, dict],
                 lib: Union[str, pd.DataFrame],
                 rep_deets: pd.DataFrame = None,
                 samp_deets: pd.DataFrame = None,
                 root_dir: os.PathLike = None):

        isStr = lambda x: type(x) is str

        if root_dir is None:
            root_dir = os.getcwd()
        self.root_dir = root_dir

        if isStr(expd):
            expd = yaml.load(open(expd))

        if isStr(count):
            count = pd.read_csv(count, '\t', index_col=0)

        self.count = count
        self.countng = count.drop('gene', 1)

        self.sample_replicates = expd['sample_reps']
        self.samples = self.sample_replicates.keys()

        expd['controls'] = process_control_map(expd['controls'], self.samples)
        self.controls = expd['controls']

        if 'labels' in expd:
            self.labels = expd['labels']
        else:
            self.labels = None

        self.expd = AttrDict(expd)

        # Might be good to have getter/setters associated with self.count that updates
        #   this and clonal counts whenever self.count changes.
        self.lncount = size_factor_normalise(count)

        self.clonal_lfc = AttrDict()
        for grp in self.controls:
            self.clonal_lfc[grp] = clonal_lfcs(
                self.lncount,
                self.controls[grp],
                self.sample_replicates
            )

        if isStr(lib):
            lib = pd.read_csv(lib)
            lib.set_index('guide', drop=False, inplace=True)
        self.lib = lib

        self.results = AttrDict()

        self.rep_deets = rep_deets
        self.samp_deets = samp_deets

    def add_mageck(self, prefix):
        res = AttrDict()
        for grp in self.controls.keys():
            res[grp] = tabulate_mageck(f"{prefix}.{grp}.")
        self.results['mageck'] = res

    def add_jacks(self, prefix):
        res = AttrDict()
        for grp in self.controls.keys():
            res[grp] = tabulate_score(f"{prefix}.{grp}.")
        self.results['jacks'] = res

    def add_results(self, expd=None):
        """Add analysis results from mageck or jacks, using info specified in
        expd."""
        if expd is None:
            expd = self.expd
        # todo will need fixing when job spec in expd is improved
        p = os.path.join(self.root_dir, expd['exp_name'], expd['analysis_name'])

        if 'skip_mageck' not in expd or not expd['skip_mageck']:
            prefix = os.path.join(p, 'mageck', 'files', f"{expd['file_prefix']}")
            self.add_mageck(prefix)

        if 'skip_jacks' not in expd or not expd['skip_jacks']:
            prefix = os.path.join(p, 'jacks_median', 'files', f"{expd['file_prefix']}")
            self.add_jacks(prefix)

# def load_exp_from_expd(expd: Union[str, bytes, os.PathLike, Dict, AttrDict], root_dir='./'):
#     if type(expd) is (str, bytes, os.PathLike):
#         with open(expd) as f:
#             expd = yaml.safe_load(f)
#         exp = CrisprExperiment(root_dir,
#                                expd['counts_fn'], )

#     # since clonal lfc is fairly expensive to calculate, lets leave it till we need it
#     @property
#     def clonal_lfc(self):

#         try:
#             return self._clonal_lfc
#         except AttributeError:
#             self._clonal_lfc = AttrDict()
#             for grp in self.controls:
#                 crispr_tools.qc.get_clonal_lfcs(
#                     self.lncount,
#                     self.controls[grp],
#                     self.sample_replicates
#                 )

#             return self._clonal_lfc




