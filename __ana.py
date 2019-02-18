
from jacks_tools import *
from crispr_tools import *
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])
import logging

pltlogger = logging.getLogger('matplotlib')
pltlogger.setLevel(logging.WARNING)

from pathlib import Path, PosixPath, WindowsPath

from typing import List

os.chdir('/Users/johnc.thomas/Dropbox/adrst/')

a_d2, a_nt = [tablulate_score('ana/take3/jacks_median/files/HS714.'+s) for s in ('D2.', 'NT.')]
scores_scatterplot(a_d2.NT.jacks_score, a_d2['6TG'].jacks_score, min_mahal=0.5)