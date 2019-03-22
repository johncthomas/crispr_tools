from .count_reads import count_reads, count_batch, map_counts
from .tools import size_factor_normalise, plot_read_violins, plot_ROC, plot_volcano, revcomp, plot_volcano_from_mageck, \
    tabulate_mageck
from .jacks_tools import tabulate_score, tabulate_jacks, plot_volcano_from_scoretable, scores_scatterplot
import pathlib
with open(pathlib.Path(__file__).parent/'version.txt') as f:
    __version__ = f.readline().replace('\n', '')