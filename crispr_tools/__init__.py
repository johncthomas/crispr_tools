from .count_reads import get_count_table_from_file_list, count_batch, map_counts
from .tools import (size_factor_normalise,  revcomp,
    tabulate_mageck, tabulate_drugz,  hart_list )
from .plotting import plot_read_violins, plot_ROC, plot_volcano,plot_volcano_from_mageck, pca_grid, scores_scatter_plotly
from .exp_class import CrisprExperiment
from .version import __version__


# import pathlib
# with open(pathlib.Path(__file__).parent/'version.txt') as f:
#     __version__ = f.readline().replace('\n', '')