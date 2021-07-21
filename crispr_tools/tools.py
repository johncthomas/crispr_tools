import os, itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import yaml

try:
    from adjustText import adjust_text
except ModuleNotFoundError:
    def adjust_text(*args, **kwargs):
        pass
from typing import List, Union, Dict, Collection, Iterable, Tuple
PathType = Union[str, bytes, os.PathLike]
from pathlib import Path, PosixPath, WindowsPath
import pandas as pd
from itertools import zip_longest
from attrdict import AttrDict
import statsmodels.api as sm
OLS = sm.regression.linear_model.OLS
import xlsxwriter

#import multiprocessing as mp

hart_list = ['AARS1', 'ABCE1', 'ABCF1', 'ACTB', 'ACTL6A', 'ACTR10', 'ACTR2',
       'ADSL', 'ADSS2', 'AHCY', 'ALG1', 'ALG14', 'ALG2', 'ANAPC2',
       'ANAPC4', 'ANAPC5', 'AQR', 'ARCN1', 'ARIH1', 'ARL2', 'ATP2A2',
       'ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5MF-PTCD1',
       'ATP5MG', 'ATP5PO', 'ATP6V0B', 'ATP6V0C', 'ATP6V1A', 'ATP6V1D',
       'ATP6V1E1', 'ATR', 'AURKB', 'BANF1', 'BIRC5', 'BUB1B', 'BUB3',
       'BUD31', 'BYSL', 'TWNK', 'C1orf109', 'CFAP298', 'NEPRO', 'SPOUT1',
       'CCDC84', 'YJU2', 'CCNA2', 'CCNH', 'CCNK', 'CCT2', 'CCT3', 'CCT4',
       'CCT5', 'CCT6A', 'CCT7', 'CCT8', 'CDC123', 'CDC16', 'CDC20',
       'CDC27', 'CDC37', 'CDC5L', 'CDC73', 'CDK1', 'CDK7', 'CDK9', 'CDT1',
       'CEBPZ', 'CENPA', 'CENPC', 'CFL1', 'CHAF1A', 'CHAF1B', 'CHEK1',
       'CHERP', 'CHMP2A', 'CHMP6', 'CIAO1', 'CINP', 'UTP4', 'CKAP5',
       'CLNS1A', 'CLP1', 'CLTC', 'CMPK1', 'CMTR1', 'CNOT3', 'COA5',
       'COPA', 'COPB1', 'COPB2', 'COPS3', 'COPS6', 'COPZ1', 'COQ4',
       'COX10', 'COX11', 'COX15', 'COX4I1', 'COX5B', 'COX6B1', 'CPSF1',
       'CPSF2', 'CPSF3', 'CPSF4', 'CRNKL1', 'CSE1L', 'CTDP1', 'CTPS1',
       'CTR9', 'CYCS', 'DAD1', 'DBR1', 'DCTN5', 'DDB1', 'DDOST', 'DDX10',
       'DDX18', 'DDX20', 'DDX21', 'DDX27', 'DDX41', 'DDX47', 'DDX49',
       'DDX55', 'DDX56', 'DGCR8', 'DHODH', 'DHPS', 'DHX15', 'DHX33',
       'DHX37', 'DHX8', 'DHX9', 'UTP25', 'DIMT1', 'DIS3', 'DKC1', 'DLST',
       'DMAP1', 'DNAJA3', 'DNAJC9', 'DNM2', 'DNMT1', 'DOLK', 'DONSON',
       'DPAGT1', 'DTL', 'DTYMK', 'DYNC1I2', 'ECD', 'EEF2', 'EFTUD2',
       'EIF2B1', 'EIF2B3', 'EIF2B5', 'EIF2S1', 'EIF2S2', 'EIF2S3',
       'EIF3A', 'EIF3B', 'EIF3C', 'EIF3D', 'EIF3G', 'EIF3I', 'EIF4A3',
       'EIF5A', 'EIF5B', 'EIF6', 'ELAC2', 'ELL', 'EPRS1', 'ERCC2',
       'ERCC3', 'ERH', 'EXOSC2', 'EXOSC3', 'EXOSC4', 'EXOSC6', 'EXOSC7',
       'EXOSC8', 'CIAO2B', 'FARS2', 'FARSA', 'FARSB', 'FAU', 'FNTA',
       'FNTB', 'FTSJ3', 'GABPA', 'GAPDH', 'GART', 'GEMIN5', 'GEMIN8',
       'GFM1', 'GGPS1', 'GINS2', 'GINS3', 'GINS4', 'GMPPB', 'GMPS',
       'RACK1', 'GNL3', 'GPN3', 'GPS1', 'GRPEL1', 'GRWD1', 'GSPT1',
       'GTF2B', 'GTF2H1', 'GTF2H2C', 'GTF2H4', 'GTF3A', 'GTF3C1',
       'GTF3C2', 'GTF3C5', 'GTPBP4', 'GUK1', 'HARS1', 'HAUS1', 'HAUS5',
       'HCFC1', 'HDAC3', 'HEATR1', 'HINFP', 'H2AC14', 'H2AC18', 'HJURP',
       'HNRNPC', 'HNRNPK', 'HNRNPL', 'HNRNPU', 'HSD17B10', 'HSPA9',
       'HSPD1', 'HUWE1', 'HYPK', 'IARS1', 'IGBP1', 'ILF3', 'IMP3', 'IMP4',
       'INTS1', 'INTS3', 'INTS8', 'INTS9', 'IPO13', 'ISCU', 'ISG20L2',
       'KANSL3', 'KARS1', 'KAT8', 'KIF11', 'KIF23', 'KPNB1', 'KRI1',
       'KRR1', 'LARS1', 'LAS1L', 'LONP1', 'LRR1', 'LSG1', 'LSM11',
       'LSM12', 'LSM2', 'LSM7', 'LUC7L3', 'MAD2L1', 'MAGOH', 'MAK16',
       'MARS1', 'MARS2', 'MASTL', 'MCM3', 'MCM3AP', 'MCM4', 'MCM5',
       'MCM7', 'MDN1', 'MED11', 'MED12', 'MED18', 'MED27', 'MED30',
       'MEPCE', 'METTL16', 'MMS22L', 'MPHOSPH10', 'MRPL57', 'MRPL18',
       'MRPL28', 'MRPL38', 'MRPL4', 'MRPL43', 'MRPL45', 'MRPL46',
       'MRPL53', 'MRPS14', 'MRPS24', 'MRPS34', 'MSTO1', 'MTG2', 'MVK',
       'MYBBP1A', 'MYC', 'NAA10', 'NAA38', 'NAA50', 'NAMPT', 'NAPA',
       'CIAO3', 'NARS1', 'NAT10', 'NCBP1', 'NCBP2', 'NDC80', 'NDUFA13',
       'NEDD8', 'NELFB', 'NHP2', 'SNU13', 'NIP7', 'NKAP', 'NLE1', 'NMD3',
       'NMT1', 'NOC4L', 'NOL10', 'NOL11', 'NOL6', 'NOL9', 'NOP16', 'NOP2',
       'NOP56', 'NOP9', 'NPLOC4', 'NSA2', 'NSF', 'NUDC', 'NUDCD3',
       'NUDT21', 'NUDT4', 'NUF2', 'NUP133', 'NUP155', 'NUP160', 'NUP214',
       'NUP85', 'NUP88', 'NUP93', 'NUS1', 'NUTF2', 'NVL', 'NXF1', 'OGDH',
       'OGT', 'LTO1', 'ORC6', 'OSGEP', 'PABPC1', 'PAFAH1B1', 'PAICS',
       'PAK1IP1', 'PCID2', 'PCNA', 'PFDN2', 'PFN1', 'PGAM1', 'PGGT1B',
       'PGK1', 'PHB', 'PHB2', 'PHF5A', 'PKMYT1', 'PLK1', 'PLRG1', 'PMPCA',
       'PMPCB', 'PNKP', 'POLA2', 'POLR1A', 'POLR1B', 'POLR1C', 'POLR2A',
       'POLR2B', 'POLR2C', 'POLR2D', 'POLR2E', 'POLR2G', 'POLR2H',
       'POLR2I', 'POLR2L', 'POLR3A', 'POLR3C', 'POLR3H', 'POLR3K',
       'POLRMT', 'POP1', 'POP5', 'PPA1', 'PPAN', 'PPAT', 'PPIL2',
       'PPP2CA', 'PTPA', 'PPP4C', 'PPWD1', 'PREB', 'PRELID1', 'PRIM1',
       'PRMT1', 'PRMT5', 'PRPF19', 'PRPF31', 'PRPF38A', 'PRPF38B',
       'PRPF4', 'PRPF8', 'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5',
       'PSMA6', 'PSMA7', 'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB7',
       'PSMC2', 'PSMC3', 'PSMC5', 'PSMC6', 'PSMD1', 'PSMD11', 'PSMD12',
       'PSMD13', 'PSMD14', 'PSMD3', 'PSMD4', 'PSMG3', 'PTPN23', 'PUF60',
       'PWP2', 'EPRS1', 'RABGGTB', 'RACGAP1', 'RAD21', 'RAD51C', 'RAD51D',
       'RAE1', 'RAN', 'RANGAP1', 'RARS2', 'RBBP6', 'RBM14', 'RBM17',
       'RBM8A', 'RBMX', 'RBX1', 'RCC1', 'RCL1', 'RFC2', 'RFC4', 'RFC5',
       'RFK', 'RHEB', 'RIOK2', 'RNF20', 'RNGTT', 'ROMO1', 'RPA1', 'RPA2',
       'RPF2', 'RPL10A', 'RPL11', 'RPL12', 'RPL13', 'RPL14', 'RPL18',
       'RPL18A', 'RPL19', 'RPL23', 'RPL24', 'RPL27', 'RPL27A', 'RPL3',
       'RPL30', 'RPL35', 'RPL35A', 'RPL36', 'RPL37A', 'RPL4', 'RPL6',
       'RPL8', 'RPLP0', 'RPLP1', 'RPLP2', 'RPP21', 'RPP38', 'RPS11',
       'RPS12', 'RPS13', 'RPS15A', 'RPS16', 'RPS18', 'RPS19', 'RPS2',
       'RPS20', 'RPS21', 'RPS23', 'RPS3', 'RPS4X', 'RPS5', 'RPS6', 'RPS7',
       'RPS8', 'RRM1', 'RRP1', 'RRP12', 'RRS1', 'RTCB', 'RUVBL2',
       'SACM1L', 'SAE1', 'SAMM50', 'SAP18', 'SARS1', 'SARS2', 'SART3',
       'SBNO1', 'SDAD1', 'SDHC', 'SEC13', 'SEH1L', 'SF1', 'SF3A2',
       'SF3A3', 'SF3B1', 'SF3B2', 'SF3B3', 'SF3B5', 'SKP1', 'SLC35B1',
       'PRELID3B', 'SLU7', 'SMC1A', 'SMC2', 'SMC4', 'SMU1', 'SNAPC1',
       'SNAPC2', 'SNAPC4', 'SNRNP200', 'SNRNP25', 'SNRNP27', 'SNRNP35',
       'SNRNP70', 'SNRPA1', 'SNRPD1', 'SNRPD2', 'SNRPD3', 'SNRPF', 'SNW1',
       'SPATA5L1', 'SPC24', 'SPC25', 'SRBD1', 'SRP19', 'SRRM1', 'SRRT',
       'SRSF1', 'SRSF2', 'SRSF3', 'SRSF7', 'SS18L2', 'SSU72', 'SUPT5H',
       'SUPT6H', 'SUPV3L1', 'SYMPK', 'SYS1', 'TAF1B', 'TAF6', 'TANGO6',
       'TARS1', 'TBCD', 'TBL3', 'TCP1', 'TELO2', 'TFAM', 'TFRC', 'THOC2',
       'THOC3', 'THOC5', 'TICRR', 'TIMM10', 'TIMM13', 'TIMM23', 'TIMM44',
       'TMEM258', 'TNPO3', 'TOMM22', 'TOMM40', 'TONSL', 'TOP1', 'TOP2A',
       'TPT1', 'TPX2', 'TRAPPC1', 'TRAPPC3', 'TRIAP1', 'TRMT112', 'TRMT5',
       'TRNAU1AP', 'TRRAP', 'TSR1', 'TTC1', 'TTC27', 'TTI1', 'TTI2',
       'TUBB', 'TUBG1', 'TUBGCP2', 'TUBGCP3', 'TUBGCP6', 'TUFM', 'TUT1',
       'TXN', 'TXNL4A', 'U2AF1', 'U2AF2', 'UBA1', 'UBA52', 'UBE2L3',
       'UBE2M', 'UBE2N', 'UBL5', 'UBTF', 'UPF1', 'UPF2', 'UQCRC1',
       'UQCRFS1', 'UROD', 'USP39', 'USP5', 'USPL1', 'UTP15', 'UTP20',
       'UTP23', 'UXT', 'VARS1', 'VARS2', 'VCP', 'VPS25', 'VPS28', 'WARS1',
       'BUD23', 'WDR12', 'WDR25', 'WDR3', 'WDR33', 'WDR43', 'WDR61',
       'WDR70', 'WDR74', 'WDR75', 'WDR77', 'WDR92', 'WEE1', 'XAB2',
       'XPO1', 'XRCC6', 'YARS1', 'YARS2', 'YRDC', 'ZBTB8OS', 'ZMAT5',
       'ZNF131', 'ZPR1', 'ZNF574']

hart_list = pd.Series(hart_list, index=hart_list)

__version__ = 'v1.4.0'

# import logging
# slowvolcLOG = logging.getLogger('slow_volc')
# slowvolcLOG.setLevel(logging.DEBUG)

#v1.3.0
# Adding Mageck tab
# fixed plot volcanos enrichment labeling




#todo; pass plot_volcano a filen string and it loads the table
# todo standardised tables for all results that allow you to pull out x, y wiht same sample and score keys
import pkg_resources

def list_not_str(thing):
    """If thing is a string, wrap it in a list."""
    if type(thing) is str:
        return [thing]
    return thing


def drop_nonumeric(tab):
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )
    return tab

def split_nonumeric(tab):
    """drop columns from DF (tab) where [0, col] type is string,
    and return as a tuple with dropped columns as a sperate DF"""
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    texttab = tab.loc[:, nonumeric].copy()
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )

    return tab, texttab

def size_factor_normalise(cnt_tab, log=True):
    """The number by which MAGeCK uses to normalise its reads.
    Returns factor for each column.

    Uses +1 counts, returns a copy of the table.
    log=True returns np.log2 of the counts."""

    cnt_tab = cnt_tab.copy()
    out_tab = cnt_tab.copy()
    cnt_tab = drop_nonumeric(cnt_tab)
    cnt_tab = cnt_tab + 1

    # cnt_tab = cnt_tab.apply(np.log2)
    # get geometric means, by guide
    gm = stats.gmean(cnt_tab, axis=1)
    # divide the counts by the gmean for each guide
    tab = cnt_tab.T / gm
    # get the median gmean for each experiment and normalise by that
    norm_tab = cnt_tab / tab.T.apply(np.median)
    #norm_tab = norm_tab.apply(round)
    if log:
        norm_tab = norm_tab.apply(np.log2)
    out_tab.loc[:, norm_tab.columns] = norm_tab
    return out_tab


def plot_read_violins(tab, samp_per_row='all', column_labs=None, log=True, size_norm=False, ax=None):
    """Takes counts table, does log2 (optionally size factor normalises)
    and produces violin density plots.

    figsize is (3*n, 6)

    Returns a plt.Axes

    per_row specifies the number of violins on each row"""

    n_samp = tab.shape[1]

    if samp_per_row == 'all':
        samp_per_row = n_samp

    rows = n_samp // samp_per_row + int((n_samp % samp_per_row) > 0)

    if ax is None:
        fig, axes = plt.subplots(rows, 1, figsize=(3 * samp_per_row, 6 * rows))
    else:
        axes = ax

    if rows == 1:
        axes = [axes]

    # remove nonnumeric columns
    tab = tab.copy()
    # drop any text columns
    nonumeric = tab.columns[tab.iloc[0, :].apply(type) == str]
    if len(nonumeric) > 0:
        tab = tab.drop(
            nonumeric,
            axis=1
        )
    # transform the data
    if size_norm:
        tab = size_factor_normalise(tab, log=log)
    elif log:
        tab = tab.__add__(1).apply(np.log2)

    # apply specified column labels
    if column_labs is None:
        column_labs = tab.columns

    for rowi, ax in enumerate(axes):
        inds = slice(rowi * samp_per_row, (1 + rowi) * samp_per_row)
        stab = tab.iloc[:, inds]

        # this is poo
        if stab.shape[0] == 0:
            break
        ax.violinplot(stab, widths=0.8)
        ax.boxplot(stab, widths=0.2)

        ax.set_xticks(range(1, samp_per_row + 1))
        clabs = list(column_labs[inds])
        while len(clabs) < samp_per_row:
            clabs.append('')
        ax.set_xticklabels(clabs, rotation=40)
    plt.tight_layout()
    if rows == 1:
        return axes[0]
    else:
        return axes

def get_ROC_values(series, things_oi):
    things_oi = things_oi[things_oi.isin(series.index)]

    col = series.dropna().sort_values()
    toi = things_oi[things_oi.isin(col.index)]
    y = np.cumsum(col.index.isin(toi)) / len(toi)
    x = np.arange(col.shape[0])/col.shape[0]
    return x, y

def plot_ROC(tab, things_oi, label=None, ax = None):
    #todo: make things_oi work as any kind of iterable (currently pd.Series)
    """Tab needs to be a series or dataframe containing values used to
    order the tab index."""

    if type(things_oi) is list:
        things_oi = pd.Series(things_oi, index=things_oi)

    if ax is None:
        _, ax = plt.subplots(1,1, figsize=(4,4))
    else:
        plt.sca(ax)
    tab = tab.copy()
    things_oi = things_oi[things_oi.isin(tab.index)]
    #print(things_oi)
    if len(tab.shape) > 1:
        for i in range(tab.shape[1]):
            if label is not None:
                lab = label[i]
            else:
                lab = tab.columns[i]
            # get ordered list of the values and subset things of interest to exclude nan values.
            col = tab.iloc[:, i].dropna().sort_values()
            toi = things_oi[things_oi.isin(col.index)]
            plt.plot(np.arange(col.shape[0])/col.shape[0],
                     np.cumsum(col.index.isin(toi))/len(toi),
                    label=lab)
    else:
        tab = tab.sort_values()
        plt.plot(np.arange(tab.shape[0])/tab.shape[0], np.cumsum(tab.index.isin(things_oi))/len(things_oi))
#     plt.plot(np.arange(kin.shape[0])/kin.shape[0], np.cumsum(kin.index.isin(gn))/len(gn),
#             label='RPE1 screen essentials')
    plt.plot([0,1], [0,1], 'k--', alpha=0.4)
    plt.xlabel('Not essential')
    plt.ylabel('Essential')
    plt.title('ROC')
    plt.legend()
    plt.tight_layout()
    return ax


def plot_volcano_from_mageck(tab, title='', label_genes=None, outfn='', ax=None, source='mageck'):
    """Take a mageck table and do a volcano plot. """

    dep = tab.loc[:, 'pos|lfc'] > 0
    enr = tab.loc[:, 'pos|lfc'] < 0
    tab.loc[dep, 'ml10_fdr'] = -np.log10(tab.loc[dep, 'pos|fdr'])
    tab.loc[enr, 'ml10_fdr'] = -np.log10(tab.loc[enr, 'neg|fdr'])
    plot_volcano(tab['pos|lfc'], tab['ml10_fdr'], title=title, other_labels=label_genes,
                 outfn=outfn, ax=ax)

#POINT = 0

def plot_volcano(lfc, fdr, tab=None, title='', label_deplet=0, label_enrich=0,
                 other_labels=None, p_thresh=0.05, outfn='', ax=None,
                 exclude_labs=('NonT', 'Rando'), plot_kw: dict = None):
    """Draw a volcano plot of lfc vs fdr. assumes fdr is -log10.

    :param lfc: str giving tab[lfc] or series with gene names as index
    :param fdr: str giving tab[fdr] or series with gene names as index
    :param tab: DataFrame containing lfc and fdr
    :param title: optional plot title
    :param label_enrich: int giving top n enriched genes to label
    :param label_deplet: int giving top n depleted genes to label
    :param other_labels: other genes to label
    :param p_thresh:
    :param outfn: if provided a .png will be written
    :param ax: optional plt.Axes instance to use
    :return: plt.Axes
    """
    # print('TP ', POINT)

    # this is silly
    lfc_lab, fdr_lab = None, None
    if tab is not None:
        lfc_lab = lfc
        fdr_lab = fdr
        lfc = tab[lfc]
        fdr = tab[fdr]

    sctkw = dict(marker='o', linestyle='none', alpha=0.4)
    if plot_kw is not None:
        sctkw.update(plot_kw)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 12))

    # todo update to use scatter and set minimum matplotlib version
    ax.plot(lfc, fdr, **sctkw)
    # plt.yscale('log')
    # plt.gca().invert_yaxis()
    ax.set_title(title)

    xmin, xmax, ymin, ymax = ax.axis()

    p_thresh = -np.log10(p_thresh)
    ax.plot([xmin, xmax], [p_thresh, p_thresh], 'k--')
    ax.plot([0, 0], [ymin, ymax], color='silver')
    # label the top and bottom most, if specified
    texts = []
    texts_done = []

    # get subtables
    # filter out excluded labels by getting Series containing only included labs,
    # turn that into a mask of the original series and combining that with the dep/enr masks.
    filtered_lfc = lfc.copy()
    if exclude_labs:
        for exclude in exclude_labs:
            filtered_lfc = filtered_lfc.loc[~filtered_lfc.index.str.contains(exclude)]

        included = lfc.index.isin(filtered_lfc.index)
        depmask = (lfc < 0) & included
        enrmask = (lfc > 0) & included
    else:
        depmask = (lfc < 0)
        enrmask = (lfc > 0)

    # get tails from ascending order
    if tab is not None:
        dep = tab.loc[depmask, :].sort_values(lfc_lab, ascending=False).sort_values(fdr_lab).tail(label_deplet)[fdr_lab]
        enr = tab.loc[enrmask, :].sort_values([fdr_lab, lfc_lab]).tail(label_enrich)[fdr_lab]
        # label the tails
        for end in dep, enr:
            for lab, an_fdr in end.items():
                if an_fdr < p_thresh:
                    continue
                texts_done.append(lab)
                texts.append(plt.text(lfc[lab], fdr[lab], lab))
    # label additional genes
    if other_labels is not None:
        for lab in other_labels:
            if lab in texts_done:
                continue
            texts.append(
                ax.text(lfc[lab], fdr[lab], lab)
            )

    if texts and adjust_text:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'), ax=ax)

    if lfc_lab is None:
        ax.set_xlabel('Log$_2$ Fold Change')
    else:
        ax.set_xlabel(lfc_lab)

    ax.set_ylabel('-log$_{10}$(FDR)')

    if outfn:
        plt.tight_layout()
        plt.savefig(outfn, dpi=150)
    return ax


revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])

def iter_files_by_prefix(prefix:Union[str, Path], req_suffix=None):
    prefix = Path(prefix)
    check_suffix = lambda s: s.endswith(req_suffix) if req_suffix is not None else True

    for fn in os.listdir(prefix.parent):
        # filter incorrect files, the '._' are mac files that are not ignored automatically on unix filesystem
        if not check_suffix(fn) or \
                prefix.parts[-1] not in fn or \
                fn.startswith('._'):
            continue
        yield fn

def tabulate_mageck(prefix):
    """
    :param prefix: Input file prefix, including path
    :return: pd.DataFrame
    """
    prefix = Path(prefix)
    tables = {}
    tab = None
    for fn in iter_files_by_prefix(prefix, '.gene_summary.txt'):

        #mtab from the mageck output, reformatted into tab
        mtab = pd.read_csv(prefix.parent / fn, sep='\t', index_col=0)
        tab = pd.DataFrame(index=mtab.index)

        # keep the pos and neg p values with better names
        tab.loc[:, 'neg_p'] = mtab['neg|p-value']
        tab.loc[:, 'pos_p'] = mtab['pos|p-value']
        tab.loc[:, 'neg_fdr'] = mtab['neg|fdr']
        tab.loc[:, 'pos_fdr'] = mtab['pos|fdr']

        tab.loc[:, 'lfc'] = mtab.loc[:, 'neg|lfc']
        # turn sep pos|neg columns into one giving only the appropriate LFC/FDR
        pos = mtab['pos|lfc'] > 0
        for stat in 'fdr', 'p-value':
            tabk = stat.replace('-value', '')
            tab.loc[pos, tabk] = mtab.loc[pos, f'pos|{stat}']
            tab.loc[~pos, tabk] = mtab.loc[~pos, f'neg|{stat}']
            tab.loc[:, f'{tabk}_log10'] = tab[tabk].apply(lambda x: -np.log10(x))

        sampnm = fn.split(prefix.stem)[1].split('.gene_s')[0]
        tables[sampnm] = tab
    if tab is None:
        raise FileNotFoundError('Failed to find any .gene_summary.txt files with prefix '+str(prefix))

    # create multiindex using the sample names and stats
    tbcolumns = pd.MultiIndex.from_product(
        [sorted(tables.keys()), ['lfc', 'fdr', 'fdr_log10', 'p', 'p_log10', 'pos_p', 'neg_p',
                                 'neg_fdr', 'pos_fdr']],
        1
    )
    table = pd.DataFrame(index=tab.index, columns=tbcolumns)
    for exp, tab in tables.items():
        # deal with missing genes
        m = tab.isna().any(1)
        tab.loc[m, :] = 1
        tab.loc[m, ['lfc', 'fdr_log10', 'p_log10']] = 0
        table[exp] = tab



    return table

def tabulate_drugz(prefix, compjoiner='→'):
    prefix=Path(prefix)
    tables = {}
    for fn in iter_files_by_prefix(prefix):
        comp = fn[len(prefix.parts[-1]):-4].replace('-', compjoiner)

        tab = pd.read_csv(os.path.join(os.path.split(prefix)[0], fn), sep='\t', index_col=0)

        # sort out the column names to be consistent with other results tables
        tab.index.name = 'gene'
        tab = tab.loc[:, ['normZ', 'pval_synth', 'fdr_synth', 'pval_supp', 'fdr_supp']]
        stats_cols = 'normZ neg_p neg_fdr pos_p pos_fdr'.split()

        tab.columns = stats_cols

        # get the minimum value significance stats
        for stat in 'p', 'fdr':
            min_stat = tab.loc[:, [f'neg_{stat}', f'pos_{stat}']].min(1)
            tab.loc[:, f'{stat}'] = min_stat
            tab.loc[:, f'{stat}_log10'] = min_stat.apply(lambda p: -np.log10(p))
        tables[comp] = tab

    tbcolumns = pd.MultiIndex.from_product(
        [sorted(tables.keys()), tab.columns],
        1
    )
    table = pd.DataFrame(index=tab.index, columns=tbcolumns)
    for exp, tab in tables.items():
        table[exp] = tab
    return table


def pca_grid(pca, hue_deet, style_deet, max_components=5, also_return_fig=False):
    """Plot components against each other in a grid of scatter plots.
    Deets MUST be in the same order as the columns used for the PCA.
    Args:
        pca: a sklearn.decomposition.PCA() object that has been fit
        hue_deet: grouping variables used for setting colours, by sample
        style_deet: as hue_deet but marker style
        max_components: components to plot, max_comp*(max_comp-1) scatter plots will be produced
        also_return_fig: When false only axes are returned, set to true to also return the Figures"""

    thing = sns.scatterplot(pca.components_[0], pca.components_[0],
                            hue=hue_deet, style=style_deet,
                            s=150)

    leg = thing.get_legend_handles_labels()
    plt.close()

    max_b = max_components - 1
    fig, axes = plt.subplots(max_components, max_b, figsize=(4.5 * max_b, 3.85 * max_components))
    # pc ind also used for subplots
    for pc_a in range(max_components):
        for pc_b in range(max_b):
            plt.sca(axes[pc_a][pc_b])

            if pc_a == pc_b:
                plt.legend(*leg)
                continue

            # pc_a/b swapped to get column PC on the x and row PC on the y
            sns.scatterplot(pca.components_[pc_b], pca.components_[pc_a],
                            hue=hue_deet, style=style_deet,
                            s=150, legend=False)
            plt.xlabel(f"PC {pc_b+1} (variance explained: {pca.explained_variance_ratio_[pc_b]*100:.3}%)")
            plt.ylabel(f"PC {pc_a+1} (variance explained: {pca.explained_variance_ratio_[pc_a]*100:.3}%)")
    plt.tight_layout()
    if also_return_fig:
        return (fig, axes)
    else:
        return axes


def write_repmap(sample_reps:Dict[str, list], ctrlmap:Dict[str, list], repmap_fn:str):
    """Write a JACKS format replicate file.

    Args:
        sample_reps: maps given sample names to replicate names in the count file.
        ctrlmap: specifies which samples are controls for which other samples.
        repmap_fn: path to file that will be written.
    Returns:
        None
        """
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


def get_req_infection(libcnt, minimum=100, quantile = 0.05):
    libcnt = libcnt.sort_values().copy()
    x = libcnt.quantile(quantile)
    if x == 0:
        return np.inf
    return sum(libcnt)/(x/minimum)

def plot_req_inf(counts, reps, qrange=(0.01,0.1), moi=0.2):
    plt.figure(figsize=(6, 10))
    for n in reps:
        ys=[]
        xs = np.linspace(*qrange, 100)
        for quantile in xs:
            ys.append(get_req_infection(counts, n, quantile = quantile))
        ys=pd.Series(ys)
        plt.plot((1-xs)*100, ys/moi, label=str(n))
        plt.xlabel("% guides > X ")
        plt.ylabel("required cell infections")
    plt.legend(title='X (guide rep)')


def load_analyses_via_expd(expd_or_yaml:Union[Dict, PathType],
                           results_types:Union[str, List[str]]='infer',
                           use_attrdict=True) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Return results from tables generated by crispr_pipeline.py.

    results_types is inferred from the yaml by default, or a single string or
    list of str giving results types can be given.

    Returned results will be keyed by results_type then control group"""

    if type(expd_or_yaml) != dict:
        with open(expd_or_yaml) as f:
            expd_or_yaml = yaml.safe_load(f)

    expd = expd_or_yaml

    #todo support DrugZ.
    #  analyses should be classes that are global and have tabulate method attached
    supported_analyses = ['jacks', 'mageck']
    # TODO consistent file name for results tables
    suffixes = dict(zip(supported_analyses, ['jacks_scores', 'mageck_table']))

    if results_types == 'infer':
        results_types = []
        for rt in supported_analyses:
            if not expd[f"skip_{rt}"]:
                results_types.append(rt)

    LocalDict = dict
    if use_attrdict:
        LocalDict = AttrDict
    all_results = LocalDict()


    for rt in results_types:
        all_results[rt] = res = LocalDict()
        for ctrl_grp in expd['controls'].keys():
            fn = f"{expd['exp_name']}/{expd['analysis_name']}/{rt}/tables/{expd['file_prefix']}.{ctrl_grp}.{suffixes[rt]}.csv"
            res[ctrl_grp] = pd.read_csv(fn, index_col=0, header=[0, 1])
    if use_attrdict:
        return AttrDict(all_results)
    return all_results


def select_clonal_lfc_by_samp(control_samp, treat_samp, lfcs, sample_reps,):
    #todo make a ClonalLFC object with this as a method

    ctrl_reps  = sample_reps[control_samp]
    treat_reps = sample_reps[treat_samp]

    # get bool masks of comps with the control reps, and treat reps
    ctrl_mask, treat_mask = [lfcs.columns.map(lambda s: s.split('-')[i] in reps)
                             for i, reps in ((0, ctrl_reps), (1, treat_reps))]

    return lfcs.loc[:, ctrl_mask.values & treat_mask.values]


def get_clonal_lfcs(lncounts, ctrl_dict:dict, sample_reps:dict,
                    lognorm=False):
    """get a DF of clonal LFCs using the ctrl/sample pairs specified by ctrl_dict.
    Assumes that clones are grouped by order of appearance in sample_reps"""
    _lfcs = {}

    if lognorm:
        lncounts = size_factor_normalise(lncounts)

    for ctrl_samp, treat_samples in ctrl_dict.items():
        treat_samples = list_not_str(treat_samples)
        for trt_samp in treat_samples:
            # we shall assume that the clones are int he same order in sample_reps
            c_reps, t_reps = sample_reps[ctrl_samp], sample_reps[trt_samp]
            for A, B in (c_reps, t_reps), (t_reps, c_reps):
                while len(A) < len(B):
                    A.append(A[0])

            treat_clone_pairs = zip_longest(B, A)
            for c, t in treat_clone_pairs:

                _lfcs[c+'-'+t] = lncounts[t] - lncounts[c]

    return pd.DataFrame(_lfcs)


### WARNING: This completely fucks things up sometimes, do not uncomment unless you're debugging
# def write_results_excel(results:Dict[str, pd.DataFrame],
#                         filename,
#                         score_name='',
#                         stat_labels:Dict[str, str]=None):
#     """Write the results from a single analysis to an Excel file.
#     By default columns fdr, fdr_log10, neg_p, pos_p are included.
#     Args:
#         results: A dict of DataFrames. The keys will be used to name the sheets
#         filename: name of the file that will be written
#         score_name: 'lfc' or 'jacks_score' depending on what was used. Use
#             stat_labels for anything else.
#         stat_labels: Custom columns can be included in the output by passing a
#             dictionary with dataframe_label:output_label mapped. You'll
#             need to include every stat of interest.
#     """
#
#     import xlsxwriter
#
#     if not score_name and stat_labels is None:
#         RuntimeWarning('No score key included')
#
#     if not filename.endswith('.xlsx'):
#         filename += '.xlsx'
#     workbook = xlsxwriter.Workbook(filename)
#
#     for sheetname, tab in results.items():
#         tab = tab.copy()
#         # ** rename the stats columns **
#         #    this could be in a separate function...
#
#         tab.index.name = 'Gene'
#         # get the score label for common scores
#         score_lab = {'lfc':'Log2(FC)', 'jacks_score':'JACKS score', '':''}[score_name]
#
#         if stat_labels:
#             good_stats = stat_labels
#         else:
#             good_stats = dict(
#                 zip(['fdr', 'fdr_log10', score_name, 'neg_p', 'pos_p'],
#                     ['FDR', '–Log10(FDR)', score_lab, 'Dropout p-value', 'Enrichment p-value'])
#             )
#
#         # orginally wrote this to work with multiindex columns, but also
#         if hasattr(tab.columns, 'levels'):
#             # get the new name for favoured stats, and record those to be removed
#             new_level = []
#             dump_stats = []
#             for c in tab.columns.levels[1]:
#                 try:
#                     new_level.append(good_stats[c])
#                 except KeyError:
#                     new_level.append(c)
#                     dump_stats.append(c)
#
#             # change the labels of the ones we care for
#             tab.columns.set_levels(new_level, 1, inplace=True)
#             # drop the values from the table, does not effect the index
#             tab.drop(dump_stats, 1, level=1, inplace=True)
#
#             # generate a new multiindex
#             new_multiindex = []
#             for col in tab.columns:
#                 if col[1] not in dump_stats:
#                     new_multiindex.append(col)
#
#             tab.columns = pd.MultiIndex.from_tuples(new_multiindex)
#
#
#         # ** Write the worksheet **
#         sheet = workbook.add_worksheet(name=sheetname)
#
#         index_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'right'})
#         header_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'center'})
#         num_format = workbook.add_format({'num_format': '0.00'})
#
#         # do the cell merge, we're gonna use the length of levels[1] to figure out the merge range
#         n_to_merge = len(tab.columns.levels[1])
#         samples = tab.columns.levels[0]
#
#         for samp_i, samp in enumerate(samples):
#             # row, col, row, col
#             sheet.merge_range(0, 1 + samp_i * n_to_merge, 0, 0 + (samp_i + 1) * n_to_merge, samp, header_format)
#
#         # write the index header
#         sheet.write(0, 0, 'Gene')
#
#         # write the stats headers
#         for coli, (c1, c2) in enumerate(tab.columns):
#             # c1 gets wrote using merge, above
#             sheet.write(1, coli + 1, c2, index_format)
#
#         for row_i, gene in enumerate(tab.index):
#             sheet.write(2 + row_i, 0, gene, index_format)
#
#         for col_i, col in enumerate(tab):
#             for row_i, val in enumerate(tab[col]):
#                 if not np.isnan(val):
#                     sheet.write(row_i + 2, col_i + 1, val, num_format)
#
#     workbook.close()


def write_stats_workbook(sheets: Dict[str, pd.DataFrame], filename=None,
                         workbook:xlsxwriter.Workbook=None,
                         close_workbook=True) -> xlsxwriter.Workbook:
    """First row and column bold text, everything else numbers with 3 d.p.

    Args:
        sheets: dict containing dataframes to be written, keys used as sheet names.
        filename: The path to which the file will be written. Not required if passing
            an open workbook
        workbook: An xlsxwriter.Workbook can be passed
        close_workbook: Close and write the workbook if True.
    """

    if workbook is None:
        if filename is None and close_workbook:
            raise RuntimeError('Provide a filename, or set close_workbook=False to avoid an Error')
        workbook = xlsxwriter.Workbook(filename)
    index_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'right'})
    header_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'center'})
    num_format = workbook.add_format({'num_format': '0.000'})

    for sheet_name, tab in sheets.items():

        sheet = workbook.add_worksheet(name=sheet_name)

        # write the index header
        sheet.write(0, 0, 'Gene')

        # write the stats headers
        for coli, c in enumerate(tab.columns):
            # print(c)
            sheet.write(0, coli + 1, c, header_format)

        for row_i, gene in enumerate(tab.index):
            sheet.write(1 + row_i, 0, gene, index_format)

        for col_i, col in enumerate(tab):
            for row_i, val in enumerate(tab[col]):
                sheet.write(row_i + 1, col_i + 1, val, num_format)
    if close_workbook:
        workbook.close()
    return workbook

def quantile_normalise(df):
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    return df.rank(method='min').stack().astype(int).map(rank_mean).unstack()



# function that takes two samples and returns ps based on control mask
def ps_from_kde(x, y, ctrl_mask, table=None, score_key='jacks_score',
                delta_lobf=False):
    if table is not None:
        x, y = [table[k][score_key] for k in (x, y)]

    if delta_lobf:
        # get the delta of y from y predicted by lobf
        fit = OLS(y, sm.add_constant(x)).fit()
        pred_y = fit.predict(sm.add_constant(x))
        diff = y - pred_y
    else:
        diff = y - x

    ctrls = diff.loc[ctrl_mask]
    kde_dist = stats.gaussian_kde(ctrls)
    dep_ps = diff.loc[~ctrl_mask].apply(
        lambda x: kde_dist.integrate_box_1d(-100, x)
    )
    enr_ps = 1 - dep_ps
    dep_fdr, enr_fdr = [sm.stats.multipletests(ps, method='fdr_bh', )[1]
                        for ps in (dep_ps, enr_ps)]
    return pd.DataFrame({"enriched_p": enr_ps,
                         "depleted_p": dep_ps,
                         "enriched_fdr": enr_fdr,
                         "depleted_fdr": dep_fdr}, index=diff.index)


def scores_scatter_plotly(x, y, fdr,
                          fdr_thresholds=(0.05,),
                          fdr_colours=('#CC7700',),
                          gene_set_masks: List[Tuple] = None, ):
    """gene_set_masks, [('label1', mask1), ('label2', mask2), ...]. Mask here either a
    pd.Series boolean mask, or list of genes that can be passed to .loc[]"""
    from plotly import graph_objects as go
    fig = go.Figure()
    fig.update_layout(template='plotly_white')
    Xall = x
    Yall = y

    # These will be populated with values and passed to go.Scatter
    xys = []
    marker_labels = []
    markers = []
    symbols = []

    # get list of genes between FDR thresholds
    #  Sort thresholds and colours, largest to smallest
    sranks = stats.rankdata(fdr_thresholds, 'ordinal')
    fdr_thresholds = sorted(fdr_thresholds, key=lambda x: sranks[fdr_thresholds.index(x)], reverse=True)
    fdr_colours = sorted(fdr_colours, key=lambda x: sranks[fdr_colours.index(x)], reverse=True)

    # rank big to small
    sranks = stats.rankdata(fdr_thresholds, 'ordinal')
    fdr_thresholds = [1] + sorted(fdr_thresholds, key=lambda x: sranks[fdr_thresholds.index(x)], reverse=True)
    fdr_colours = ['#d9d9d9'] + sorted(fdr_colours, key=lambda x: sranks[fdr_colours.index(x)], reverse=True)

    # get masks giving a gene's fdr bracket
    sig_masks = []
    genes_in_mask = []
    for i in range(len(fdr_thresholds) - 1):
        m = (fdr_thresholds[i + 1] < fdr) & (fdr <= fdr_thresholds[i])
        m = m[m].index
        sig_masks.append(m)
        genes_in_mask.extend(m)
    m = fdr < fdr_thresholds[-1]
    sig_masks.append(m[m].index)

    # populate the fdr related Scatter values
    for m, c in zip(sig_masks, fdr_colours):
        xys.append((Xall.loc[m], Yall.loc[m]))
        symbols.append('circle-open')
        markers.append(dict(color=c, size=7))

    marker_labels.append(f"FDR > {fdr_thresholds[1]}")
    for thresh in fdr_thresholds[1:]:
        marker_labels.append(f"FDR < {thresh}")
    # finished with FDR

    # populate gene set scatter values
    if gene_set_masks is not None:
        for lab, m in gene_set_masks:
            xys.append((Xall.loc[m], Yall.loc[m]))
            symbols.append('circle')
            markers.append(dict(size=4))
            marker_labels.append(lab)

    # plot it finally
    for (x, y), data_label, marker, symbol in zip(xys, marker_labels, markers, symbols):
        fig.add_trace(
            go.Scatter(
                x=x, y=y,
                name=data_label,
                mode='markers',
                marker_symbol=symbol,
                marker=marker,
                text=x.index,
                customdata=fdr.loc[x.index],
                hovertemplate='%{text}:<br>  %{x:.2f},%{y:.2f}<br>  FDR=%{customdata:.2f}',

            )
        )

    return fig

def write_plotly_html(results_table:pd.DataFrame,
                      xy_key:str,
                      samplePairs:Iterable=None,
                      gene_set_masks:Tuple[str, pd.Series]=None,
                      fdr_thresholds=(0.05,),
                      fdr_colours=('#CC7700',),
                      out_fmt_str='',
                      samp_labels=Dict[str, str],
                      samp_pair_key_str='{}-{}',
                      show_fig=False):
    """
    Produces set of HTMLs applying scores_scatter_plotly using a table, but I
    don't really get the format of the input table and need to go back and look
    at where I originally used this.
    TODO: how do you use write_plotly_html

    Args:
        results_table: DF with multiindex columns, (sample, stat). Must include
            neg_fdr & pos_fdr stat columns
        xy_key: The stat name (in results_table) that will be used as x/y values
        samplePairs: All X/Y samples that will have html produced, all pairs by
            default.
        gene_set_masks: List of labelling strings and bool masks that indicate
            subsets of genes that will be plotted separately as small dots within
            the main points.
        out_fmt_str: Path to which files will be written, containing one {} to
            which the comparison will be written.
        samp_labels: Dictionary of sample labels
        samp_pair_key_str:


    """

    if samplePairs is None:
        if type(results_table) is pd.DataFrame:
            samps = results_table.columns.levels[0]
        else:
            samps = results_table.keys()
        samplePairs = itertools.permutations(samps, 2)

    for rsamp, psamp in samplePairs:

        smpk = samp_pair_key_str.format(rsamp, psamp)

        Xall, Yall = results_table[smpk][xy_key]

        enr, dep = [results_table[smpk][pk] for pk in ('pos_fdr', 'neg_fdr')]

        # get the lowest fdr for hovertext
        fdr_lowest = enr.copy()
        # wherever dep is smaller, overwrite the enr fdr value
        m = dep < enr
        fdr_lowest[m] = dep[m]

        axis_labels = [samp_labels[s] if samp_labels else s for s in (rsamp, psamp)]

        fig = scores_scatter_plotly(Xall, Yall, fdr_lowest, fdr_thresholds, fdr_colours, gene_set_masks)

        # comp = f"{rsamp}-{psamp}"
        if out_fmt_str:
            fn = out_fmt_str.format(rsamp, psamp)
            fig.write_html(
                fn,
                include_plotlyjs='directory',
            )
        if show_fig:
            fig.show()


def plot_rank_vs_score(score:pd.Series, sig:pd.Series, n_labels, sig_threshold=0.05, step=0.015):
    """Score is the thing they will be ranked on,
    sig series is just used for drawing the significance line
    todo: make sig optional."""

    def get_rank_of_thresholds(score, sig, threshold=0.05):

        # first dropouts, then enrichments
        threshold_yvalue = []
        for endmask in (score < 0), (score > 0):
            yend = sig.loc[endmask]
            thresh_rank = (yend < threshold).sum()
            threshold_yvalue.append(
                score.loc[endmask].abs().sort_values(ascending=False).iloc[thresh_rank]
            )

        threshold_yvalue[0] = -threshold_yvalue[0]
        return threshold_yvalue

    plt.figure(figsize=(7, 10))
    #tab = res_drugz[k]
    normz_rank = score.rank()
    plt.scatter(normz_rank, score)
    low, hi = get_rank_of_thresholds(score, sig, sig_threshold)
    for thresh in (low, hi):
        plt.plot([0, score.shape[0]], [thresh, thresh], 'k--', alpha=0.6)
        plt.text(score.shape[0] / 2, thresh + 0.1, f'{sig_threshold * 100}% FDR')

    for isneg, endmask in (True, (score < 0)), (False, (score > 0)):
        ax = plt.gca()
        axis_to_data = (ax.transAxes + ax.transData.inverted()).transform
        scoreend = score.loc[endmask]

        top_genes = scoreend.loc[(sig.loc[endmask] < sig_threshold)].abs().sort_values(ascending=False).head(n_labels).index
        try:
            extreme_value = scoreend.loc[top_genes[0]]
        except IndexError:
            continue

        for gi, gn in enumerate(top_genes):
            x = normz_rank[gn]
            y = score[gn]


            if isneg:
                text_x, text_y = axis_to_data([0.25, 0.05 + gi * step])
                alignment = 'left'

            else:
                text_x, text_y = axis_to_data([0.75, 1 - (gi * step)])
                alignment = 'right'

            plt.text(
                text_x, text_y, gn,
                fontsize=9,
                bbox={'facecolor': 'white', 'alpha': 1, 'pad': 0.1},
                ha=alignment,
            )
            plt.plot([x + 50, text_x, ], [y, text_y], 'k-')

def load_counts(file_path, index_col='guide', gene_col='gene') -> (pd.DataFrame, pd.Series):
    """Load count table, putting gene column into separate Series and dropping it
    from counts. Returns (cnt, guide_gene)
    """
    cnt = pd.read_csv(file_path, sep='\t', index_col=index_col)
    guide_gene = cnt[gene_col]
    cnt.drop(gene_col, 1, inplace=True)
    return cnt, guide_gene

def load_replicate_sample_details(xlsx) -> (pd.DataFrame, pd.DataFrame):
    """From a screen details Excel file (sheet='Sample details), return DFs
    containg replicate and sample details"""
    rep_deets = pd.read_excel(xlsx, sheet_name='Sample details', )
    rep_deets = rep_deets.set_index('Sequence name', drop=False)
    samp_deets = rep_deets.drop_duplicates('Sample').set_index('Sample', drop=False).drop(['Replicate', 'Sequence name'], 1)
    return rep_deets, samp_deets