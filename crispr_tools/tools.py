import os, itertools
import numpy as np
from scipy import stats
import re


from typing import List, Union, Dict, Collection, Iterable, Tuple
PathType = Union[str, bytes, os.PathLike]
from pathlib import Path, PosixPath, WindowsPath
import pandas as pd
from itertools import zip_longest
from attrdictionary import AttrDict

import statsmodels.api as sm
OLS = sm.regression.linear_model.OLS
import xlsxwriter

from crispr_tools.data_classes import (
    analysis_score_names,
    statistic_labels
)
ARROW = '→'
#import multiprocessing as mp

def maybe_its_gz(filename) -> str:
    """If filename doesn't exist, look for filename+'.gz'. If the gz version exists,
    return that filename, otherwise return the original filename."""
    filename = str(filename)
    if not os.path.isfile(filename):
        filegz = filename+'.gz'
        if os.path.isfile(filegz):
            return filegz
        else:
            return filename
    else:
        return filename

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

# import logging
# slowvolcLOG = logging.getLogger('slow_volc')
# slowvolcLOG.setLevel(logging.DEBUG)

#v1.3.0
# Adding Mageck tab
# fixed plot volcanos enrichment labeling

def is_olfactory(gene_symbol):
    m = re.match('^OR\\d+[a-zA-Z]+\\d+$', gene_symbol)
    if m is None:
        return False
    return True


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


def abundance_normalise(c:pd.DataFrame, pseudocount=1, log2=True, drop_cols=('gene',)):
    """Returns table with column counts of equal sum.

    Pseudocount of zero will be a problem.
    """
    if drop_cols:
        for col in drop_cols:
            try:
                c = c.drop(col, 1)
            except:
                pass

    lnc = (c/c.sum())*c.sum().median() + pseudocount
    if log2:
        lnc = lnc.apply(np.log2)
    lnc.head()
    return lnc

def ROC_values(values:pd.Series, things_oi:pd.Series):
    """
    :param values: Series, index with labels, values used to do the ROC
    :param things_oi: Labels found in values.index that are "true"
    :return: x y values for plot, y being proportion of ordered value labels
        found in things_oi
    """
    if type(things_oi) in (np.ndarray, list, tuple):
        things_oi = pd.Series(things_oi, index=things_oi)
    col = values.dropna().sort_values()
    toi = things_oi[things_oi.isin(col.index)]
    y = np.cumsum(col.index.isin(toi)) / len(toi)
    x = np.arange(col.shape[0])/col.shape[0]
    return x, y


revcomp = lambda s: ''.join([dict(zip('ACTGN', 'TGACN'))[nt] for nt in s[::-1]])

def iter_files_by_prefix(prefix:Union[str, Path], req_suffix=None, allow_gz=True):
    prefix = Path(prefix)
    if req_suffix is not None:
        if allow_gz:
            check_suffix = lambda s: (s.endswith(req_suffix) or s.endswith(req_suffix+'.gz'))
        else:
            check_suffix = lambda s: s.endswith(req_suffix)
    else:
        check_suffix = lambda s: True

    for fn in os.listdir(prefix.parent):
        if (
            check_suffix(fn)  and
            (prefix.parts[-1] in fn) and
            (not fn.startswith('._')) # mac files not ignored by unix
        ):
            yield fn

def tabulate_mageck(prefix, compjoiner=ARROW) -> dict[str, pd.DataFrame]:
    """One DF per comparison.

    Args:
        prefix: full path including the string added to output files.
        compjoiner: symbol used to join comparisons, that forms part of the filename."""
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

        sampnm = fn.split(prefix.stem, 1)[1].split('.gene_s')[0].replace('-', compjoiner).replace('.', '')
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
        m = tab.isna().any(axis=1)
        tab.loc[m, :] = 1
        tab.loc[m, ['lfc', 'fdr_log10', 'p_log10']] = 0
        table[exp] = tab

    return table

def tabulate_drugz(prefix, compjoiner=ARROW) -> dict[str, pd.DataFrame]:
    """One DF per comparison.

    Args:
        prefix: full path including the string added to output files.
        compjoiner: symbol used to join comparisons, that forms part of the filename.
            Will be replaced by hyphen. If it's already hyphen it doesn't matter."""
    prefix=Path(prefix)
    tables = {}
    for fn in iter_files_by_prefix(prefix):
        # crop the fn removing the filename part of the prefix, and the file extension
        #  this is the comparison
        comp = fn.split(prefix.stem, 1)[1]
        comp = re.sub("\\.tsv.*$", "", comp)
        comp = comp.replace("-", compjoiner)

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



def load_analyses_via_expd(expd_or_yaml:Union[Dict, PathType],
                           results_types:Union[str, List[str]]='infer',
                           use_attrdict=True) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Return results from tables generated by crispr_pipeline.py.

    results_types is inferred from the yaml by default, or a single string or
    list of str giving results types can be given.

    Returned results will be keyed by results_type then control group"""

    import yaml
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

    ctrl_reps  = sample_reps[control_samp]
    treat_reps = sample_reps[treat_samp]

    # get bool masks of comps with the control reps, and treat reps
    ctrl_mask, treat_mask = [lfcs.columns.map(lambda s: s.split('-')[i] in reps)
                             for i, reps in ((0, ctrl_reps), (1, treat_reps))]

    return lfcs.loc[:, ctrl_mask.values & treat_mask.values]


def clonal_lfcs(lncounts, ctrl_dict:dict, sample_reps:dict,
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


def write_stats_workbook(sheets: Dict[str, pd.DataFrame], filename=None,
                         num_format='0.000',
                         workbook:xlsxwriter.Workbook=None,
                         close_workbook=True,
                         freeze_indexes:Union[bool, dict]=True,) -> xlsxwriter.Workbook:
    """First row and column bold text, everything else numbers with 3 d.p.

    Args:
        sheets: dict containing dataframes to be written, keys used as sheet names.
        filename: The path to which the file will be written. Not required if passing
            an open workbook
        num_format: Format of numbers, defaults to 3 decimal places.
        workbook: An xlsxwriter.Workbook instance (optional).
        close_workbook: Close and write the workbook if True.
        freeze_indexes: Freeze the first column and row in each sheet. Set True to freeze
            all, or pass a dict with sheetName->bool to specify sheets.
    """

    if workbook is None:
        if filename is None and close_workbook:
            raise RuntimeError('Provide a filename, or set close_workbook=False to avoid an Error')
        workbook = xlsxwriter.Workbook(filename, )
    index_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'right'})
    header_format = workbook.add_format({'bold': True, 'num_format': '@', 'align': 'center'})
    num_format = workbook.add_format({'num_format': num_format})

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
                try:
                    sheet.write(row_i + 1, col_i + 1, val, num_format)
                except TypeError:
                    if pd.isna(val):
                        pass
                    else:
                        raise

        # freeze index row/column
        try:
            freeze = freeze_indexes[sheet_name]
        except:
            freeze = freeze_indexes
        if freeze:
            sheet.freeze_panes(0, 0)

    if close_workbook:
        workbook.close()
    return workbook


def convert_stats_csv_to_excel(
        intable: Union[str, pd.DataFrame],
        analysis_type: str,
        outpath: str = None,
        selected_comparisons: List[str] = False,

):
    """Read or recieve a table output by a tablulate_[method] function,
    where .columns.levels[0] is comparisons, and levels[1] are
    statistics. Output an Excel workbook with one tab per comparison,
    and renamed (to be clearer) statistic names.

    Args:
        intable: a DF or path to a DF
        analysis_type: (long) name of the analysis being loaded
        outpath: where the xlsx will be written, if None return {comp:DF}
            of formated tables.
        selected_comparisons: Write only specific comparisons from intable
    """

    addarrow = lambda s: s.replace('-', f' {ARROW} ')

    scorek = analysis_score_names[analysis_type]

    if type(intable) is not pd.DataFrame:
        fulltable = pd.read_csv(intable, header=[0, 1], index_col=0)
    else:
        fulltable = intable

    if selected_comparisons is False:
        selected_comparisons = fulltable.columns.levels[0]

    outtables = {}
    scores = fulltable.xs(scorek, axis=1, level=1)
    scores.columns = scores.loc[:, selected_comparisons].columns.map(addarrow)
    outtables[f"All {statistic_labels[scorek]}"] = scores
    for c in selected_comparisons:
        table = fulltable[c]
        table = table.loc[:, [scorek] + ['fdr', 'neg_p', 'pos_p']].sort_values('neg_p')
        table.columns = table.columns.map(statistic_labels)

        outtables[addarrow(c)] = table

    if outpath is not None:
        write_stats_workbook(
            outtables,
            outpath,
            num_format='0.00',
        )
    else:
        return outtables

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




def load_counts(file_path, index_col='guide', gene_col='gene') -> (pd.DataFrame, pd.Series):
    """Load count table, putting gene column into separate Series and dropping it
    from counts. Returns (cnt, guide_gene)
    """
    cnt = pd.read_csv(file_path, sep='\t', index_col=index_col)
    guide_gene = cnt[gene_col]
    cnt.drop(gene_col, axis=1, inplace=True)
    return cnt, guide_gene

def write_counts(file_path, counts:pd.DataFrame,
                 guide_gene:pd.Series, gene_col_name='gene'):
    cnt = counts.copy()
    cnt.insert(0, gene_col_name, guide_gene)
    cnt.to_csv(file_path, sep='\t')


def load_replicate_sample_details(xlsx) -> (pd.DataFrame, pd.DataFrame):
    """From a screen details Excel file (sheet='Sample details), return DFs
    containg replicate and sample details"""
    rep_deets = pd.read_excel(xlsx, sheet_name='Sample details', )
    rep_deets = rep_deets.set_index('Sequence name', drop=False)
    samp_deets = rep_deets.drop_duplicates('Sample').set_index('Sample', drop=False).drop(['Replicate', 'Sequence name'], 1)
    return rep_deets, samp_deets

def is_ntctrl(gene):
    g = gene

    return ((g is np.nan) or (not g) or g.startswith('ctrl_') or
            any([g.lower().startswith(x.lower()) for x in ('NonTargeting', 'Luciferase', 'Scrambled', 'control')]) or
            (g.startswith('nt') and g[2:].isnumeric()) or
            g.startswith('ntctrl') or g.startswith('NEG_CONTROL') or
            g.startswith('NonTargetingControlGuideForHuman') or
            (g.lower() in ['luciferase', 'ntctrl', 'non', 'negative_control', 'uncertain', 'lacz', 'egfp']) or
            'scambled' in g.lower() or (g.startswith('NonTargeting')))
