#!/usr/bin/env python
import os, sys, subprocess
import typing

import pandas as pd
from collections import Counter
import gzip
import argparse
from pathlib import Path, PosixPath, WindowsPath
from typing import Union, Tuple, List, Dict
import datetime
from typing import List
import multiprocessing
from functools import partial
# make it easy to find print statements used for debugging
debugprint = print

from crispr_tools.tools import maybe_its_gz

import logging
LOG = logging.getLogger(__name__, )
LOG.setLevel(logging.INFO) # need to set the logger to the lowest level...
log_formatter = logging.Formatter('%(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(log_formatter)
LOG.addHandler(stream_handler)

"""Functions for counting and mapping reads from multiple FASTQ files
to some kind of sequence library, e.g. a CRISPR guide library. 
Merging counts assumes that FQ filenames are in the format 
{sample_name}_L00?_R1_001.fastq[.gz]
and that the sequence to be mapped is in the same position in every read
(including the entire read being the sequence).
Produces dereplicated sequence counts (one file per sample) and then a single
file containing reads mapped to guide/gene from a library file."""

__version__ = '1.9'
# 1.9   fixed edge case when the guide id column in the library contains duplicated values
# 1.8.5.1 better way to not error if directory exists
# 1.8.5 don't error if directory exists
# 1.8.4 create directory from prefix if not exist
# 1.8.3 enabled custom library headers, pass using -j -n -g flags, defaults to
#       old behaviour (-j guide -n gene -g seq) if not provided
# 1.8.2 added checkpointing. If rawcounts file(s) already exist, then they won't
#       be re-generated. Does not check actual contents of existing files.
# 1.8.1 bugfix in merge samples. Won't work with runs that have more than 9 lanes
# 1.8.0 added back fuzzy matching, more efficiently this time
# 1.7.0 adding permenant logging
# 1.6.0 removed fuzzy counting
# 1.5.1 removed some unneccesary print calls
# 1.5.0 added support for FastA
# 1.4.1 efficiency when merging the files in map_reads
# 1.4.0 added fuzzy matching
# 1.3.4 bug in mapping function
# 1.3.3 fixed mapping so it properly truncates the filename to sample name
# 1.3.2 Bugfixes, documentation
# 1.3.1 returned column order of map_counts now 'guide', 'gene', [everything else]
#       fixed splitter issue
# v1.3.0 added map_counts
# v1.2.4 added printed info explaining sample merges.
# v1.2.3 changed func name to count_reads, added import to __init__, renamed to crispr_tools
# v1.2.2 bugfixes, started using pathlib, reimplimented ability to pass dir
# v1.2.1
#   added option to pass substring to slicer which will be used to start cutting
#   no longer checks line number/4 == sequences counted
# v1.2 filename prefix, removed directory walking, added merge_samples, removed merge_lanes
# v1.1 added argument parsing, options to merge lanes

#todo just ouput a single file
#todo use logging instead of stealing print
#todo more information about matches per file, per sample mapping stats
#todo unmerged counts with identical filenames overwrite each other currently in the
#todo check library seq len and file seq len match.


FILE_FMT = {'a':['.fna', '.fasta', '.fna.gz', '.fasta.gz'],
              'q': ('.fastq', '.fastq.gz', '.fq', '.fq.gz')}

ALL_FMT = []
for fmt in FILE_FMT.values():
    ALL_FMT.extend(fmt)


def read_fastq(file_obj):
    """assumes the 2nd line, and every 4th thereafter, is sequence"""
    pos = 3
    for line in file_obj:
        if pos == 4:
            pos = 0
            yield  line.strip()
        pos += 1


def read_fasta(file_obj, description_character='>'):
    """Read FastA file, sequence defined as lines following record description lines
    that begin with '>' by default. """
    while True:
        # not sure why I have to explicitly deal with StopIter here and not in read_fastq
        #   maybe .__next__() works better?
        try:
            line = next(file_obj)

        except StopIteration:
            return
        seq = []
        while line[0] != description_character:
            seq.append(line.strip())
            try:
                line = next(file_obj)
            except StopIteration as e:
                break
        s = ''.join(seq)
        if s:
            yield s


def count_reads_from_file(fn, slicer=(None, None), s_len=None, s_offset=0, file_format='infer') -> Counter:
    """Count single fastA/Q file's reads. Use count_batch for multiple files.
    Slicer should be a tuple of indicies.
    Returns a Counter.
    """
    if type(slicer[0]) == int:
        chop = slice(*slicer)
    seqs_count = Counter()
    failed_count = 0
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rt')
    else:
        f = open(fn)

    if file_format == 'infer':
        for fm, suffixes in FILE_FMT.items():

            if any([fn.endswith(suf) for suf in suffixes]):

                file_format = fm
        # else:
        #     raise RuntimeError('Unknown file format for '+fn, 'Please specify in args.')

    reader = {'q':read_fastq(f), 'a':read_fasta(f)}[file_format]

    # parse out the sequences from the FastQ
    for s in reader:
        # get the relevant part of the sequence
        if s_len is None:
            s = s[chop]
        else:
            try:
                i = s.index(slicer)+s_offset
                s = s[i:i+s_len]
            except:
                failed_count+=1
                continue

        seqs_count[s] += 1
    f.close()
    LOG.info(f"{fn} {len(seqs_count)} unique sequences, {sum(seqs_count.values())} reads")
    if s_len is not None:
        LOG.info(f"{failed_count} sequences did not contain subsequence")
    return seqs_count



def get_file_list(files_dir) -> List[os.PathLike]:
    """Pass single string or list of stings, strings that are files go on the
    file list, directories will have all files within (not recursive) will be
    added to the file list.

    A list of Path obj are returned."""
    LOG.debug(f'get_file_list in: {files_dir}')
    file_list = []
    # We want a list, not single string
    if type(files_dir) in (str, PosixPath, WindowsPath):
        files_dir = [files_dir]
    # convert to Path
    files = [Path(f) for f in files_dir]

    # get single list of Path, containing all listed files and files in listed dir
    # dir within `files` will be ignored.
    for fndir in files:
        if os.path.isdir(fndir):
            for fn in os.listdir(fndir):
                fn = fndir/fn
                if os.path.isfile(fn):
                    file_list.append(fn)
        else:
            file_list.append(fndir)
    file_list = [fn for fn in file_list if fn.name[0] != '.']
    LOG.debug(f"get_file_list out: {file_list}")
    return file_list

def file_exists_or_zero(fn):
    """Determine whether the output rawcounts file already exists and contains
    data. Returns True if the file exists and is not zero bytes. Returns False
    if the file does not exist or it does exist but is zero bytes. Used for
    checkpointing purposes. Count steps will skip if the output rawcounts file
    exists and is not size zero bytes (function returns True)."""
    if os.path.exists(fn) and os.path.getsize(fn) != 0:
        return True
    else:
        return False


def count_batch(fn_or_dir, slicer, fn_prefix='', seq_len=None, seq_offset=0, fn_suffix='.rawcount',
                fn_split='_R1_', merge_samples=False, just_go=False, quiet=False,
                file_type = 'infer', no_log_file=True, debug=False):
    """Write a table giving the frequency of all unique sequences from a fastq
    files (or fastq.gz).

    Output filenames are derived from input filenames. Outfn are returned.

    Merge_samples looks for file names containing "_L001_" and assumes all
    files with the same prefix are the same sample.

    Arguments:
        fn_or_dir:
            A list of files or dir. All files in given dir that end
            with .fastq or .fastq.gz or .fq will be counted.

        slicer (M,N)|subseq:
            Slice indicies to truncate sequences (zero indexed,
            end exclusive). Comma-sep ints.
            OR a subsequence from which slicing will happen.

        seq_len (int):
            Required when slicer is a sequence, determines output sequence
            lengths.

        fn_prefix:
            Prefix added to output files, can include absolute or
            relative paths.

        fn_suffix:
            Suffix added to output files, .txt added after. Default `rawcount`

        fn_split:
            String used to split filenames and form output file prefix.
            Default `_R1_`. Doesn't do anything if --merge-samples is used.

        merge_samples:
            Merge counts from files with identical sample names. Be
            careful not to double count decompressed & compressed files. Bool.

        allowed_extensions  ('fastq', 'fastq.gz', '.fq')

        no_log_file:
            Set to false to prevent a log file from being written


    """
    kwargs = locals()
    del kwargs['fn_or_dir']
    LOG.info(f"Kwargs:\n\t{kwargs}")
    if quiet:
        stream_handler.setLevel(logging.WARNING)
    elif debug:
        stream_handler.setLevel(logging.DEBUG)
        LOG.setLevel(logging.DEBUG)

    if not no_log_file:
        t = datetime.datetime.now()
        ts = t.strftime('%y%m%dh%Hm%Ms%S_%f')[:-2]
        file_handler = logging.FileHandler(f"{fn_prefix}.logfile.{ts}.txt")
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(log_formatter)
        LOG.addHandler(file_handler)

    t = os.path.split(fn_prefix)[0]
    t = os.path.expanduser(t)
    if t:
        if not os.path.isdir(t):
          LOG.info(f"Making directory {t}")
          os.makedirs(t)

    # accepts list of mixed files and dir, ends up as list of files
    file_list = get_file_list(fn_or_dir)

    # filter the file list
    # strings are easier to work with at this point
    file_list = [str(f) for f in file_list]

    #banned_files = [any([f.endswith(suf) for suf in allowed_extensions  if f not in file_list]]
    banned_files = [f for f in file_list if not any([f.endswith(suf) for suf in ALL_FMT])]

    if file_type == 'infer' and banned_files:
        LOG.warning(
            ('WARNING: The following files might not be fastq/a and will be ignored:\n'+
            ', '.join(banned_files))
        )
        file_list =  [f for f in file_list if f not in banned_files]

    # map filenames to samples
    if merge_samples:
        samples = set([f.split('_L001_')[0].split('/')[-1] for f in file_list if '_L001_' in f])

        file_dict = {s:[f for f in file_list if s+'_L00' in f] for s in samples}
        message_samples = ['Samples found:']
        for k,v in file_dict.items():
            message_samples.append(k)
            for f in v:
                message_samples.append('\t'+f)
        LOG.info('\n'.join(message_samples))
    else:
        LOG.info('input files:\n'+'\n'.join(file_list))

    if type(slicer[0]) is int:
        lengthstring = 'Length={}'.format(slicer[1]-slicer[0])
    else:
        lengthstring = 'Length={}'.format(seq_len)
    if not just_go:
        input(lengthstring+'\nPress enter to go...')
    else:
        LOG.info(lengthstring)

    # called in the main loop
    
    def get_outfn(fn_base):
        if fn_prefix:
            if fn_prefix[-1] == '/':
                fs = '{}{}{}.txt'
            else:
                fs = '{}.{}{}.txt'
            outfn = fs.format(fn_prefix, fn_base, fn_suffix)
        else:
            outfn = '{}.{}.txt'.format(fn_base, fn_suffix)
        return(outfn)
    
    def write_count(a_cnt, fn_base):
        outfn = get_outfn(fn_base)
        with open(outfn, 'w') as f:
            for s, n in a_cnt.most_common():
                f.write('{}\t{}\n'.format(s,n))
        return outfn
    # used if library mapping being done
    out_files = []

    # main loop(s)
    if merge_samples:
        LOG.debug(str(file_dict))
        for samp, fn_or_dir in file_dict.items():
            outfn = get_outfn(samp)
            out_files.append(outfn)
            if not file_exists_or_zero(outfn):
                cnt = Counter()
                for fn in fn_or_dir:
                    cnt += count_reads_from_file(fn, slicer, seq_len, seq_offset)
                write_count(cnt, samp)
            else:
                LOG.info('Output counts file already exists. Skipping: ' + outfn)
    else:
        LOG.debug(str(file_list))
        for fn in file_list:
            samp = fn.split(fn_split)[0].split('/')[-1]
            outfn = get_outfn(samp)
            out_files.append(outfn)
            if not file_exists_or_zero(outfn):
                cnt = count_reads_from_file(fn, slicer, seq_len, seq_offset, file_type)
                out_files.append(
                    write_count(cnt, samp)
                )
            else:
                LOG.info('Output counts file already exists. Skipping: ' + outfn)

    return out_files

def get_count_table_from_file_list(file_list:List[Path], splitter='.raw',
                                   remove_prefix=True) -> pd.DataFrame:
    """Put the contents of a list of raw count filepaths into a DF.
    Args:
        file_list: list of files
        splitter: The substring that follows the sample name, used to extract
        the sample name.
        remove_prefix: If True, <sample_name>.split('.')[1] is performed.
    Returns:
        pd.DataFrame indexed by sequence with integer counts from each file.
        Column names derived from the file names.
    """
    rawcnt = {} #  will be cast to DF
    for fn in file_list:
        fn = Path(fn)
        # filtering of fn done before here
        sn = fn.name.split(splitter)[0]
        if remove_prefix:
            sn = sn.split('.')[1]
        #LOG.info(f'sample header: {sn}')
        rawcnt[sn] = pd.read_csv(fn, index_col=0, header=None, sep='\t')[1]
        # rawcnt = pd.concat([rawcnt, samp], 1)
    LOG.info(f'sample headers: {list(rawcnt.keys())}')
    return pd.DataFrame(rawcnt).fillna(0).astype(int)


def map_allowing_mismatch(query_sequences:typing.Collection[str],
                          lib_sequences:typing.Collection[str],
                          processors=4) -> Dict[str, str]:
    """Return dictionary of query sequences that uniquely match, with 1 or 0 mismatches,
    a library sequence. Queries that don't match, or are 1 mismatch from two library
    sequences are not returned."""
    # use a Go program to map close seqs quickly
    # write temp files
    with open('tmp_unmapped.txt', 'w') as f:
        f.write('\n'.join(query_sequences))

    with open('tmp_libseqs.txt', 'w') as f:
        f.write('\n'.join(lib_sequences))

    from pkg_resources import resource_filename
    fuzzy_prog = resource_filename(__name__, "fuzzyTwoLists")

    output = subprocess.run(
        [fuzzy_prog,
         'tmp_unmapped.txt', 'tmp_libseqs.txt', 'tmp.out', '1', f'{processors}'],
        capture_output=True
    )

    LOG.info(output.stdout.decode())

    # output is in a tab sep list of query->library seq
    # only containing positive hits
    now_mapped = {}
    with open('tmp.out') as f:
        for line in f:
            line = line.strip()
            a, b = line.split('\t')
            now_mapped[a] = b
    return now_mapped


def map_counts(fn_or_dir:Union[str, List[str]], lib:Union[str, pd.DataFrame],
               seqhdr='seq', guidehdr='guide', genehdr='gene',
               report=False, splitter='.raw',
               remove_prefix=True, out_fn=None,
               allow_mismatch=False, processors=4) -> pd.DataFrame:
    """
    Map guide sequences in a set of files containing guide sequences and abundance
    using a library file.

    Returns a DataFrame indexed by guide name, and columns giving guide abundance
    per replicate.

    Args:
        fn_or_dir: Directory, or list of directories/files from which rawcount
            files will be processed.
        lib: library file or DF containing the guide sequences and guide.
            If DF all 3 columns must be present (index is ignored).
        seqhdr: library column header indicating guide sequences
        guidehdr: library column header indicating guide names
        genehdr: library column header indicating gene names
        report: When True, print stats about the mapping
        splitter: column headers in output table will be
            <prefix>.<filename>.split(splitter)[0]
        remove_prefix: If True, <sample_name>.split('.')[1]
        out_fn: Write DF to filename, if provided
        allow_mismatch: Allow one mismatch when mapping
        processors: number of proc used when allow_mismatch is True

    Returns:
        pd.DataFrame
    """

    if type(lib) in (str, PosixPath, WindowsPath):
        lib = str(lib)
        if lib.endswith('.csv') or lib.endswith('.csv.gz'):
            sep = ','
        else:
            sep = '\t'
        lib = maybe_its_gz(lib)
        lib = pd.read_csv(lib, sep=sep)

    for hdr in (seqhdr, guidehdr, genehdr):
        if hdr not in lib.columns:
            raise RuntimeError(
                f'Library header doesnt match expeted: "{hdr}" not found.\n',
                f'Library header: {lib.columns}\n',
                f'Expected headers: seqhdr={seqhdr}, guidehdr={guidehdr}, genehdr={genehdr}'
            )

    lib.set_index(seqhdr, drop=False, inplace=True)

    # else the library should already be in a useable form.

    # get a single table
    file_list = get_file_list(fn_or_dir)
    rawcnt = get_count_table_from_file_list(file_list, splitter, remove_prefix)

    # map the matching counts
    cnt = rawcnt.reindex(lib[seqhdr], fill_value=0)

    if allow_mismatch:
        # get unmapped sequences
        m =rawcnt.index.isin(lib[seqhdr])
        unmapped_seqs = rawcnt.loc[~m].index
        now_mapped = map_allowing_mismatch(unmapped_seqs, lib[seqhdr], processors)

        # sum up the abundance of the newly mapped reads
        raw_now_mapped = rawcnt.loc[now_mapped.keys()]
        raw_now_mapped.index = raw_now_mapped.index.map(now_mapped)
        gained_counts = raw_now_mapped.groupby(raw_now_mapped.index).sum()

        # report
        frac_gained = gained_counts.sum().sum() / cnt.sum().sum()
        LOG.info(f"Percentage abundance gained by fuzzy matching: {frac_gained:.2%}")

        # add the gained counts, reindexing to deal with duplicate seqs in cnt
        cnt += gained_counts.reindex(cnt.index, fill_value=0)

    # Use the guide names as the index as the *should* always be unique,
    #   even if the sgRNA are not
    cnt.index = lib[guidehdr]

    if report:
        prop = cnt.sum().sum()/rawcnt.sum().sum()
        per_samp_prop = (cnt.sum()/rawcnt.sum()*100).apply(lambda n: f"{n:.3}%")
        LOG.info("{:.3}% of reads map.".format(prop*100))
        no_hits = (cnt == 0).all(1).sum()
        LOG.info("{:.3}% ({}) of library guides not found.".format(
            no_hits/cnt.shape[0]*100, no_hits
        ))
        LOG.info(
            f"Per replicate mapping (%):\n{per_samp_prop}"
        )

    # Add gene column, guide
    cnt.insert(0, 'gene', lib.set_index(guidehdr)[genehdr])
    cnt = cnt.drop_duplicates()

    if out_fn:
        cnt.to_csv(out_fn, sep='\t')
    return cnt


if __name__ == '__main__':
    print('count_reads.py version', __version__)

    #print('v', __version__)
    parser = argparse.ArgumentParser(description='Count unique sequences in FASTQs. Assumes filenames are {sample_name}_L00?_R1_001.fastq[.gz]')
    parser.add_argument('files', nargs='+',
                        help="A list of files or dir. All files in given dir that end with valid suffixes (inc. gzips) will be counted.")

    parser.add_argument('-s', '--slice', metavar='M,N',required=True,
                        help='Slice indicies to truncate sequences (zero indexed, not end-inclusive). Comma-sep numbers. Required.',)
    parser.add_argument('-p', '--prefix', metavar='FN_PREFIX', required=True,
                        help="Prefix added to output files, can include absolute or relative paths.")
    parser.add_argument('--suffix', default='.rawcount', metavar='FN_SUFFIX',
                        help="Suffix added to output files, .txt will always be added after. Default `.rawcount`")
    parser.add_argument('--fn-split', default='_R1_', metavar='STR',
                        help="String used to split filenames and form output file prefix. Default `_R1_`." \
                             "Doesn't do anything if --merge-samples is used.")
    parser.add_argument('--merge-samples', action='store_true', default=False,
                        help="Merge counts from files with identical sample names.")
    parser.add_argument('--just-go', action='store_true', default=False, help="Don't wait for confirmation.")
    parser.add_argument('--quiet', action='store_true', default=False, help="Don't print helpful messages, enables --just-go.")
    parser.add_argument('--library', default=None, metavar='LIB_PATH',
                        help="Pass a library file and a single mapped reads count file will be output with" \
                             "the name [FN_PREFIX].counts.tsv")
    parser.add_argument('-g', '--seqhdr', default='seq',
                        help="Name  of column in the library file that contains sequences.")
    parser.add_argument('-j', '--guidehdr', default='guide',
                        help="Name  of column in the library file that contains sequences.")
    parser.add_argument('-n', '--genehdr', default='gene',
                        help="Name  of column in the library file that contains sequences.")
    parser.add_argument('--allow-mismatch', action='store_true', default=False,
                        help='Allow one mismatch when mapping reads to library.')
    parser.add_argument('--file_type', metavar='FILE_TYPE', choices={'a', 'q', 'infer'}, default='infer',
                        help='If your files end in something weird, use "a" for fastA or "q" for fastQ. '
                             'Cant handle mixed file types.'
                            'Can handle .gz. By default its infered by the file names.')
    parser.add_argument('--no-log-file', action='store_true', help='Prevent a log file from being written')
    parser.add_argument('--debug', action='store_true',
                        help='Print/log additional debug messages.')

    clargs = parser.parse_args()

    if clargs.library:
        if not os.path.isfile(clargs.library):
            raise RuntimeError(f"Library file not found: {clargs.library}")

        lib = pd.read_table(clargs.library, index_col=clargs.guidehdr, sep=None, engine="python")
        if lib.index.duplicated().any():
            print('WARNING: Duplicated guide names found in library, first one will be kept, all others discarded.')

    assert all([os.path.isfile(f) for f in clargs.files])

    # slices list of input files, or dir
    slicer = [int(n) for n in clargs.slice.split(',')]
    
    os.makedirs(os.path.dirname(clargs.prefix), exist_ok=True)

    written_fn = count_batch(fn_or_dir=clargs.files,
                             slicer=slicer,
                             fn_prefix=clargs.prefix,
                             fn_suffix=clargs.suffix,
                             fn_split=clargs.fn_split,
                             merge_samples=clargs.merge_samples,
                             just_go=clargs.just_go,
                             quiet=clargs.quiet,
                             no_log_file=clargs.no_log_file,
                             file_type=clargs.file_type)

    if clargs.library:
        map_counts(written_fn, clargs.library, report=True, remove_prefix=True,
                out_fn=clargs.prefix+'.counts.tsv', splitter=clargs.suffix,
                allow_mismatch=clargs.allow_mismatch, seqhdr=clargs.seqhdr,
                guidehdr=clargs.guidehdr, genehdr=clargs.genehdr)

