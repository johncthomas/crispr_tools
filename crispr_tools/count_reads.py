#!/usr/bin/env python
import os, sys
import pandas as pd
from collections import Counter
import gzip
import argparse
from pathlib import Path, PosixPath, WindowsPath
from typing import Union, Tuple, List, Dict
try:
    from Levenshtein import distance as levendist
except ImportError:
    print('Levenshtein package not present, fuzzy matching not available')
    def levendist():
        raise ModuleNotFoundError('Levenshtein package not installed')
from typing import List
import multiprocessing
from functools import partial
# make it easy to find print statements used for debugging
debugprint = print

"""Functions for counting and mapping reads from multiple FASTQ files
to some kind of sequence library, e.g. a CRISPR guide library. 
Assumes that FQ filenames are in the format {sample_name}_L???_R1_001.fastq[.gz]
and that the sequence to be mapped is in the same position in every read.
It's probably not very useful if either of these things are false.

Produces dereplicated sequence counts (one file per sample) and then a single
file containing reads mapped to guide/gene from a library file."""

__version__ = '1.5.1'
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
#todo order reps by

FILE_FMT = {'a':['.fna', '.fasta', '.fna.gz', '.fasta.gz'],
              'q': ('.fastq', '.fastq.gz', '.fq', '.fq.gz')}

ALL_FMT = []
for fmt in FILE_FMT.values():
    ALL_FMT.extend(fmt)


def read_fastq(file_obj):
    for line in file_obj:
        if line[0] == '@':
            s = file_obj.__next__()
            s = s.replace('\n', '')
            yield s

def read_fasta(file_obj):

    while True:
        # not sure why I have to manually deal with StopIter here and not in read_fastq
        #   maybe .__next__() works better?
        try:
            line = next(file_obj)

        except StopIteration:
            return
        seq = []
        while line[0] != '>':
            seq.append(line.strip())
            try:
                line = next(file_obj)
            except StopIteration as e:
                break
        s = ''.join(seq)
        if s:
            yield s


def count_reads(fn, slicer=(None, None), s_len=None, s_offset=0, file_format='infer') -> Counter:
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
    print(fn, len(seqs_count), 'unique sequences, ', sum(seqs_count.values()), 'reads')
    if s_len is not None:
        print(failed_count, 'sequences did not contain subsequence')
    return seqs_count


def fuzzy_match(query_seq, lib_seqs, max_dist=3):
    """Check that exactly one library sequence is sufficiently close
    to the given sequence.
    If so, return the library sequence that is close,
    otherwise return None."""

    closeness = lib_seqs.map(lambda x: levendist(x, query_seq))
    # check the closest seq(s) are sufficiently close
    if min(closeness) > max_dist:
        return None
    # check only one is close enough
    close_mask = closeness == min(closeness)
    if sum(close_mask) == 1:
        return lib_seqs.loc[close_mask].values[0]
    else:
        # explicitly...
        return None


def get_file_list(files_dir) -> List[os.PathLike]:
    """Pass single string or list of stings, strings that are files go on the
    file list, directories will have all files within (not recursive) will be
    added to the file list.

    A list of Path obj are returned."""

    file_list = []
    # just in case trying to make a list of a single string...
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
    #print(file_list)
    return file_list


def count_batch(fn_or_dir, slicer, fn_prefix='', seq_len=None, seq_offset=0, fn_suffix='.rawcount',
                fn_split='_R1_', merge_samples=False, just_go=False, quiet=False,
                file_type = 'infer'):
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
    """
    global print
    if quiet:
        print = lambda x: None
        just_go = True
    else:
        print = print

    t = os.path.split(fn_prefix)[0]
    t = os.path.expanduser(t)
    if t:
        assert os.path.isdir(t)

    # accepts list of mixed files and dir, ends up as list of files
    file_list = get_file_list(fn_or_dir)

    # filter the file list
    # strings are easier to work with at this point
    file_list = [str(f) for f in file_list]
    #print(file_list)

    #banned_files = [any([f.endswith(suf) for suf in allowed_extensions  if f not in file_list]]
    banned_files = [f for f in file_list if not any([f.endswith(suf) for suf in ALL_FMT])]

    if file_type == 'infer' and banned_files:
        print(
            'WARNING, the following files might not be fastq/a and will cause the pipeline to crash:\n'
            ', '.join(banned_files)
        )
        file_list =  [f for f in file_list if f not in banned_files]

    # map filenames to samples
    if merge_samples:
        samples = set([f.split('_L001_')[0].split('/')[-1] for f in file_list if '_L001_' in f])

        file_dict = {s:[f for f in file_list if s in f] for s in samples}
        print('Samples found:',)
        for k,v in file_dict.items():
            print(k)
            for f in v:
                print('\t'+f)
    else:
        print('input files')
        print('\n'.join(file_list))

    if type(slicer[0]) is int:
        lengthstring = 'Length={}'.format(slicer[1]-slicer[0])
    else:
        lengthstring = 'Length={}'.format(seq_len)
    if not just_go:
        input(lengthstring+'\nPress enter to go...')
    else:
        print(lengthstring, '\n')

    # called in the main loop
    def write_count(a_cnt, fn_base):
        if fn_prefix:
            if fn_prefix[-1] == '/':
                fs = '{}{}{}.txt'
            else:
                fs = '{}.{}{}.txt'
            outfn = fs.format(fn_prefix, fn_base, fn_suffix)
        else:
            outfn = '{}.{}.txt'.format(fn_base, fn_suffix)
        with open(outfn, 'w') as f:
            for s, n in a_cnt.most_common():
                f.write('{}\t{}\n'.format(s,n))
        return outfn
    # used if library mapping being done
    out_files = []

    # main loop(s)
    if merge_samples:
        #debugprint(file_dict)
        for samp, fn_or_dir in file_dict.items():
            cnt = Counter()
            for fn in fn_or_dir:
                cnt += count_reads(fn, slicer, seq_len, seq_offset)
            out_files.append(
                write_count(cnt, samp)
            )
    else:
        #debugprint(file_list)
        for fn in file_list:
            cnt = count_reads(fn, slicer, seq_len, seq_offset, file_type)
            out_files.append(
                write_count(cnt, fn.split(fn_split)[0].split('/')[-1])
            )

    return out_files

def get_count_table_from_file_list(file_list:List[Path], splitter='.raw', remove_prefix=True) -> pd.DataFrame:
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
        print('sample header:', sn)
        rawcnt[sn] = pd.read_csv(fn, index_col=0, header=None, sep='\t')[1]
        # rawcnt = pd.concat([rawcnt, samp], 1)

    return pd.DataFrame(rawcnt).fillna(0).astype(int)


def map_counts(fn_or_dir:Union[str, List[str]], lib:Union[str, pd.DataFrame],
               seqhdr='seq', guidehdr='guide', genehdr='gene',
               drop_unmatched=False, report=False, splitter='.raw',
               remove_prefix=True, out_fn=None, fuzzy=False, fuzzy_threads=None,
               ignore_singletons=False) -> pd.DataFrame:
    """
    Map guide sequences in a set of files containing guide sequences and abundance
    using a library file.

    Returns a DataFrame indexed by guide name, and columns giving guide abundance
    per replicate.

    Args:
        fn_or_dir: Directory, or list of directories/files from which rawcount
            files will be processed.
        lib: library file or DF containing the guide sequences and guide
        seqhdr: library column header indicating guide sequences
        guidehdr: library column header indicating guide names
        genehdr: library column header indicating gene names
        drop_unmatched: When True reads matching no sequence in the library
            file are dropped.
        report: When True, print stats about the mapping
        splitter: column headers in output table will be
            <prefix>.<filename>.split(splitter)[0]
        remove_prefix: drop the prefix from column headers
        out_fn: Write DF to filename, if provided
        fuzzy: If True use fuzzy matching... Don't set to True
        fuzzy_threads: ...
        ignore_singletons: Don't attempt to map sequences that only appear
            once in all the rawcount files. Speeds things up, singletons
            represent ~40% of unique sequences but there's only one of them.

    Returns:
        pd.DataFrame
    """

    if fuzzy and fuzzy_threads is None:
        fuzzy_threads = multiprocessing.cpu_count() - 1


    if type(lib) in (str, PosixPath, WindowsPath):
        lib = str(lib)
        if lib.endswith('.csv'):
            sep = ','
        else:
            sep = '\t'
        lib = pd.read_csv(lib, sep)
        if seqhdr in lib.columns:
            lib.set_index(seqhdr, drop=False, inplace=True)
        else:
            lib.set_index(lib.columns[0])
    # else the library should be in a useable form.

    # write a single table
    file_list = get_file_list(fn_or_dir)
    file_list = [f for f in file_list if splitter in str(f)]

    rawcnt = get_count_table_from_file_list(file_list, splitter, remove_prefix)

    # we should have a table indexed by the sequences
    if fuzzy:
        if type(fuzzy) == int:
            maxdist = fuzzy
        else:
            maxdist = 3
        if ignore_singletons:
            rawcnt = rawcnt.loc[rawcnt.sum(1) > 1]

        seqs = pd.Series(rawcnt.index, index=rawcnt.index)
        chunk_sz = seqs.shape[0] // fuzzy_threads

        # asynchronously match
        fz_results = []
        with multiprocessing.Pool(fuzzy_threads) as pool:
            for chunkn in range(fuzzy_threads):
                subseqs = seqs.iloc[chunkn * chunk_sz:(1 + chunkn) * chunk_sz]
                fz_results.append(
                    pool.apply_async(
                        subseqs.apply,
                        args=(fuzzy_match,),
                        kwds=dict(
                            lib_seqs=pd.Series(lib.index),
                            max_dist=maxdist
                        )
                    )
                )

            for r in fz_results:
                matched = r.get().dropna()
                rawcnt.loc[matched.index, 'fuzzy'] = matched

            rawcnt = rawcnt.groupby('fuzzy').sum()
            rawcnt.index.name = 'seq'

    # the absent guides
    missing = lib.loc[~lib.index.isin(rawcnt.index), :].index

    # get the present guides
    matches = rawcnt.loc[rawcnt.index.isin(lib.index), :].index
    if report:
        prop = rawcnt.loc[matches, :].sum().sum()/rawcnt.sum().sum()
        print("{:.3}% of reads map.".format(prop*100))
        print("{:.3}% ({}) of library guides not found.".format(
            missing.shape[0] / lib.shape[0] *100, missing.shape[0]
        ))
    #cnt = rawcnt.loc[matches, :].copy()
    rawcnt.loc[matches, 'guide'] = lib.loc[matches, guidehdr]
    rawcnt.loc[matches, 'gene'] = lib.loc[matches, genehdr]

    missingdf = pd.DataFrame(index=missing, columns=rawcnt.columns)
    missingdf.loc[:, :] = 0
    missingdf.loc[missing, 'guide'] = lib.loc[missing, guidehdr]
    missingdf.loc[missing, 'gene'] = lib.loc[missing, genehdr]
    rawcnt = rawcnt.append(missingdf)

    if drop_unmatched:
        cnt = rawcnt.loc[lib.index, :].copy()
    else:
        cnt = rawcnt

    # sort out columns
    cols = list(cnt.columns)
    cnt = cnt.reindex([guidehdr, genehdr] + cols[:-2], axis='columns',)
    cnt.set_index(guidehdr, inplace=True)

    if out_fn:
        cnt.to_csv(out_fn, sep='\t')
    return cnt





# os.chdir('/Users/johnc.thomas/thecluster/jct61/counts/nomask')
# map_counts(
#     'tst', '/Users/johnc.thomas/thecluster/jct61/crispr_libraries/Kinase_gRNA_library_no_duplicates.csv',
#     drop_unmatched=True, report=True, out_fn='tst/tstout2.tsv'
#
# )




class MAIN:
    pass
if __name__ == '__main__':
    print('count_reads.py version', __version__)
    cpus = multiprocessing.cpu_count()
    if cpus > 10:
        cpus = 10

    #print('v', __version__)
    parser = argparse.ArgumentParser(description='Count unique sequences in FASTQs. Assumes filenames are {sample_name}_L00?_R1_001.fastq[.gz]')
    parser.add_argument('files', nargs='+',
                        help="A list of files or dir. All files in given dir that end with valid suffixes (inc. gzips) will be counted.")

    parser.add_argument('-s', metavar='M,N',
                        help='Slice indicies to truncate sequences (zero indexed, not end-inclusive). Comma-sep numbers. Required.',
                        required=True)
    parser.add_argument('--suffix', default='.rawcount', metavar='FN_SUFFIX',
                        help="Suffix added to output files, .txt will always be added after. Default `.rawcount`")
    parser.add_argument('-p',  metavar='FN_PREFIX',
                        help="Prefix added to output files, can include absolute or relative paths.")
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
    parser.add_argument('--fuzzy-matching', type=int, default=False, metavar='N',
                        help='Count sequences with up to N differences from the library (excludes ambiguous reads' \
                        'that fuzzily match multiple lib seqs)')
    parser.add_argument('--cpus', type=int, default=cpus, help="Number of processes when doing fuzzy matching. "\
                        "Default set to what's available or 10, whichever's lower.")
    parser.add_argument('--file_type', metavar='FILE_TYPE', choices={'a', 'q', 'infer'}, default='infer',
                        help='If your files end in something weird, use "a" for fastA or "q" for fastQ. Cant handle mixed file types.'
                        'Can handle .gz. By default its infered by the file names.')
    clargs = parser.parse_args()

    if clargs.library:
        assert os.path.isfile(clargs.library)

    assert all([os.path.isfile(f) for f in clargs.files])

    # slices list of input files, or dir
    slicer = [int(n) for n in clargs.s.split(',')]

    written_fn = count_batch(clargs.files, slicer, clargs.p, None, 0, clargs.suffix, clargs.fn_split, clargs.merge_samples,
                clargs.just_go, clargs.quiet,)

    if clargs.library:
        map_counts(written_fn, clargs.library, drop_unmatched=True, report=True, remove_prefix=True,
                out_fn=clargs.p+'.counts.tsv', splitter=clargs.suffix,
                   fuzzy=clargs.fuzzy_matching, fuzzy_threads=clargs.cpus)

