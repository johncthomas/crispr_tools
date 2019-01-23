import os, sys
from collections import Counter
import gzip
import argparse
from pathlib import Path

__version__ = '1.2.4'
# v1.2.4 added printed info explaining sample merges.
# v1.2.3 changed func name to count_reads, added import to __init__, renamed to crispr_tools
# v1.2.2 bugfixes, started using pathlib, reimplimented ability to pass dir
# v1.2.1
#   added option to pass substring to slicer which will be used to start cutting
#   no longer checks line number/4 == sequences counted
# v1.2 filename prefix, removed directory walking, added merge_samples, removed merge_lanes
# v1.1 added argument parsing, options to merge lanes

#todo just ouput a single file what
#todo use **kwargs to pass parsed args to count_batch


def count_reads(fn, slicer=(None, None), s_len=None, s_offset=0, ):
    """Count single .fq reads. Use count_batch for multiple files.
    Slicer should be a tuple of indicies.
    Returns a Counter.

    s_len and s_offset are depreciated, use Cutadapt instead.
    """
    if type(slicer[0]) == int:
        chop = slice(*slicer)
    seqs = Counter()
    lns = 0
    failed_count = 0
    if fn.endswith('.gz'):
        f = gzip.open(fn, 'rt')
    else:
        f = open(fn)
    for line in f:
        lns += 1
        if line[0] == '@':
            s = f.__next__()
            s = s.replace('\n', '')
            lns+=1
            if s_len is None:
                s = s[chop]
            else:
                try:
                    i = s.index(slicer)+s_offset
                    s = s[i:i+s_len]
                except:
                    failed_count+=1
                    continue
            seqs[s] += 1
    f.close()
    print(fn, len(seqs), 'sequences counted.')
    if s_len is not None:
        print(failed_count, 'sequences did not contain subsequence')
    return seqs

def count_batch(files, slicer, fn_prefix='', seq_len=None, seq_offset=0, fn_suffix='rawcount',
                fn_split='_R1_', merge_samples=False, just_go=False, quiet=False,
                allowed_extensions = ('.fastq', 'fastq.gz', '.fq')):
    """Write a table giving the frequency of all unique sequences from a fastq
    files (or fastq.gz).

    Output filenames are derived from input filenames.

    Merge_samples looks for file names containing "_L001_" and assumes all
    files with the same prefix are the same sample.

    Arguments:
        files:
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

    if quiet:
        print = lambda x: None
        just_go = True
    else:
        print = print

    # accepts list of mixed files and dir, ends up as list of files
    file_list = []
    if type(files) is str:
        files = [files]
    files = [Path(f) for f in files]
    for fndir in files:
        if fndir.is_dir():
            for fn in os.listdir(fndir):
                # some syntactic sugar from pathlib
                fn = fndir / fn
                if os.path.isfile(fn):
                    file_list.append(fn)
        else:
            file_list.append(fndir)

    # filter the file list
    # strings are easier to work with at this point
    file_list = [str(f) for f in file_list]
    file_list = [ f for f in file_list if any([f.endswith(suf) for suf in allowed_extensions])]

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
                fs = '{}{}.{}.txt'
            else:
                fs = '{}.{}.{}.txt'
            outfn = fs.format(fn_prefix, fn_base, fn_suffix)
        else:
            outfn = '{}.{}.txt'.format(fn_base, fn_suffix)
        with open(outfn, 'w') as f:
            for s, n in a_cnt.most_common():
                f.write('{}\t{}\n'.format(s,n))

    # main loop(s)
    if merge_samples:
        for samp, files in file_dict.items():
            cnt = Counter()
            for fn in files:
                cnt += count_reads(fn, slicer, seq_len, seq_offset)
            write_count(cnt, samp)
    else:
        for fn in file_list:
            cnt = count_reads(fn, slicer, seq_len, seq_offset)
            write_count(cnt, fn.split(fn_split)[0].split('/')[-1])


        # # merge lanes if required
        # # iteratively look for files containing _L002_, _L003_, etc
        # lane_i = 2
        # while merge_lanes:
        #     next_fn = fn.replace('_L001_', '_L{0:03d}_'.format(lane_i))
        #     if os.path.isfile(next_fn):
        #         cnt += count_seqs(next_fn, slicer)
        #     else: # if the next file doesn't exist, we are done.
        #         break
        #     lane_i += 1

        # build output filename and write



if __name__ == '__main__':
    print('v', __version__)
    parser = argparse.ArgumentParser(description='Count unique sequences in FASTQs. Assumes filenames are {sample_name}_L00?_R1_001.fastq[.gz]')
    parser.add_argument('files', nargs='+',
                        help="A list of files or dir. All files in given dir that end with .fastq or .fastq.gz or .fq will be counted.")
    parser.add_argument('-s', metavar='M,N',
                        help='Slice indicies to truncate sequences (zero indexed, not end-inclusive). Comma-sep numbers. Required.',
                        required=True)
    parser.add_argument('-f', default='rawcount', metavar='FN_SUFFIX',
                        help="Suffix added to output files, .txt will always be added after. Default `rawcount`")
    parser.add_argument('-p', default='', metavar='FN_PREFIX',
                        help="Prefix added to output files, can include absolute or relative paths.")
    parser.add_argument('--fn-split', default='_R1_', metavar='STR',
                        help="String used to split filenames and form output file prefix. Default `_R1_`. Doesn't do anything if --merge-samples is used.")
    parser.add_argument('--merge-samples', action='store_true', default=False,
                        help="Merge counts from files with identical sample names. Be careful not to double count decompressed & compressed files.")
    parser.add_argument('--just-go', action='store_true', default=False, help="Don't wait for confirmation.")
    parser.add_argument('--quiet', action='store_true', default=False, help="Don't print helpful messages, enables --just-go.")
    clargs = parser.parse_args()

    # slices list of input files, or dir
    slicer = [int(n) for n in clargs.s.split(',')]

    count_batch(clargs.files, slicer, clargs.p, None, 0, clargs.f, clargs.fn_split, clargs.merge_samples, clargs.just_go, clargs.quiet)

