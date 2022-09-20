#from crispr_pipeline import *
from shutil import rmtree
import os



def test_pipeline(**kwargs):
    from crispr_pipeline import process_arguments, run_analysis
    import yaml, pathlib
    print(kwargs)
    p = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    print(p)
    testp = p/'tests'
    # outp = str(Path(p.parent/'tests'/'ran_test'))
    unprocessed_args = yaml.safe_load(open(testp/'test_pipeline.yaml'))
    from crispr_pipeline import process_arguments, run_analysis
    args = process_arguments(unprocessed_args)
    for k, v in kwargs.items():
        args[k] = v
    #print(args)
    run_analysis(**args)


def test_count_reads_read_fasta(fn):
    import count_reads

    with open(fn) as f:
        seqs = [s for s in count_reads.read_fasta(f)]
    print(len(seqs))

    print(seqs[:5])

def test_count_reads_count_reads(fn):
    from count_reads import count_reads_from_file

    res = count_reads_from_file(fn, (24, 44), )
    print(res.most_common(10))

def test_count_reads_count_batch(fn):
    from count_reads import count_batch
    count_batch(fn, (24,44), 'test_batch', )



# if __name__ == '__main__':
#     test_pipeline(charts_only=False)
#     print('That worked, hopefully')
fn = '/Users/johnc.thomas/mnt/thecluster/jct61/tst.fna'
test_count_reads_count_batch(fn)
