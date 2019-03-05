from pipeline import *
from shutil import rmtree
import os
def test_pipeline(**kwargs):
    p = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    #print(p)
    testp = str(p.parent/'tests')
    outp = str(Path(p.parent/'tests'/'ran_test'))
    if os.path.isdir(outp):
        rmtree(outp)
    run_analysis(testp+'/test_counts.tsv', testp+'/test_pipeline.xlsx', outp, 'itsatest',
             **kwargs)

if __name__ == '__main__':
    test_pipeline(charts_only=False)
    print('That worked, hopefully')