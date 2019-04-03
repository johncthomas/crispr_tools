from crispr_pipeline import *
from shutil import rmtree
import os



def test_pipeline(**kwargs):
    print(kwargs)
    p = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    print(p)
    testp = p/'tests'
    # outp = str(Path(p.parent/'tests'/'ran_test'))

    args = process_arguments(testp / 'test_pipeline.yaml')
    for k, v in kwargs.items():
        args[k] = v
    #print(args)
    run_analysis(**args)

if __name__ == '__main__':
    test_pipeline(charts_only=False)
    print('That worked, hopefully')