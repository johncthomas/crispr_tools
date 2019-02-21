from pipeline import *
def test_pipeline(**kwargs):
    p = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    testp = str(p.parent/'tests')
    run_analysis(testp+'/test_counts.tsv', testp+'/test_pipeline.xlsx', testp+'/ran_test', 'itsatest',
             labeldep=10, labelenr=10,  **kwargs)

if __name__ == '__main__':
    test_pipeline(charts_only=False)