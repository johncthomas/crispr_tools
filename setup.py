from setuptools import setup, find_packages
# from os import path
# import os

#todo remove shapely requirement from mahal calc
#todo sort out drugz, make it part of this repository.

with open('./crispr_tools/version.py') as f:
    v = f.readline()
    # account for either string char
    if '"' in v:
        sep = '"'
    else:
        sep = "'"
    __version__ = v.split(sep)[-2]


setup(
    name = 'crispr_tools',
    version = __version__,
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    #package_data={'files':['files/Hart2017_TableS2_core_genes.txt']},
    install_requires = [
        'numpy', 'scipy', 'pandas',
        'statsmodels', 'attrdictionary', 'xlsxwriter', 'openpyxl',
        'scikit-learn', 'pyaml'
    ],
    extras_require={
        'plotting':['matplotlib', 'adjustText', 'seaborn'],
        'depreciated':['yaml']
    },
    
    python_requires = '>3.8',
    scripts=['crispr_tools/crispr_pipeline.py', 'crispr_tools/count_reads.py'],
    include_package_data=True,
)


