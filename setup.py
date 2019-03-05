from setuptools import setup, find_packages
from os import path

#todo use yaml instead of xlsx and remove xlrd requirement

setup(
    name = 'crispr_tools',
    version = '1.7.6b6',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 'pandas', 'jacks', 'xlrd', 'shapely'
    ],
    python_requires = '>=3.6'

    #scripts=['count_reads']
    #module=['count_reads']
)
