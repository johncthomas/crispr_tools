from setuptools import setup, find_packages
from os import path
setup(
    name = 'crispr_tools',
    version = '1.7.2b1',
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 'pandas',
    ],

    #scripts=['count_reads']
    #module=['count_reads']
)
