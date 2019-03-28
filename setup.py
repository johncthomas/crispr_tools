from setuptools import setup, find_packages
# from os import path
# import os

#todo remove shapely requirement from mahal calc

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
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 'pandas', 'shapely', 'pyyaml', 'statsmodels'
    ],
    python_requires = '>=3.6',
    scripts=['crispr_tools/crispr_pipeline.py', 'crispr_tools/count_reads.py'],
    include_package_data=True,
)
