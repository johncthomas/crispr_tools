from setuptools import setup, find_packages
# from os import path
# import os

#todo use yaml instead of xlsx and remove xlrd requirement
#todo remove shapely requirement from mahal calc

with open('./crispr_tools/version.txt') as f:
    __version__ = f.readline().replace('\n', '')

setup(
    name = 'crispr_tools',
    version = __version__,
    author = 'John C. Thomas',
    author_email = 'jcthomas000@gmail.com',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 'pandas', 'jacks', 'xlrd', 'shapely'
    ],
    python_requires = '>=3.6',
)
