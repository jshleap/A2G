from setuptools import setup, find_packages
from A2G.__version__ import version
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with open(path.join(this_directory, 'requirements.txt'), encoding='utf-8'
          ) as r:
    requirements = r.read().strip().split()

setup(
    name='A2G',
    version=version,
    packages=find_packages() + ['A2G'],
    url='https://github.com/jshleap/A2G',
    license='GNU v3',
    author='jshleap',
    author_email='jshleap@gmail.com',
    scripts=['bin/A2G', 'bin/A2G_mpi'],
    description='Accurarate amplicon alignment to gene consensus',
    python_requires='>=3.6',
    install_requires=[requirements],
    long_description=long_description,
    long_description_content_type='text/markdown'
)