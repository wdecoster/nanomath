# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open('nanomath/version.py').read())

setup(
    name='nanomath',
    version=__version__,
    description='A few simple math function for other Oxford Nanopore processing scripts',
    long_description='A few simple math function for other Oxford Nanopore processing scripts.',
    url='https://github.com/wdecoster/nanomath',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore sequencing plotting quality control',
    packages=find_packages(),
    install_requires=['pandas', 'numpy',],
    package_dir={'nanomath': 'nanomath'})
