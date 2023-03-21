# nanomath

This module provides a few simple math and statistics functions for other scripts processing Oxford Nanopore sequencing data

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40wouter_decoster)](https://twitter.com/wouter_decoster)
[![install with conda](https://anaconda.org/bioconda/nanomath/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanomath)
[![install with Debian](https://www.debian.org/logos/button-mini.png)](https://tracker.debian.org/pkg/python-nanomath)

## FUNCTIONS

* Calculate read N50 from a set of lengths `get_N50(readlenghts)`  
* Remove extreme length outliers from a dataset `remove_length_outliers(dataframe, columname)`  
* Calculate the average Phred quality of a read `ave_qual(qualscores)`  
* Write out the statistics report after calling readstats function `write_stats(dataframe, outputname)`  
* Compute a number of statistics, return a dictionary `calc_read_stats(dataframe)`  

As of **v1.3.0**, nanomath calculates the average quality differently, by first converting per-read phred scale averages to error rates, take the average, and converting back ([nanostat#40](<https://github.com/wdecoster/nanostat/issues/40>))

## INSTALLATION

```bash
pip install nanomath
```

or  
[![install with conda](https://anaconda.org/bioconda/nanomath/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanomath)

```
conda install -c bioconda nanomath
```

## CONTRIBUTORS

[@alexomics](https://github.com/alexomics) for fixing the indentation of the printed stats

## CITATION

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty149/4934939).
