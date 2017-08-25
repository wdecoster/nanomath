# nanomath
This module provides a few simple math and statistics functions for other scripts processing Oxford Nanopore sequencing data

[![Twitter URL](https://img.shields.io/twitter/url/https/twitter.com/wouter_decoster.svg?style=social&label=Follow%20%40wouter_decoster)](https://twitter.com/wouter_decoster)
[![install with conda](https://anaconda.org/bioconda/nanomath/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanomath)
[![Build Status](https://travis-ci.org/wdecoster/nanomath.svg?branch=master)](https://travis-ci.org/wdecoster/nanomath)
[![Code Health](https://landscape.io/github/wdecoster/nanomath/master/landscape.svg?style=flat)](https://landscape.io/github/wdecoster/nanomath/master)


## FUNCTIONS
* Calculate read N50 from a set of lengths `getN50(readlenghts)`  
* Remove extreme length outliers from a dataset `removeLengthOutliers(dataframe, columname)`  
* Calculate the average Phred quality of a read `aveQual(qualscores)`  
* Write out the statistics report after calling readstats function `writeStats(dataframe, outputname)`  
* Compute a number of statistics, return a dictionary `readstats(dataframe)`  


## INSTALLATION
```bash
pip install nanomath
```
or  
[![install with conda](https://anaconda.org/bioconda/nanomath/badges/installer/conda.svg)](https://anaconda.org/bioconda/nanomath)
```
conda install -c bioconda nanomath
```

## STATUS
[![Build Status](https://travis-ci.org/wdecoster/nanomath.svg?branch=master)](https://travis-ci.org/wdecoster/nanomath)
[![Code Health](https://landscape.io/github/wdecoster/nanomath/master/landscape.svg?style=flat)](https://landscape.io/github/wdecoster/nanomath/master)
