# nanomath
This module provides a few simple math and statistics functions for other scripts processing Oxford Nanopore sequencing data


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
