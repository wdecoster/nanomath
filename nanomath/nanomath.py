# wdecoster
import numpy as np
import pandas as pd
from .version import __version__


def getN50(readlengths):
	'''
	Get read N50.
	Based on https://github.com/PapenfussLab/Mungo/blob/master/bin/fasta_stats.py
	'''
	return readlengths[np.where(np.cumsum(readlengths)>=0.5*np.sum(readlengths))[0][0]]


def removeLengthOutliers(df, columnname):
	'''
    Remove records with length-outliers above 3 standard deviations from the median
    '''
	return df[df[columnname] < (np.median(df[columnname]) + 3 * np.std(df[columnname]))]


def aveQual(quals):
	'''	Calculation function: Receive the integer quality scores of a read and return the average quality for that read'''
	return sum(quals) / len(quals)


def readstats(readlengths, qualities):
	'''
	For an array of readlengths, return a dictionary containing:
	- the number of reads
	- the total number of bases sequenced
	- the median length
	- the mean length
	- the top 5 longest reads and their quality
	- the maximum average basecall quality
	- the fraction and number of reads above > Qx (use a set of cutoffs depending on the observed quality scores)
	'''
	res = dict()
	res["NumberOfReads"] = readlengths.size
	res["TotalBases"] = np.sum(readlengths)
	res["MedianLength"] = np.median(readlengths)
	res["MeanLength"] = np.mean(readlengths)
	indices_top_L = np.argpartition(readlengths, -5)[-5:]
	res["MaxLengthsAndQ"] = [zip(list(readlengths[indices_top]),  list(qualities[indices_top]))]
	indices_top_Q = np.argpartition(qualities, -5)[-5:]
	res["MaxQualsAndL"] = [zip(list(readlengths[indices_top]),  list(qualities[indices_top]))]
	qualgroups = [q for q in range(5,50,5) if q < res["MaxQual"] + 5]
	res["QualGroups"] = dict()
	for q in qualgroups:
		numberAboveQ = np.sum(readlengths > q)
		res["QualGroups"][q] = (numberAboveQ, numberAboveQ / res["NumberOfReads"])
