# wdecoster
'''
This module provides a few simple math and statistics functions
for other scripts processing Oxford Nanopore sequencing data




## FUNCTIONS
* Calculate read N50 from a set of lengths
get_N50(readlenghts)
* Remove extreme length outliers from a dataset
remove_length_outliers(dataframe, columname)
* Calculate the average Phred quality of a read
ave_qual(qualscores)
* Write out the statistics report after calling readstats function
write_stats(dataframe, outputname)
* Compute a number of statistics, return a dictionary
calc_read_stats(dataframe)
'''

import numpy as np
from math import log
import sys


def get_N50(readlengths):
    '''
    Get read N50.
    Based on https://github.com/PapenfussLab/Mungo/blob/master/bin/fasta_stats.py
    '''
    return readlengths[np.where(np.cumsum(readlengths) >= 0.5 * np.sum(readlengths))[0][0]]


def remove_length_outliers(df, columnname):
    '''
    Remove records with length-outliers above 3 standard deviations from the median
    '''
    return df[df[columnname] < (np.median(df[columnname]) + 3 * np.std(df[columnname]))]


def ave_qual(quals):
    '''
    Calculation function:
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    '''
    return -10 * log(sum([10**(q / -10) for q in quals]) / len(quals), 10)


def median_qual(quals):
    '''
    Calculation function:
    Receive the integer quality scores of a read and return the median quality for that read
    '''
    return np.median(quals)


def calc_read_stats(datadf):
    '''
    For an array of readlengths, return a dictionary containing:
    - the number of reads
    - the total number of bases sequenced
    - the median length
    - the mean length
    - the top 5 longest reads and their quality
    - the maximum average basecall quality
    - the fraction and number of reads above > Qx
    (use a set of cutoffs depending on the observed quality scores)
    - the number of active channels
    '''
    readlengths = np.array(datadf["lengths"])

    res = dict()
    if "runIDs" in datadf:
        res["runIDs"] = list(np.unique(datadf["runIDs"]))
    if "channelIDs" in datadf:
        res["ActiveChannels"] = np.unique(datadf["channelIDs"]).size
    res["NumberOfReads"] = readlengths.size
    if res["NumberOfReads"] < 10:
        return None
    if "percentIdentity" in datadf:
        res["ave-pID"] = np.mean(datadf["percentIdentity"])
        res["med-pID"] = np.median(datadf["percentIdentity"])
    if "aligned_lengths" in datadf:
        res["TotalAlignedBases"] = np.sum(datadf["aligned_lengths"])
    res["TotalBases"] = np.sum(readlengths)
    res["MedianLength"] = np.median(readlengths)
    res["MeanLength"] = np.mean(readlengths)
    indices_top_L = np.argpartition(readlengths, -5)[-5:]
    if "quals" in datadf:
        qualities = np.array(datadf["quals"])
        res["MaxLengthsAndQ"] = [list(e)
                                 for e in zip(
            list(readlengths[indices_top_L]),
            list(qualities[indices_top_L]))]
        indices_top_Q = np.argpartition(qualities, -5)[-5:]
        res["MaxQualsAndL"] = [list(e)
                               for e in zip(
            list(readlengths[indices_top_Q]),
            list(qualities[indices_top_Q]))]
        qualgroups = [q for q in range(5, 50, 5) if q < np.amax(qualities) + 5]
        res["QualGroups"] = dict()
        for q in qualgroups:
            numberAboveQ = np.sum(qualities > q)
            res["QualGroups"][q] = (numberAboveQ, numberAboveQ / res["NumberOfReads"])
    return res


def write_stats(datadf, outputfile):
    '''
    Call calculation function and write to file
    '''
    stat = calc_read_stats(datadf)
    if stat:
        if outputfile == 'stdout':
            output = sys.stdout
        else:
            output = open(outputfile, 'wt')
        output.write("General summary\n")
        output.write("Number of reads:\t{:,}\n".format(stat["NumberOfReads"]))
        output.write("Total bases:\t{:,}\n".format(stat["TotalBases"]))
        if "TotalAlignedBases" in stat:
            output.write("Total bases aligned:\t{:,}\n".format(stat["TotalAlignedBases"]))
        output.write("Median read length:\t{:,}\n".format(stat["MedianLength"]))
        output.write("Mean read length:\t{:,}\n".format(round(stat["MeanLength"], 2)))
        output.write("Readlength N50:\t{:,}\n".format(
            getN50(np.sort(datadf["lengths"]))))
        output.write("\n")
        if "quals" in datadf:
            output.write("Top 5 read lengths and their average basecall quality score\n")
            for length, qual in sorted(stat["MaxLengthsAndQ"], key=lambda x: x[0], reverse=True):
                output.write("Length: {:,}bp\tQ: {}\n".format(length, round(qual, 2)))
            output.write("\n")
            output.write("Top 5 average basecall quality scores and their read lengths\n")
            for length, qual in sorted(stat["MaxQualsAndL"], key=lambda x: x[1], reverse=True):
                output.write("Length: {:,}bp\tQ: {}\n".format(length, round(qual, 2)))
            output.write("\n")
            output.write("Number of reads and fraction above quality cutoffs\n")
            for q in sorted(stat["QualGroups"].keys()):
                output.write("Q{}:\t{:,} ({}%)\n".format(
                    q, stat["QualGroups"][q][0], round(100 * stat["QualGroups"][q][1], 2)))
        if "ave-pID" in stat:
            output.write("\nAverage percent identity:\t{:0.2f}\n".format(stat["ave-pID"]))
            output.write("Median percent identity:\t{:0.2f}\n".format(stat["med-pID"]))
        if "ActiveChannels" in stat:
            output.write("\nData produced using {:,} active channels.\n".format(
                stat["ActiveChannels"]))
        if "runIDs" in stat:
            output.write("\nData was produced in run(s) with ID:\n{}\n".format(
                "\n".join(stat["runIDs"])))
    else:
        sys.stderr.write("Number of reads too low for meaningful statistics.\n")


# To ensure backwards compatilibity, for a while, keeping exposed function names duplicated:
getN50 = get_N50
removeLengthOutliers = remove_length_outliers
aveQual = ave_qual
readstats = calc_read_stats
writeStats = write_stats
