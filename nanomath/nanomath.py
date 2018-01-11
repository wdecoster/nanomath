# wdecoster
"""
This module provides a few simple math and statistics functions
for other scripts processing Oxford Nanopore sequencing data




# FUNCTIONS
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
"""

import numpy as np
from math import log
import sys


class Stats(object):
    def __init__(self, df, name):
        self.name = name
        self.number_of_reads = len(df)
        self.number_of_bases = np.sum(df["lengths"])
        if "aligned_lengths" in df:
            self.number_of_bases_aligned = np.sum(df["aligned_lengths"])
        self.median_read_length = np.median(df["lengths"])
        self.mean_read_length = np.mean(df["lengths"])
        self.n50 = get_N50(np.sort(df["lengths"]))
        if "percentIdentity" in df:
            self.average_identity = np.mean(df["percentIdentity"])
            self.median_identity = np.median(df["percentIdentity"])
        if "channelIDs" in df:
            self.active_channels = np.unique(df["channelIDs"]).size
        if "quals" in df:
            self.mean_qual = np.mean(df["quals"])
            self.median_qual = np.median(df["quals"])
            self.top5_lengths = get_top_5(df, "lengths")
            self.top5_quals = get_top_5(df, "quals")
            self.reads_above_qual = [reads_above_qual(df, q) for q in range(5, 30, 5)]


def get_N50(readlengths):
    """Calculate read length N50.

    Based on https://github.com/PapenfussLab/Mungo/blob/master/bin/fasta_stats.py
    """
    return readlengths[np.where(np.cumsum(readlengths) >= 0.5 * np.sum(readlengths))[0][0]]


def remove_length_outliers(df, columnname):
    """Remove records with length-outliers above 3 standard deviations from the median."""
    return df[df[columnname] < (np.median(df[columnname]) + 3 * np.std(df[columnname]))]


def ave_qual(quals):
    """Calculate average basecall quality of a read.

    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale

    Return None for ZeroDivisionError
    """
    if quals:
        return -10 * log(sum([10**(q / -10) for q in quals]) / len(quals), 10)
    else:
        return None


def median_qual(quals):
    """Receive the integer quality scores of a read and return the median quality for that read."""
    return np.median(quals)


def get_top_5(df, col):
    res = df.sort_values(col, ascending=False) \
        .head(5)[["lengths", "quals"]] \
        .reset_index(drop=True) \
        .itertuples(index=False, name=None)
    return [str(i) + " (" + str(j) + ")" for i, j in res]


def reads_above_qual(df, qual):
    numberAboveQ = np.sum(df["quals"] > qual)
    return "{} ({})".format(numberAboveQ, 100 * (numberAboveQ / len(df.index)))


def feature_list(stats, feature, index=None):
    if index:
        return '\t'.join([s.__dict__[feature][index] for s in stats])
    else:
        return '\t'.join([s.__dict__[feature] for s in stats])


def write_stats(datadfs, outputfile, names=[]):
    """Call calculation functions and write stats file.

    This function takes a list of DataFrames,
    and will create a column for each in the tab separated output.
    """
    if outputfile == 'stdout':
        output = sys.stdout
    else:
        output = open(outputfile, 'wt')

    stats = [Stats(df, name) for df, name in zip(datadfs, names)]
    features = {
        "General summary": "name",
        "Number of reads": "number_of_reads",
        "Total bases": "number_of_bases",
        "Total bases aligned": "number_of_bases_aligned",
        "Median read length": "median_read_length",
        "Mean read length": "mean_read_length",
        "Read length N50": "n50",
        "Average percent identity": "average_identity",
        "Median percent identity": "median_identity",
        "Active channels": "active_channels",
        "Mean read quality": "mean_qual",
        "Median read quality": "median_qual",
    }
    long_features = {
        "Top 5 longest reads and their mean basecall quality score": "top5_lengths",
        "Top 5 highest mean basecall quality scores and their read lengths": "top5_quals",
        "Number of reads and fraction above quality cutoffs": "reads_above_qual",
    }
    for f in features.keys():
        try:
            output.write("{}:\t{}\n".format(f, feature_list(stats, features[f])))
        except AttributeError:
            pass
    if all(["quals" in df for df in datadfs]):
        for lf in long_features.keys():
            output.write(lf + "\n")
            for i in range(5):
                output.write("{}:\t{}\n".format(i + 1, feature_list(stats, features[f], index=i)))
