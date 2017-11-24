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
import pandas as pd
from math import log
import sys
from itertools import chain


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
    try:
        return -10 * log(sum([10**(q / -10) for q in quals]) / len(quals), 10)
    except ZeroDivisionError:
        return None


def median_qual(quals):
    """Receive the integer quality scores of a read and return the median quality for that read."""
    return np.median(quals)


def format_nums(nums):
    return '\t'.join(["{:>12,.2f}".format(n).rstrip('0').rstrip('.') for n in nums])


def format_names(names):
    return '\t'.join(["{:>12}".format(n) for n in names])


def format_pairs_block(pairs):
    res = ""
    for j, l in enumerate(pairs, start=1):
        res += "{:>30}\t".format(str(j) + ":")
        for i in range(len(l)):
            if i % 2:
                res += (" ({:,.2f})\t".format(float(l[i])))
            else:
                res += ("{:>8,.2f}".format(float(l[i]))).rstrip('0').rstrip('.')
        res += "\n"
    return res


def format_pairs_line(pairs):
    pairs = list(chain.from_iterable(pairs))
    res = "\t"
    for i in range(len(pairs)):
        if i % 2:
            res += (" ({:,.1f})\t".format(float(pairs[i])))
        else:
            res += ("{:>8,.0f}".format(float(pairs[i]))).rstrip('0').rstrip('.')
    res += "\n"
    return res


def get_top_5(df, col):
    return df.sort_values(col, ascending=False).head(5)[["lengths", "quals"]].reset_index(drop=True)


def reads_above_qual(df, qual):
    numberAboveQ = np.sum(df["quals"] > qual)
    return numberAboveQ, 100 * (numberAboveQ / len(df.index))


def write_stats(datadfs, outputfile, names=[]):
    """Call calculation functions and write stats file.

    This function takes a list of DataFrames,
    and will create a column for each in the tab separated output.
    """
    if outputfile == 'stdout':
        output = sys.stdout
    else:
        output = open(outputfile, 'wt')
    output.write("General summary\t" + " " * 15 + "{}\n".format(format_names(names)))
    output.write("{:>30}\t{}\n".format("Number of reads:",
                                       format_nums([len(df) for df in datadfs])))
    output.write("{:>30}\t{}\n".format("Total bases:", format_nums(
        [np.sum(df["lengths"]) for df in datadfs])))
    if "aligned_lengths" in datadfs[0]:
        output.write("{:>30}\t{}\n".format("Total bases aligned:", format_nums(
            [np.sum(df["aligned_lengths"]) for df in datadfs])))
    output.write("{:>30}\t{}\n".format("Median read length:", format_nums(
        [np.median(df["lengths"]) for df in datadfs])))
    output.write("{:>30}\t{}\n".format("Mean read length:", format_nums(
        [np.mean(df["lengths"]) for df in datadfs])))
    output.write("{:>30}\t{}\n".format("Read length N50:", format_nums(
        [getN50(np.sort(df["lengths"])) for df in datadfs])))
    if "percentIdentity" in datadfs[0]:
        output.write("{:>30}\t{}\n".format("Average percent identity:", format_nums(
            [np.mean(df["percentIdentity"]) for df in datadfs])))
        output.write("{:>30}\t{}\n".format("Median percent identity:", format_nums(
            [np.median(df["percentIdentity"]) for df in datadfs])))
    if "channelIDs" in datadfs[0]:
        output.write("{:>30}\t{}\n".format("Active channels:", format_nums(
            [np.unique(df["channelIDs"]).size for df in datadfs])))
    output.write("\n")

    if "quals" in datadfs[0]:
        output.write("Top 5 longest reads and their mean basecall quality score\n")
        output.write("{}\n".format(
            format_pairs_block([tuple(l) for l in pd.concat(
                [get_top_5(df, "lengths") for df in datadfs],
                axis=1).values])))
        output.write("Top 5 highest mean basecall quality scores and their read lengths\n")
        output.write("{}\n".format(
            format_pairs_block([tuple(l) for l in pd.concat(
                [get_top_5(df, "quals") for df in datadfs],
                axis=1).values])))
        output.write("Number of reads and fraction above quality cutoffs\n")
        for q in range(5, 30, 5):
            output.write("{:>30}{}".format("Q" + str(q) + ":", format_pairs_line(
                [reads_above_qual(df, q) for df in datadfs])))


# To ensure backwards compatilibity, for a while, keeping exposed function names duplicated:
getN50 = get_N50
removeLengthOutliers = remove_length_outliers
aveQual = ave_qual
writeStats = write_stats
