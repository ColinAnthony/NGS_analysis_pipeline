#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator

__author__ = 'colin'


def line(df, out_string, item, glyc_site, ab_time, bnab_time):
    '''
    :param df: (dataframe) containing the data to plot
    :param out_string: filepath and file name prefix
    :param item: (str) the glycan site to plot
    :param n: (str) the suffix for the file name (excluding file extension)
    :return: None, writes figure to file
    '''

    # force Arial font
    mpl.rc('font', serif='Arial')
    mpl.rcParams['font.family'] = 'Arial'

    outfile = out_string + "_" + glyc_site + "png"
    print("outfile is :", outfile)

    headers = list(df)

    # set axis limits
    xmax = max(df[headers[0]]) + 10
    xmin = 0
    ymax = 100
    ymin = 0

    fig, ax = plt.subplots(1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xticks(list(range(xmin, xmax, 20)), fontsize=14)
    plt.yticks(list(range(ymin, int(ymax), 20)), fontsize=14)
    plt.ylim(ymin, ymax+1)
    plt.xlim(xmin, xmax)

    # set the axis labels
    plt.xlabel("Weeks Post Infection", fontsize=16, labelpad=10)
    plt.ylabel("Frequency of " + str(item) + " (%)", fontsize=16, labelpad=10)

    # plot the data
    plt.plot(df[headers[0]], df[item], linewidth=1, c="#cd6155", marker="o", markersize=4, zorder=1)

    # add intervention/antibody time annotations
    if ab_time is not None:
        plt.text(ab_time, ymax, ' ssNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
        ax.axvline(x=ab_time, color='black', ls='dotted', lw=1, zorder=1)

    if bnab_time is not None:
        plt.text(bnab_time, ymax, ' bNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
        plt.axvline(x=bnab_time, color='black', ls='dotted', lw=1, zorder=1)

    # set figure size, and write to file
    w = 6.875
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)
    plt.savefig(outfile, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')


def main(infile, outpath, lower, ab_time, bnab_time):

    name = os.path.split(infile)[-1].replace("_glycan_freq.csv", "")
    outfile = os.path.join(outpath, name)

    # get the data into a dataframe
    data = pd.read_csv(infile, sep=',', header=0, parse_dates=True, na_values=[' '])
    df = pd.DataFrame(data)
    ndf = df.fillna(method='ffill')

    headers = list(ndf)

    # plot figure for each site with change greater than "lower" %
    upper = 100 - lower
    for item in headers[1:]:
        # only plot sites which change
        if max(df[item]) > lower and min(df[item]) < upper:
            glyc_site = item.replace(" ", "_")
            line(ndf, outfile, item, glyc_site, ab_time, bnab_time)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots the glycan sites which show variation in frequency > 1 %')
    parser.add_argument('-i', '--infile', type=str,
                        help='The input .csv file', required=True)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=True)
    parser.add_argument('-s', '--sensitivity', type=int, default=1,
                        help='degree of % change (ie: to plot if site changes by 1 % use: -s 1', required=False)
    parser.add_argument('-t', '--ab_time', type=int, required=False,
                        help='The time to mask, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', type=int, required=False,
                        help='The second time point to mask, ie: start of bnnAb "-b 50"')

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    sensitivity = args.sensitivity
    ab_time = args.ab_time
    bnab_time = args.bnab_time

    main(infile, outpath, sensitivity, ab_time, bnab_time)
