#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
import matplotlib.cm as cm


__author__ = 'colin (and David Matten: david.matten@uct.ac.za)'


def main(infile, outpath, group_col, time_header, selection_bound, y_lim_min, y_lim_max, ab_time, bnab_time):
    print("plotting frequency data.")
    df = pd.read_csv(infile, sep=',', header=0)
    if outpath == None:
        outpath = os.getcwd()
    try:
        y_lim_max = float(y_lim_max)
    except TypeError as e:
        print(e)
        print("Could not convert y limit max to a float.\nExiting.")
        sys.exit()
    try:
        y_lim_min = float(y_lim_min)
    except TypeError as e:
        print(e)
        print("Could not convert y limit min to a float.\nExiting.")
        sys.exit()
    try:
        selection_bound= float(selection_bound)
    except TypeError as e:
        print(e)
        print("Could not convert selection bound to a float.\nExiting.")
        sys.exit()


    for pos, df_grp in df.groupby(group_col):
        df_grp.drop(group_col, axis=1, inplace=True)

        transp = df_grp.transpose()

        transp = transp[(transp > selection_bound).any(axis=1)]
        transp = transp[(transp < (100.0-selection_bound)).any(axis=1)]
        df_new = transp.transpose()
        # print(df_new)
        if time_header not in df_new.columns:
            df_new[time_header] = df_grp[time_header]

        my_colors = ["#C68A3B", "#7394CA", "#579F6A", "#C47BB2", "#D36C6E", "#53A7A6", "#8D909B", "#D6C463", "#A98B6D",
                     "grey"]
                        #["#d73027", "#4575b4",  "#5aae61",  "#fdae61",  "#fee090",  "#abd9e9",  "#b2abd2",
                        # "#80cdc1",  "#de77ae",  "#c2a5cf"]

        if len(df_new.columns) > 2:
            xmax = max(df_new[time_header])

            fig, ax = plt.subplots(1, 1)
            plt.axis('off')
            ax.set_axis_off()
            plt.axes(frameon=False)
            ax.axes.get_xaxis().set_visible(True)
            ax.axes.get_yaxis().set_visible(True)
            ax.set_facecolor('white')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(True)
            ax.spines['left'].set_visible(True)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            ax.get_xaxis().set_minor_locator(mpl.ticker.MultipleLocator(base=5))

            # plot the aa freq for the given position
            df_new.plot(x=time_header, title="Position " + str(int(pos)), color=my_colors, marker='o',
                        markersize=5, linewidth=2, zorder=3)
            # plt.minorticks_on()
            fig.patch.set_visible(False)
            # ax.axis('off')
            # ax.grid(True, which='both')
            plt.xticks(list(range(0, int(xmax), 20)), fontsize=10)
            plt.yticks(np.arange(0, y_lim_max, 20), fontsize=10)
            # plt.yticks(list(np.linspace(0.0, xmax, num=6)), fontsize=10)
            plt.ylim(0, 100.05)
            # plt.yticks(np.arange(0, 1.1), fontsize=10)
            plt.grid(b=True, which='major', color='w', linewidth=1.0)
            plt.grid(b=True, which='minor', color='w', linewidth=0.5)

            # plot the antibody annotations
            if ab_time is not None:
                plt.text(ab_time, ymax, ' ssNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
                plt.axvline(x=ab_time, color='black', ls='dotted', lw=1, zorder=2)

            if bnab_time is not None:
                plt.text(bnab_time, ymax, ' bNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
                plt.axvline(x=bnab_time, color='black', ls='dotted', lw=1, zorder=1)

            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.get_xaxis().set_minor_locator(mpl.ticker.MultipleLocator(base=2))
            ax.get_yaxis().set_minor_locator(mpl.ticker.MultipleLocator(base=2))
            plt.tick_params(axis='x', pad=6, direction="out", length=4, width=1.0, which='major')
            plt.tick_params(axis='y', pad=6, direction="out", length=4, width=1.0, which='major')
            plt.tick_params(axis='x', pad=2, direction="out", length=2, width=0.5, which='minor')
            plt.tick_params(axis='y', pad=2, direction="out", length=2, width=0.5, which='minor')

            plt.xlabel("Weeks Post Infection", fontsize=16, labelpad=10)
            plt.ylabel("Frequenecy", fontsize=16, labelpad=12)
            plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., borderpad=1, labelspacing=0.5,
                       handletextpad=0)
            plt.title("Position " + str(int(pos)), ha='center', fontsize=20)
            out_fn = os.path.join(outpath, "Site_" + str(int(pos)) + "_rel.png")
            f = plt.gcf()
            f.set_size_inches(15.0, 10.0)
            plt.savefig(out_fn, dpi=300)#, facecolor='white')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot amino acid frequencies in each position, from an input .csv file'
                                                 'described format in the repo readme.')
    parser.add_argument('-i', '--inpath', type=str,
                        help='The input .csv file', required=False)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=False)
    parser.add_argument('-grp', '--group_col', type=str,
                        help='The heading of the column which should be used to group on. Typically '
                        'this is position., default="Pos"', default="HXB2_position", required=False)
    parser.add_argument('-th', '--time_header', type=str, help='The heading of the column which should be used for '
                        'the x axis of the plots. Typically this time., default="Time"', default="Time", required=False)
    parser.add_argument('-sb', '--selection_bound', type=str, help='The threshold value to use when selecting '
                        'frequencies of amino acids from the input file, for plotting. If a value of X is supplied here'
                        ' then only values between X and 1-X will be included in the plot. So, for more sensitive '
                        'plots, use a lower value (eg: 0.0001) and less sensitive plots, use a larger value (eg: '
                        ' 0.01), default=0.01. This argument should not be specified if absolute values are being used.'
                        , default=1.00, required=False)
    parser.add_argument('-ymin', '--y_limit_min', type=float, default=-0.00, help='The min value on the y axis., '
                        'default=-0.05', required=False)
    parser.add_argument('-ymax', '--y_limit_max', type=float, default=100.05, help='The max value on the y axis., '
                        'default=1.1', required=False)
    parser.add_argument('-t', '--ab_time', type=int, required=False,
                        help='The time to mask, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', type=int, required=False,
                        help='The second time point to mask, ie: start of bnnAb "-b 50"')


    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    group_col = args.group_col
    time_header = args.time_header
    selection_bound = args.selection_bound
    ymin = args.y_limit_min
    ymax = args.y_limit_max
    ab_time = args.ab_time
    bnab_time = args.bnab_time


    main(infile, outpath, group_col, time_header, selection_bound, ymin, ymax, ab_time, bnab_time)
