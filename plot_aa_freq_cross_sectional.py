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

import seaborn as sns

__author__ = 'colin (and David Matten: david.matten@uct.ac.za)'


def main(infile, outpath, ab_time, bnab_time):
    sns.set(style="whitegrid", palette="muted")

    print("plotting frequency data.")
    data = pd.read_csv(infile, sep=',', header=0)
    df = pd.DataFrame(data)

    if outpath == None:
        outpath = os.getcwd()

    selection_bound = 1
    x_header = "participant"
    group_col = "HXB2_position"
    y_lim_min = 0
    y_lim_max = 105
    # xmin = 0

    for pos, df_grp in df.groupby(group_col):

        # manipulate the dataframe
        df_grp.drop(group_col, axis=1, inplace=True)
        transp = df_grp.transpose()
        transp = transp[(transp > selection_bound).any(axis=1)]
        transp = transp[(transp < (100.0-selection_bound)).any(axis=1)]
        df_new = transp.transpose()

        if x_header not in df_new.columns:
            df_new[x_header] = df_grp.loc[:, x_header]

        my_colors = ["#C68A3B", "#7394CA", "#579F6A", "#C47BB2", "#D36C6E", "#53A7A6", "#8D909B", "#D6C463", "#A98B6D",
                     "grey"]

        # if the site is not conserved
        if len(df_new.columns) > 2:
            df_new["participant"] = df_new["participant"].astype(str)

            fig, ax = plt.subplots(1, 1)

            plt.axis('off')
            plt.axes(frameon=False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # ax.spines['bottom'].set_visible(True)
            # ax.spines['left'].set_visible(True)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()

            # set the axis limits and axis tick intervals
            names_list = [str(x) for x in set(df["participant"])]
            plt.yticks(np.arange(y_lim_min, y_lim_max, 20), fontsize=10)
            plt.ylim(0, 100.5)

            y_headers = [x for x in list(df_new)[:-1]]

            num_cols = len(y_headers)
            # x_tick = list(range(1, len(names_list) + 1))
            # plt.xticks(x_tick)
            # ax.set_xticks(range(1, num_cols + 1))
            # plt.axes().set_xticklabels(names_list, rotation='vertical', fontsize=10)


            # if pos == 428:
            #     rawr = pd.melt(df_new, id_vars=['participant'], value_vars=['E', 'K'])
            #     print(rawr)
            # sns.swarmplot(x="participant", y="value", hue="variable", data=rawr)
            # plt.show()
            # r = input("next?")

            df_new.plot(x="participant", title="Position " + str(int(pos)), color=my_colors, marker='o',
                        markersize=14, linewidth=0, zorder=3)

            # plot axis labels
            plt.xlabel("Participant", fontsize=16, labelpad=10)
            plt.ylabel("Frequenecy", fontsize=16, labelpad=12)
            plt.title("Position " + str(int(pos)), ha='center', fontsize=20)
            # add figure legend
            plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., borderpad=1, labelspacing=0.5,
                       handletextpad=0)

            # plot the antibody annotations
            if ab_time is not None:
                plt.text(ab_time, y_lim_max, ' ssNAb', horizontalalignment='left', verticalalignment='center',
                         fontsize=10)
                plt.axvline(x=ab_time, color='black', ls='dotted', lw=1, zorder=2)

            if bnab_time is not None:
                plt.text(bnab_time, y_lim_max, ' bNAb', horizontalalignment='left', verticalalignment='center',
                         fontsize=10)
                plt.axvline(x=bnab_time, color='black', ls='dotted', lw=1, zorder=1)

            # set figure size
            f = plt.gcf()
            w = 6.875
            h = 4
            f.set_size_inches(w, h)

            # write figure to file
            out_fn = os.path.join(outpath, "Site_" + str(int(pos)) + "_rel.png")
            plt.savefig(out_fn, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')
            plt.close('all')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot amino acid frequencies for each position, from a .csv file'
                                                 'format output from calc_entropy_aa_glycan_freq.py.')

    parser.add_argument('-i', '--infile', type=str,
                        help='The input .csv file', required=True)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=True)
    parser.add_argument('-t', '--ab_time', default=None, type=int, required=False,
                        help='The time to mask, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', default=None, type=int, required=False,
                        help='The second time point to mask, ie: start of bnnAb "-b 50"')

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath

    ab_time = args.ab_time
    bnab_time = args.bnab_time

    main(infile, outpath, ab_time, bnab_time)
