#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

__author__ = 'colin (and David Matten: david.matten@uct.ac.za)'


def main(infile, outpath, ab_time, bnab_time):

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
        print("pos: {}".format(pos))
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

        aa_dict = {'A': '#C8C8C8',
                   'R': '#145AFF',
                   'N': '#00DCDC',
                   'D': '#E60A0A',
                   'C': '#E6E600',
                   'Q': '#00DCDC',
                   'E': '#E60A0A',
                   'G': '#EBEBEB',
                   'H': '#8282D2',
                   'I': '#0F820F',
                   'L': '#0F820F',
                   'K': '#145AFF',
                   'M': '#E6E600',
                   'F': '#3232AA',
                   'P': '#DC9682',
                   'S': '#FA9600',
                   'T': '#FA9600',
                   'W': '#B45AB4',
                   'Y': '#3232AA',
                   'V': '#0F820F',
                   '-': '#707070',
                   'X': '#000000'}

        # aa_dict = {
        #         'A': '#99ff99',
        #         'R': '#0000ff',
        #         'N': '#ff6666',
        #         'D': '#660000',
        #         'C': '#ffff66',
        #         'Q': '#ff3333',
        #         'E': '#660000',
        #         'G': '#000000',
        #         'H': '#6666ff',
        #         'I': '#006600',
        #         'L': '#336633',
        #         'K': '#3333cc',
        #         'M': '#cc9933',
        #         'F': '#666633',
        #         'P': '#666666',
        #         'S': '#ff6633',
        #         'T': '#cc3300',
        #         'W': '#666600',
        #         'Y': '#996633',
        #         'V': '#ff99ff',
        #         '-': '#707070',
        #         'X': '#ff00ff'}

        # if the site is not conserved
        if len(df_new.columns) > 2:

            fig, ax = plt.subplots()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()

            # set the axis limits and axis tick intervals
            plt.yticks(np.arange(y_lim_min, y_lim_max, 20), fontsize=10)
            plt.ylim(0, y_lim_max)
            plt.xticks(rotation=90)

            headers = list(df_new)
            melt_values = headers[:-1]
            plot_df = pd.melt(df_new, id_vars=['participant'], value_vars=melt_values)
            plot_df.rename(columns={"variable": "residue"}, inplace=True)

            sns.set_style("white")
            sns.set_style("ticks")
            sns.axes_style("white")
            # if pos == 454 or pos == 472 or pos == 480:
            #print(plot_df)
            sns.stripplot(x="participant", y="value", data=plot_df, size=20, hue="residue", palette=my_colors, alpha=0.80)
            sns.despine(left=False, bottom=False)

            # plot axis labels
            plt.xlabel("Participant", fontsize=16, labelpad=10)
            plt.ylabel("Frequency (%)", fontsize=16, labelpad=12)
            plt.title("Position " + str(int(pos)), ha='center', fontsize=20)

            # add figure legend
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

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
            out_fn = os.path.join(outpath, "Site_" + str(pos) + "_rel.png")
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
