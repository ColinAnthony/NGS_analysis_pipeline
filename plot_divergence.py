#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
import matplotlib.cm as cm
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from pprint import pprint


__author__ = 'colin'


def transform_df_by_vl(headers, df, vl_file, participant):
    '''
    :param headers: column header to plot on x axis data (Time (weeks))
    :param df: dataframe with headers
    :param dcn: dictionary of Time (weeks):viral load
    :param vl: (bool) transform frequency data by viral load
    :return: new data frame, with columns for freq in counts and freq in copies of virus
    '''
    time = headers[0]
    item = headers[1]

    test_df = df.copy(deep=True)
    test_df.drop(['sequence_id'], inplace=True, axis=1)
    ndf = test_df.groupby(time)[item].value_counts().reset_index(name='count')

    total = test_df.groupby(time).agg({time: 'count'})
    total.rename(columns={time: time, time: 'total'}, inplace=True)
    total = total.reset_index()
    df1 = ndf.merge(total[[time, 'total']], on=[time])
    df1["frequency"] = ndf['count'] / df1['total'] * 100

    if vl_file is not None:
        # read in the viral load data
        vl_data = pd.read_csv(vl_file, sep=',', header=0, parse_dates=True)
        vl_df = pd.DataFrame(vl_data)
        vl_df.fillna(method='ffill', inplace=True)

        # get only those rows for the participant
        grouped_vl = vl_df.groupby("participant_id")
        participant_vl = grouped_vl.get_group(participant)
        participant_vl['Time'] = participant_vl.loc[:, 'Time'].astype(int)
        participant_vl['viral_load'] = participant_vl.loc[:, 'viral_load'].astype(str).astype(int)
        participant_vl.drop(["participant_id"], inplace=True, axis=1)
        df2 = df1.merge(participant_vl, on=[time])
        df2["frequency"] = df2['count'] / df2['total'] * 100
        df2['freq_viral_copies'] = (df2["frequency"] / 100) * df2['viral_load']
        df2.to_csv("wrangled_df.csv", sep=',')

        return df2
    else:
        return df1


def divergence_plotter(headers, df, name, outpath, ab_time, bnab_time, av_heads, av_df, vl_file):
    '''
    :param headers: x axis header
    :param item: y axis header
    :param df: dataframe containing the data
    :param name: sample name ie: CAP177
    :param vl_file:
    :return: prints graphs to file
    '''
    x_header = headers[0]
    y_header = headers[1]

    n = y_header.replace(" ", "_")
    if "_" in name:
        title = name.replace("_", " ")
    else:
        title = name

    outname = name + "_" + n + ".png"
    outfile = os.path.join(outpath, outname)
    print("outfile is:", outfile)

    ymax = int(max(df[y_header]) * 1.05)
    xmax = int(max(df[x_header]) + 20)
    ymin = 0
    xmin = 0
    # df.sort_values(by=[y_header], inplace=True, ascending=True)
    if vl_file is not None:
        df.sort_values(by=["freq_viral_copies"], inplace=True, ascending=True)
    else:
        df.sort_values(by=["frequency"], inplace=True, ascending=True)
    # set axes
    fig, ax = plt.subplots(1, 1)
    ax.axes.get_yaxis().set_visible(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_facecolor('white')

    # plot the data
    if vl_file is None:
        ax.scatter(df[x_header], df[y_header], alpha=0.6, s=df["frequency"]*20, edgecolor='black', lw=0.5, zorder=2)
    else:
        ax.scatter(df[x_header], df[y_header], alpha=0.6, s=df["freq_viral_copies"]/10, edgecolor='black', lw=0.5,
                   zorder=2)

    ## add legend
    # c=df[item],
    # add figure legend
    # Maj_var_patch = mlines.Line2D([], [], color='#b15a29', marker='.',
    #                       markersize=15, label='Major Variant')
    # Min_var_patch = mlines.Line2D([], [], color='#97b9cc', marker='.',
    #                               markersize=15, label='Minor Variant')
    # plt.legend(handles=[Maj_var_patch, Min_var_patch], bbox_to_anchor=(1.1, 1.05), frameon=False)

    # add average line
    ax.scatter(av_df[av_heads[0]], av_df[av_heads[1]], alpha=1, s=10, c="#000000", lw=0.1, zorder=3)
    plt.plot(av_df[av_heads[0]], av_df[av_heads[1]], linewidth=1.0, c="#000000", zorder=3)
    plt.fill_between(av_df[av_heads[0]], av_df[av_heads[1]] - av_df[av_heads[2]], av_df[av_heads[1]] + av_df[av_heads[2]],
                     color="#b2b5ba", alpha=1, zorder=1)

    #plt.title(title, fontsize=18)
    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax)
    plt.yticks(list([x / 10 for x in range(int(ymin * 10), int(ymax * 10), 20)]), fontsize=14)
    ax.get_yaxis().set_minor_locator(mpl.ticker.MultipleLocator(base=5))

    plt.ylabel("Divergence (%)", fontsize=24, labelpad=14)
    plt.xlabel("Weeks post infection", fontsize=24, labelpad=12)
    plt.xticks(list(range(xmin,xmax, 20)), fontsize=14)

    # add annotations for nAb time points
    if ab_time is not None:
        plt.text(ab_time, ymax, ' ssNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
        ax.axvline(x=ab_time, color='black', ls='dotted', lw=1, zorder=1)

    if bnab_time is not None:
        plt.text(bnab_time, ymax, ' bNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
        plt.axvline(x=bnab_time, color='black', ls='dotted', lw=1, zorder=1)

    w = 6.875
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')


def main(infile, vl_file, name, ab_time, bnab_time, outpath):
    '''
    :param infile: csv file with format like that produced by loop_stats.py or calc_divergence.py
    :param name: string of prefix for graph files
    :param vl: (bool) transform frequency data by viral load
    :return: writes graphs to file depending on what columns are present in csv file
    '''

    mpl.rc('font', serif='Arial')
    mpl.rcParams['font.family'] = 'Arial'
    data = pd.read_csv(infile, sep=',', header=0, parse_dates=True)
    df = pd.DataFrame(data)
    df.fillna(method='ffill', inplace=True)
    headers = list(df)
    participant = name.split("_")[0]

    ndf = transform_df_by_vl(headers, df, vl_file, participant)

    outfile1 = os.path.join(outpath, name + "_av_distance.csv")
    sum_av_dist = df.groupby([headers[0]]).agg({headers[1]: ['mean', 'std']})
    av_vals = pd.DataFrame(sum_av_dist.to_records())
    av_heads = list(av_vals)
    av_vals.rename(columns={av_heads[0]: av_heads[0], av_heads[1]: 'mean', av_heads[2]: 'std'}, inplace=True)
    av_vals.to_csv(outfile1 + "_" + str(headers[1]) + "_av_loop_stats.csv", sep=",", index=False)
    av_headers = list(av_vals)

    divergence_plotter(headers, ndf, name, outpath, ab_time, bnab_time, av_headers, av_vals, vl_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots loop stats from csv file (produced by loop_stats.py)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i1', '--infile1', type=str, required=True,
                        help='The input csv file')
    parser.add_argument('-i2', '--infile2', type=str, default=None, required=False,
                        help='The csv file with columns for participant, time point, viral load. Without this flag, '
                             'scatter point sizes will be scaled by frequency only, not adjusted by viral load')
    parser.add_argument('-n', '--name', type=str, required=True,
                        help='The name of the patient: ie: "CAP177"')
    parser.add_argument('-t', '--ab_time', type=int, required=False,
                        help='The time to mask, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', type=int, required=False,
                        help='The second time point to mask, ie: start of bnnAb "-b 50"')
    parser.add_argument('-o', '--outpath', type=str, required=True,
                        help='The path to where the outfile will be written ')

    args = parser.parse_args()
    infile1 = args.infile1
    infile2 = args.infile2
    name = args.name
    ab_time = args.ab_time
    bnab_time = args.bnab_time
    outpath = args.outpath

    main(infile1, infile2, name, ab_time, bnab_time, outpath)
