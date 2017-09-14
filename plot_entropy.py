#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import argparse
import pandas as pd
from matplotlib import pyplot as plt
import sys
import numpy as np
from matplotlib import ticker
import os
from scipy import arange
__author__ = 'colin'


def manipulate_df(df):
    headers = list(df)
    df = df.pivot(headers[0], headers[1], headers[2])
    ndf = df.transpose()
    return ndf


def cMaper(num):
    i = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap',
         'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r',
         'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r',
         'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples',
         'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r',
         'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Vega10',
         'Vega10_r', 'Vega20', 'Vega20_r', 'Vega20b', 'Vega20b_r', 'Vega20c', 'Vega20c_r', 'Wistia', 'Wistia_r', 'YlGn',
         'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn',
         'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cool', 'cool_r',
         'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth',
         'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',
         'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot',
         'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno',
         'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r',
         'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r',
         'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', 'summer_r', 'terrain', 'terrain_r', 'viridis',
         'viridis_r', 'winter', 'winter_r']

    j = sorted(i, key=str.lower)
    if num is None:
        for k in j:
            print("for colour {0} use index {1}".format(i.index(k), k))

        print("Use one of these colours with the flag '-c num' (for recommended colour map ('hot'), use: -c 124")
        sys.exit()
    else:
        return i[num]


def heatmap(ndf, outfile, c, ab_time, bnab_time, cbar_title, show):

    nh = list(ndf)
    m = ndf.max(axis=None, skipna=None, level=None, numeric_only=None)
    highest = max(m)
    highest = round(highest, 1)
    step = highest/5
    ylabels = list(arange(0, highest, step))
    las = ylabels[-1] + step
    ylabels.append(las)
    fig, ax = plt.subplots()

    heatmap = ax.pcolor(ndf, cmap=c, alpha=1) #alpha=0.9

    fig = plt.gcf()
    fig.set_size_inches(8, 15)
    ax.set_frame_on(False)
    ax.set_yticks(np.arange(ndf.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(ndf.shape[1]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_bottom()

    labels = nh
    ax.set_xticklabels(labels, minor=False, fontsize=6)

    plt.tick_params(axis='both', which='major', labelsize=6)
    ax.set_xticklabels(labels, fontsize='small')
    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=6)

    ax.set_yticklabels(ndf.index, minor=False, fontsize=6)

    plt.xticks(rotation=90)
    ax.grid(False)

    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    if ab_time is not None:
        sn = nh.index(ab_time)
        plt.text(sn, -1.20, 'ssNAb', horizontalalignment='left', verticalalignment='center', fontsize=6)
        ax.axvline(x=sn, color='#ffffff', ls='dotted', lw=2, zorder=1)

    if bnab_time is not None:
        bn = nh.index(bnab_time)
        plt.text(bn, -1.20, 'bNAb', horizontalalignment='left', verticalalignment='center', fontsize=6)
        plt.axvline(x=bn, color='#ffffff', ls='dotted', lw=2, zorder=1)

    cbar = plt.colorbar(heatmap, aspect=50)
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_yticklabels(ylabels)
    cbar.set_label(cbar_title, rotation=270, labelpad=40, fontsize=16)
    tick_locator = ticker.MaxNLocator(nbins=6)
    cbar.locator = tick_locator
    cbar.update_ticks()

    plt.ylabel("Position", fontsize=16, labelpad=24)
    plt.xlabel("Weeks Post Infection", fontsize=16, labelpad=20)

    h = 8
    w = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')
    if show == True:
        plt.show()
    print("end")


def main(infile, outpath, name, colour, ab_time, bnab_time, show):
    print(infile)
    outfile = os.path.join(outpath, name + '.png')

    c = cMaper(colour)

    data = pd.read_csv(infile, index_col=None)
    df = pd.DataFrame(data)
    head = list(df)
    cbar_title = head[2].replace("_", " ")
    ndf = manipulate_df(df)

    heatmap(ndf, outfile, c, ab_time, bnab_time, cbar_title, show)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots the JSD/entropy as a heatmap',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input entropy .csv file from jensen_shannon_entropy', required=True)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the name of your outfile', required=True)
    parser.add_argument('-c', '--colour', default=124, type=int,
                        help='The name of the colour map to use (leave out to get a list of available otions'
                             'recommended colour scheme is "hot" ie:-c 124 )', required=False)
    parser.add_argument('-t', '--ab_time', default=None, type=int, required=False,
                        help='The time to highlight, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', default=None, type=int, required=False,
                        help='The second time point to highlight, ie: start of bnnAb "-b 50"')
    parser.add_argument('-s', '--show', default=False, action='store_true',
                        help='show tree', required=False)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    colour = args.colour
    ab_time = args.ab_time
    bnab_time = args.bnab_time
    show = args.show

    main(infile, outpath, name, colour, ab_time, bnab_time, show)
