from __future__ import print_function
from __future__ import division
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


__author__ = 'colin'


def transform_df_by_vl(headers, df, vl_file, participant):
    """
    :param headers: column header to plot on x axis data (Time (weeks))
    :param df: dataframe with headers
    :param dcn: dictionary of Time (weeks):viral load
    :param vl: (bool) transform frequency data by viral load
    :return: new data frame, with columns for freq in counts and freq in copies of virus
    """
    time = headers[0]
    item = headers[1]

    test_df = df.copy(deep=True)

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
        df2 = test_df.merge(participant_vl, on=[time])
        df2['freq_viral_copies'] = df2["Frequency"] * df2['viral_load']
        df2.to_csv("wrangled_df.csv", sep=',')

        return df2
    else:
        return test_df


def divergence_plotter(headers, df, name, outfile, ab_time, bnab_time, vl_file):
    """
    :param headers: x axis header
    :param item: y axis header
    :param df: dataframe containing the data
    :param name: sample name ie: CAP177
    :param vl_file:
    :return: prints graphs to file
    """

    ymax = 101
    xmax = 30
    ymin = 0
    xmin = 0
    selection_bound = 10
    x_step = int(xmax/10.0)
    # fiter the dataframe to only those haplotypes which are greater than the % selection bound and up to time xmax

    # set haplotypes as column headers
    piv_df = df.pivot(index=headers[0], columns=headers[1])

    # reset time as column not index
    piv_df.reset_index(level=0, inplace=True)

    # filter by desited xaxis time point
    filt_df = piv_df[(piv_df.Time <= xmax)]

    # filter by those column with xmax > selection and xmin < 100 - selection
    new_df = filt_df.loc[:, (filt_df.max() >= selection_bound) & (filt_df.min() <= 100 - selection_bound)]
    new_df.fillna(value=0)

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

    new_df.set_index("Time", inplace=True)

    new_df.columns = new_df.columns.get_level_values(1)
    new_df.plot(kind='line', style='.-')

    plt.legend(frameon=False, framealpha=False,bbox_to_anchor=(1.01, 0.901))

    plt.ylim(ymin, ymax + 1)
    plt.xlim(xmin, xmax)
    plt.yticks(list(range(ymin, ymax, 20)), fontsize=14)
    plt.xticks(list(range(xmin, xmax, x_step)), fontsize=14)

    plt.ylabel("Frequency (%)", fontsize=24, labelpad=14)
    plt.xlabel("Weeks post infection", fontsize=24, labelpad=12)

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
    # plt.show()
    plt.savefig(outfile, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')


def main(infile, vl_file, name, ab_time, bnab_time, outpath):
    """

    :param vl_file: csv file with format like that produced by loop_stats.py or calc_divergence.py
    :param ab_time: (int) time point that nAbs emerge
    :param bnab_time: (int) time point that bnAbs emerge
    :param outpath: (str) path to output folder
    :return: writes graphs to file depending on what columns are present in csv file
    """

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    mpl.rc('font', serif='Arial')
    mpl.rcParams['font.family'] = 'Arial'
    data = pd.read_csv(infile, sep=',', header=0, parse_dates=True)
    df = pd.DataFrame(data)
    df.fillna(method='ffill', inplace=True)
    headers = list(df)
    participant = name.split("_")[0]

    ndf = transform_df_by_vl(headers, df, vl_file, participant)
    headers = list(ndf)
    outfile = os.path.join(outpath, name + "_multi_site_site.png")

    # sum_av_dist = df.groupby([headers[0]]).agg({headers[1]: ['mean', 'std']})
    # av_vals = pd.DataFrame(sum_av_dist.to_records())
    # av_heads = list(av_vals)
    # av_vals.rename(columns={av_heads[0]: av_heads[0], av_heads[1]: 'mean', av_heads[2]: 'std'}, inplace=True)
    # av_vals.to_csv(outfile1 + "_" + str(headers[1]) + "_av_loop_stats.csv", sep=",", index=False)
    # av_headers = list(av_vals)

    divergence_plotter(headers, ndf, name, outfile, ab_time, bnab_time, vl_file)


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
