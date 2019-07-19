import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator


__author__ = 'colin anthony'


def transform_df_by_vl(headers, item, df, vl_file, participant):
    """

    :param headers: column header to plot on x axis data (Time (weeks))
    :param df: dataframe with headers
    :param dcn: dictionary of Time (weeks):viral load
    :param vl: (bool) transform frequency data by viral load
    :return: new data frame, with columns for freq in counts and freq in copies of virus
    """

    time = headers[0]
    df_test = df.copy(deep=True)
    df_test.drop(['sequence_id'], inplace=True, axis=1)
    ndf = df_test.groupby(time)[item].value_counts().reset_index(name='count')
    total = df_test.groupby(time).agg({time: 'count'})
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
        participant_vl['Time'] = participant_vl['Time'].astype(int)
        participant_vl['viral_load'] = participant_vl['viral_load'].astype(str).astype(int)
        participant_vl.drop(["participant_id"], inplace=True, axis=1)
        df2 = df1.merge(participant_vl, on=[time])
        df2["frequency"] = df2['count'] / df2['total'] * 100
        df2['freq_viral_copies'] = (df2["frequency"] / 100) * df2['viral_load']
        df2.to_csv("wrangled_df.csv", sep=',')

        return df2
    else:
        return df1


def plot_loop_stats(headers, item, df, outfile, ab_time, bnab_time, av_heads, avdf, vl_file):
    """
    :param headers: x axis header
    :param item: y axis header
    :param df: dataframe containing the data
    :param outfile: path and name of outfile (minus extension)
    :param ab_time: int of time point label 1
    :param bnab_time: int of time point label 2
    :return: writes figures to file
    """
    n = item.replace(" ", "_")
    outname = outfile + "_" + n
    print("outfile is :", outname)

    xmax = int(max(df[headers]) + 10)
    xmin = 0
    if "length" in item:
        print("length")
        ymax = int(max(df[item]) + 5)
    else:
        ymax = int(max(df[item]) + 1)

    if int(min(df[item])) == 0:
        ymin = 0
    else:
        ymin = int(min(df[item]) - 1)

    if vl_file is not None:
        df.sort_values(by=["freq_viral_copies"], inplace=True, ascending=False)
    else:
        df.sort_values(by=["frequency"], inplace=True, ascending=True)

    fig, ax = plt.subplots(1, 1)
    ax.axes.get_yaxis().set_visible(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.get_xaxis().set_minor_locator(mpl.ticker.MultipleLocator(base=5))

    # plt.title(sname, fontsize=18)
    plt.ylim(ymin, ymax)
    plt.xlim(xmin, xmax)
    plt.ylabel(item.replace("_", " "), fontsize=12, labelpad=14)
    plt.xlabel("Weeks post infection", fontsize=12, labelpad=6)
    plt.xticks(list(range(xmin, xmax, 20)), fontsize=10)
    plt.yticks(list(range(ymin, ymax, 2)), fontsize=10)
    ax.get_yaxis().set_minor_locator(mpl.ticker.MultipleLocator(base=1))

    # plot the data
    if vl_file is None:
        ax.scatter(df[headers], df[item], alpha=0.6, s=df["frequency"]*10, edgecolor='black', lw=0.5, zorder=2)
    else:
        ax.scatter(df[headers], df[item], alpha=0.6, s=df["freq_viral_copies"]/100, edgecolor='black', lw=0.5, zorder=2)
        # c=df[item],

    # plot the average and std
    plt.plot(avdf[av_heads[0]], avdf['mean'], linewidth=0.1, alpha=1, markersize=1, c="#000000", marker="o", zorder=3)
    plt.fill_between(avdf[av_heads[0]], avdf[av_heads[1]] - avdf[av_heads[2]], avdf[av_heads[1]] + avdf[av_heads[2]],
                     color="#b2b5ba", alpha=1, zorder=1)
    # add annotations for nAb time points
    if ab_time is not None:
        plt.text(ab_time, ymax, ' ssNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
        ax.axvline(x=ab_time, color='black', ls='dotted', lw=1, zorder=1)

    if bnab_time is not None:
        plt.text(bnab_time, ymax, ' bNAb', horizontalalignment='left', verticalalignment='center', fontsize=10)
        plt.axvline(x=bnab_time, color='black', ls='dotted', lw=1, zorder=1)

    # set fig size and write output
    w = 6.875
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)
    plt.savefig(outname + '.png', ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')


def main(infile, vl_file, name, outpath, ab_time, bnab_time):
    """
    :param infile: csv file with format like that produced by loop_stats.py
    :param name: string of prefix for graph files
    :param outpath: the output path for plots to be created at.
    :param ab_time: int of time point label 1
    :param bnab_time: int of time point label 2
    :param vl: (bool) transform frequency data by viral load
    :return: writes graphs to file depending on what columns are present in csv file
    """
    print(infile)

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    outfile1 = os.path.join(outpath, name)
    mpl.rc('font', serif='Arial')
    mpl.rcParams['font.family'] = 'Arial'
    data = pd.read_csv(infile, sep=',', header=0, parse_dates=True)
    df = pd.DataFrame(data)
    df.fillna(method='ffill', inplace=True)
    headers = list(df)
    participant = name.split("_")[0]

    for item in headers[1:]:
        if item != 'sequence_id':
            print(item)
            ndf = transform_df_by_vl(headers, item, df, vl_file, participant)
            sum_av_dist = df.groupby([headers[0]]).agg({item: ['mean', 'std']})
            av_vals = pd.DataFrame(sum_av_dist.to_records())
            av_heads = list(av_vals)
            av_vals.rename(columns={av_heads[0]: av_heads[0], av_heads[1]: 'mean', av_heads[2]: 'std'}, inplace=True)
            av_vals.to_csv(outfile1 + "_" + str(item) + "_av.csv", sep=",", index=False)
            av_headers = list(av_vals)
            plot_loop_stats(headers[0], item, ndf, outfile1, ab_time, bnab_time, av_headers, av_vals, vl_file)

    print('Done')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots loop stats from csv file (produced by loop_stats.py)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i1', '--infile1', type=str, required=True,
                        help='The input csv file')
    parser.add_argument('-i2', '--infile2', type=str, default=None, required=False,
                        help='The csv file with columns for participant, time point, viral load. Without this flag, '
                             'scatter point sizes will be scaled by frequency only, not adjusted by viral load')
    parser.add_argument('-o', '--outpath', required=True, type=str,
                        help='path for output plots.')
    parser.add_argument('-n', '--name', type=str, required=True,
                        help='The name of the patient: ie: "CAP177"')
    parser.add_argument('-t', '--ab_time', type=int, required=False,
                        help='The time to mask, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', type=int, required=False,
                        help='The second time point to mask, ie: start of bnnAb "-b 50"')

    args = parser.parse_args()
    infile1 = args.infile1
    infile2 = args.infile2
    name = args.name
    outpath = args.outpath
    ab_time = args.ab_time
    bnab_time = args.bnab_time

    main(infile1, infile2, name, outpath, ab_time, bnab_time)
