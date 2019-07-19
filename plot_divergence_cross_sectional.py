import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl


__author__ = 'colin'


def transform_df_by_freq(headers, df):
    """
    :param headers: column header to plot on x axis data (Time (weeks))
    :param df: dataframe with headers
    :return: new data frame, with columns for freq in counts and freq in copies of virus
    """
    pid = headers[0]
    dist = headers[1]

    test_df = df.copy(deep=True)
    test_df.drop(['sequence_id'], inplace=True, axis=1)
    ndf = test_df.groupby(pid)[dist].value_counts().reset_index(name='count')

    total = test_df.groupby(pid).agg({pid: 'count'})
    total.rename(columns={pid: pid, pid: 'total'}, inplace=True)
    total = total.reset_index()
    df1 = ndf.merge(total[[pid, 'total']], on=[pid])
    df1["frequency"] = ndf['count'] / df1['total'] * 100

    return df1


def divergence_plotter(headers, df, outfile):
    """
    :param headers: x axis header
    :param item: y axis header
    :param df: dataframe containing the data
    :param name: sample name ie: CAP177
    :param vl_file:
    :return: prints graphs to file
    """

    print("outfile is:", outfile)

    x_header = headers[0]
    y_header = headers[1]

    ymax = int(max(df[y_header]) * 1.15)
    ymin = 0

    # set axes
    fig, ax = plt.subplots()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.ylim(ymin, ymax)

    ax.scatter(df[x_header], df[y_header], c="#ffcc66", linewidth=0.5, edgecolor='k',
                s=df["frequency"], alpha=0.60)

    ax.set_ylabel('Divergence from major variant (%)', fontsize=12)
    ax.set_xlabel('Placebo', fontsize=12)
    # ax.set(xticklabels=[])
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90, fontsize=6)

    w = 6.875
    h = 4
    f = plt.gcf()
    f.set_size_inches(w, h)

    plt.savefig(outfile, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')


def main(infile, treatment_grp, name, outpath):
    """
    :param infile: csv file with format like that produced by loop_stats.py or calc_divergence.py
    :param name: string of prefix for graph files
    :param vl: (bool) transform frequency data by viral load
    :return: writes graphs to file depending on what columns are present in csv file
    """
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    outfile = os.path.join(outpath, name + "_divergence.png")
    mpl.rc('font', serif='Arial')
    mpl.rcParams['font.family'] = 'Arial'
    df = pd.read_csv(infile, sep=',', header=0, parse_dates=True)
    df.fillna(method='ffill', inplace=True)
    headers = list(df)
    ndf = transform_df_by_freq(headers, df)

    treatments_df = pd.read_csv(treatment_grp, sep=',', header=0, parse_dates=True)
    treatments = treatments_df.set_index('pub_id')['treatment'].to_dict()
    ndf["treatment"] = ndf["participant"].map(treatments)

    ndf_headers = list(ndf)
    divergence_plotter(ndf_headers, ndf, outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots loop stats from csv file (produced by loop_stats.py)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--infile', type=str, required=True,
                        help='The input csv file')
    parser.add_argument('-s', '--treatment_grp', type=str,
                        help='The csv file with the seq id and treatment assignment', required=True)
    parser.add_argument('-n', '--name', type=str, required=True,
                        help='The name of the patient: ie: "CAP177"')
    parser.add_argument('-o', '--outpath', type=str, required=True,
                        help='The path to where the outfile will be written ')

    args = parser.parse_args()
    infile = args.infile
    treatment_grp = args.treatment_grp
    name = args.name
    outpath = args.outpath

    main(infile, treatment_grp, name, outpath)
