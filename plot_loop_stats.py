#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator

__author__ = 'colin.anthony001@gmail.com'


def returnvldict(patient):
    patient = patient.lower()
    vldict = {'cap008': {'2': 207000,
                        '3': 373000,
                        '4': 647000,
                        '5': 316000,
                        '6': 265000,
                        '8': 373000,
                        '10': 368000,
                        '13': 188000,
                        '15': 161000,
                        '19': 73100,
                        '22': 137000,
                        '25': 98400,
                        '30': 46500,
                        '37': 27600,
                        '56': 39300,
                        '69': 42700,
                        '82': 18700,
                        '94': 25300,
                        '106': 47700,
                        '119': 47800,
                        '132': 17700,
                        '144': 23900,
                        '158': 40700,
                        '171': 212000, },
                  'cap088': {'5': 75500,
                    '6': 61300,
                    '7': 93400,
                    '8': 38700,
                    '11': 12800,
                    '13': 16500,
                    '15': 50500,
                    '17': 21100,
                    '22': 45300,
                    '26': 64300,
                    '30': 180000,
                    '34': 135000,
                    '38': 25400,
                    '42': 38700,
                    '46': 101000,
                    '54': 40300,
                    '57': 136000,
                    '68': 23900,
                    '81': 72400,
                    '94': 84600,
                    '108': 123000,
                    '121': 97000,
                    '135': 142000,
                    '146': 177000,
                    '158': 108000,
                    '173': 330000,
                    '187': 735000,
                    '199': 426000,
                    '213': 314000,
                    '226': 3670000,
                    '238': 4480000,
                    '253': 984,
                    '264': 55, },
                'cap177': {'2': 359000,
                    '4': 698000,
                    '6': 162000,
                    '7': 98800,
                    '9': 60200,
                    '11': 57600,
                    '13': 50300,
                    '19': 89700,
                    '25': 32300,
                    '28': 152000,
                    '32': 109000,
                    '36': 119000,
                    '41': 57200,
                    '46': 26100,
                    '50': 18000,
                    '54': 42100,
                    '67': 108000,
                    '80': 59200,
                    '93': 19100,
                    '107': 8780,
                    '120': 52100,
                    '133': 43800,
                    '146': 63300,
                    '159': 48300,
                    '172': 41300, },
                'cap228': {'7': 3840,
                    '8': 2360,
                    '9': 2600,
                    '10': 4560,
                    '11': 899,
                    '13': 1330,
                    '16': 2680,
                    '18': 2920,
                    '22': 1770,
                    '26': 817,
                    '30': 685,
                    '34': 1150,
                    '37': 486,
                    '43': 786,
                    '47': 1520,
                    '52': 1520,
                    '56': 399,
                    '68': 825,
                    '80': 399,
                    '96': 2780,
                    '110': 6600,
                    '122': 2620,
                    '134': 1330,
                    '147': 3160,
                    '161': 14600,
                    '174': 10800,
                    '187': 22300,
                    '203': 10800,
                    '226': 13500,
                    '239': 13000,
                    '256': 34000,
                    '266': 17600,
                    '278': 10715,
                    '290': 14078,
                    '304': 21938,
                    '316': 29853,
                    '331': 35158,
                    '343': 62816,
                    '357': 218962,
                    '370': 228027,
                    '377': 102110,
                    '381': 1860,
                    '389': 347,
                    '402': 55,
                    '428': 'NaN',
                    '455': 19,
                    '482': 89406,
                    '496': 65297,
                    '507': 1302,
                    '520': 217,
                    '534': 19,
                    '553': 24,
                    '581': 25, },
                'cap255': {'8': 196000,
                    '9': 170000,
                    '10': 244000,
                    '11': 367000,
                    '13': 96800,
                    '15': 106000,
                    '17': 106000,
                    '19': 87400,
                    '23': 49200,
                    '27': 59500,
                    '31': 57800,
                    '36': 59300,
                    '39': 11000,
                    '43': 54100,
                    '47': 74500,
                    '51': 18200,
                    '55': 36400,
                    '68': 29000,
                    '80': 23600,
                    '93': 22200,
                    '106': 5580,
                    '120': 7910,
                    '135': 5110,
                    '149': 11500,
                    '163': 7760,
                    '176': 8830,
                    '188': 15700,
                    '202': 0,
                    '272': 0,
                    '298': 0,
                    '324': 0,
                    '353': 0,
                    '377': 0,
                    '403': 0, },
                'cap256': {'6': 56500,
                    '7': 116000,
                    '8': 116000,
                    '9': 116000,
                    '11': 105000,
                    '13': 51600,
                    '15': 2390000,
                    '17': 141000,
                    '23': 751000,
                    '27': 750000,
                    '30': 626000,
                    '34': 425000,
                    '38': 310000,
                    '42': 494000,
                    '48': 223000,
                    '53': 178000,
                    '59': 257000,
                    '69': 186000,
                    '82': 150000,
                    '94': 55000,
                    '106': 54200,
                    '119': 26200,
                    '134': 52800,
                    '145': 21200,
                    '159': 55000,
                    '176': 102000,
                    '193': 26800,
                    '206': 18500,
                    '227': 53500,
                    '244': 0,
                    '259': 0,
                    '285': 0,
                    '314': 0,
                    '340': 0,
                    '370': 0,
                    '392': 0, },
                'cap257': {'7': 276000,
                    '9': 306000,
                    '10': 269000,
                    '12': 173000,
                    '14': 63400,
                    '16': 32000,
                    '18': 57300,
                    '23': 14800,
                    '27': 16900,
                    '30': 14900,
                    '34': 21200,
                    '38': 10800,
                    '42': 20300,
                    '46': 11600,
                    '50': 8260,
                    '54': 10000,
                    '67': 69300,
                    '80': 26200,
                    '93': 44400,
                    '107': 18100,
                    '122': 6720,
                    '135': 20900,
                    '149': 8230,
                    '161': 15100,
                    '174': 3220,
                    '191': 16800,
                    '201': 102000,
                    '213': 28300,
                    '226': 106000,
                    '240': 52800,
                    '250': 1130,
                    '254': 1029,
                    '275': 0,
                    '301': 0,
                    '327': 0,
                    '353': 0,
                    '378': 0, },
                'cap357': {'7': 14900,
                    '8': 45500,
                    '9': 28100,
                    '10': 9080,
                    '12': 8310,
                    '14': 7220,
                    '17': 3750,
                    '19': 31000,
                    '23': 8490,
                    '28': 7650,
                    '34': 22100,
                    '43': 23300,
                    '48': 44000,
                    '52': 399,
                    '58': 11596,
                    '70': 9166,
                    '87': 2092,
                    '96': 5728,
                    '111': 3058,
                    '123': 3409,
                    '134': 2036,
                    '147': 17222,
                    '162': 15065,
                    '175': 22202,
                    '186': 6343,
                    '199': 7564,
                    '214': 9891,
                    '227': 10622, }
              }

    if str(patient).lower() in vldict.keys():
        return vldict[patient]


def transform_df_by_vl(headers, df, dcn, vl_file):
    '''
    :param headers: column header to plot on x axis data (Time (weeks))
    :param df: dataframe with headers
    :param dcn: dictionary of Time (weeks):viral load
    :return: new data frame, with columns for freq in counts and freq in copies of virus
    '''

    time = headers[0]
    item = headers[1]
    df.drop(['Sequence_ID'], inplace=True, axis=1)
    ndf = df.groupby(time)[item].value_counts().reset_index(name='count')

    total = df.groupby(time).agg({time: 'count'})
    total.rename(columns={time: time, time: 'total'}, inplace=True)
    total = total.reset_index()

    df1 = ndf.merge(total[[time, 'total']], on=[time])
    df1["frequency"] = ndf['count'] / df1['total'] * 100

    if vl_file is not None:
        vl_data = pd.read_csv(vl_file, sep=',', header=0, parse_dates=True)
        vl_df = pd.DataFrame(vl_data)
        ## get only those rows for the participant

        # vl_df = pd.DataFrame.from_dict(dcn, orient='index')
        # vl_df.reset_index(inplace=True)
        # vl_df.rename(columns={'index': 'Time_(weeks)', 0: 'viral_load'}, inplace=True)
        # vl_df['Time_(weeks)'] = vl_df['Time_(weeks)'].astype(int)
        # vl_df.sort_values(by=["Time_(weeks)"], inplace=True, ascending=True)

        df2 = df1.merge(vl_df, on=[time])
        df2["frequency"] = df2['count'] / df2['total'] * 100
        df2['freq_viral_copies'] = (df2["frequency"] / 100) * df2['viral_load']
        df2.to_csv("wrangled_df.csv", sep=',')

        return df2
    else:
        return df1


def bubble(headers, item, df, outfile, ab_time, bnab_time, av_heads, avdf, vl_file):
    '''
    :param headers: x axis header
    :param item: y axis header
    :param df: dataframe containing the data
    :param outfile: path and name of outfile (minus extension)
    :param ab_time: int of time point label 1
    :param bnab_time: int of time point label 2
    :return: writes figures to file
    '''
    n = item.replace(" ", "_")
    outname = outfile + "_" + n
    print("outfile is :", outname)

    xmax = int(max(df[headers]) + 10)
    xmin = 0
    ymax = int(max(df[item]) + 1)
    if int(min(df[item])) == 0:
        ymin = 0
    else:
        ymin = int(min(df[item]) - 1)

    df.sort_values(by=[item], inplace=True, ascending=True)
    df.sort_values(by=["frequency_copies"], inplace=True, ascending=False)
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
        ax.scatter(df[headers], df[item], alpha=0.8, s=df["frequency"]/100, edgecolor='black', lw=0.5, zorder=2)
    else:
        ax.scatter(df[headers], df[item], alpha=0.6, s=df["freq_viral_copies"]*400, edgecolor='black', lw=0.5, zorder=2)
        # c=df[item],

    # plot the average and std
    plt.plot(avdf[av_heads[0]], avdf['mean'], linewidth=0.1, markersize=1, c="#000000", marker="o", zorder=3)
    plt.fill_between(avdf[av_heads[0]], avdf[av_heads[1]] - avdf[av_heads[2]], avdf[av_heads[1]] + avdf[av_heads[2]],
                     color="#000000", alpha=0.1, zorder=1)
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
    '''
    :param infile: csv file with format like that produced by loop_stats.py
    :param name: string of prefix for graph files
    :param outpath: the output path for plots to be created at.
    :param ab_time: int of time point label 1
    :param bnab_time: int of time point label 2
    :param vl: (bool) transform frequency data by viral load
    :return: writes graphs to file depending on what columns are present in csv file
    '''
    print(infile)

    outfile1 = os.path.join(outpath, name)
    mpl.rc('font', serif='Arial')
    mpl.rcParams['font.family'] = 'Arial'
    data = pd.read_csv(infile, sep=',', header=0, parse_dates=True)
    df = pd.DataFrame(data)
    headers = list(df)
    d = returnvldict(name.split("_")[0])

    for item in headers[1:]:
        ndf = transform_df_by_vl(headers[0], df, d, vl_file)
        sumavdist = df.groupby([headers[0]]).agg({item: ['mean', 'std']})
        avs = pd.DataFrame(sumavdist.to_records())
        av_heads = list(avs)
        avs.rename(columns={av_heads[0]: av_heads[0], av_heads[1]: 'mean', av_heads[2]: 'std'}, inplace=True)
        avs.to_csv(outfile1 + "_" + str(item) + "_av_loop_stats.csv", sep=",", index=False)
        av_h = list(avs)
        bubble(headers[0], item, ndf, outfile1, ab_time, bnab_time, av_h, avs, vl_file)

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
