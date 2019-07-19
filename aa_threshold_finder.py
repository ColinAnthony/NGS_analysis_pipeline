#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import pandas as pd


__author__ = 'David Matten. Assistance from Colin Anthony'


def get_previous_freq(position_df, aa, wpi, wpis, time_heading):
    wpi_index = wpis.index(wpi)
    if wpi_index == 0:
        previous_wpi = wpis[0]
    else:
        previous_wpi = wpis[wpi_index - 1]
    previous_row = position_df[position_df[time_heading] == previous_wpi]
    previous_aa_freq = previous_row[aa].values[0]
    return previous_aa_freq


def main(infile, outpath):  # , ab_time, bnab_time
    print("Finding time points where amino acid frequencies cross thresholds.")
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    out_name = os.path.split(infile)[-1].replace(".csv", "threshold_kinetics.csv")
    outfile = os.path.join(outpath, out_name)
    df = pd.read_csv(infile, sep=',', header=0)
    print(df)

    time_heading = 'Time'
    position_heading = 'HXB2_position'

    all_col_heads = list(df.columns.values)
    all_AAs = all_col_heads[:]
    all_AAs.remove(position_heading)
    all_AAs.remove(time_heading)  # This includes gaps.

    positions = list(set(list(df['HXB2_position'])))
    positions.sort()

    wpis = list(set(list(df['Time'])))
    wpis.sort()

    threshold = 10

    with open(outfile, "w") as fw:
        fw.write("position,WPI,AA_crossing_threshold_going_UP,AA_freq,this_aa_previously_was,original\n")
        for position in positions:
            position_df = df[df[position_heading] == position]
            earliest_row = position_df[position_df[time_heading] == wpis[0]]

            earliest_row = earliest_row.drop([time_heading, position_heading], axis=1)
            original_max_aa = earliest_row.idxmax(axis=1).values[0]

            for wpi in wpis:
                # print("Considering positions: {}".format(position))
                # print("Considering wpi: {}".format(wpi))

                this_row = df[df['HXB2_position'] == position]
                this_row = this_row[this_row['Time'] == wpi]

                for aa in all_AAs:
                    this_aa_freq = this_row[aa].values[0]
                    this_aa_previous_freq = get_previous_freq(position_df, aa, wpi, wpis, time_heading)
                    if (this_aa_freq > threshold) and (this_aa_previous_freq < threshold):
                        fw.write("{}, {}, {}, {}, {}, {}\n".format(position, wpi, aa, this_aa_freq,
                                                                   this_aa_previous_freq, original_max_aa))

    print("threshold calculations completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='From amino acid frequency .csv files produced from this repositories '
                                                 'calc_entropy_aa_glycan_freq.py, here we find where amino acid '
                                                 'frequencies cross a threshold from one time point to the next.')

    parser.add_argument('-in', '--infile', type=str,
                        help='The input .csv file', required=True)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=True)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath

    main(infile, outpath)
