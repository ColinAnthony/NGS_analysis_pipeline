#!/usr/bin/python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import operator
import os
import sys
import collections
import argparse
import regex
from Bio import SeqIO
from glob import glob
from glob import glob
from smallBixTools import smallBixTools as st


__author__ = 'David Matten'


def main(in_fasta, outpath, indx):
    print("")
    dct = st.fasta_to_dct(in_fasta)
    dct_of_dcts = {}

    # Build a dictionary of dictionaries. Inner dictionary belongs to 1 patient.
    # All sequences for 1 patient end up in 1 dictionary.
    for k, v in dct.items():
        pt_id = k.split("_")[indx]
        if pt_id in dct_of_dcts.keys():
            dct_of_dcts[pt_id][k] = v
        else:
            dct_of_dcts[pt_id] = {k: v}

    majority_variants = {}
    for pt_id, inner_dct in dct_of_dcts.items():
        majority_variants[pt_id] = ""
        tmp_dct = {}
        for k, v in inner_dct.items():
            if v in tmp_dct.keys():
                tmp_dct[v] += 1
            else:
                tmp_dct[v] = 1
        majority_variants[pt_id] = max(tmp_dct.items(), key=operator.itemgetter(1))[0]
    # dictionary of pt_id: sequence
    #print(majority_variants)

    outfile = os.path.join(outpath, "normalized_distances_from_majority_variant.csv")
    with open(outfile, "w") as handle:
        handle.write("participant,Normalised_hamming_distance_adjusted_(changes_per_100_bases),sequence_id\n")

        # for every sequence for a patient, calculate the normalized distance to the majority variant.
        for pt_id, inner_dct in dct_of_dcts.items():
            this_majority = majority_variants[pt_id]
            for k, v in inner_dct.items():
                dist = st.customdist(v, this_majority)
                norm_dist = dist / len(v)
                normadjustdist_perc = round(norm_dist * 100, 2)

                handle.write(",".join([str(x) for x in [pt_id, normadjustdist_perc, k]]) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate genetic distance from a reference sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--in_fasta', type=str, required=True,
                        help='The path to the input fasta file (ie: /path/to/fasta_file.fasta')
    parser.add_argument('-out', '--outpath', type=str, required=True,
                        help='The path to the where the outfile should be written')
    parser.add_argument('-index', '--index', type=int, required=True,
                        help='We split the sequence id by underscore. Which field should we use to group sequences by '
                             'patient? zero indexed')

    args = parser.parse_args()

    inpath = args.in_fasta
    outpath = args.outpath
    indx = args.index

    main(inpath, outpath, indx)
