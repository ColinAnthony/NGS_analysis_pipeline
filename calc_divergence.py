#!/usr/bin/python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import os
import sys
import collections
import argparse
from itertools import groupby


__author__ = 'colin'


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for k, v in my_gen:
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = str(v).replace("~", "_")

    return dct


def gethxb2(dict):
    """
    :param dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    """
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dict.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dict[k]
            print("Found hxb2 ref. seq. Its full name is: ", hxb2_key)
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")

    return str(hxb2_key), str(hxb2_seq)


def customdist(s1, s2):

    if len(s1) != len(s2):
        print("sequences must be the same length")
        sys.exit()

    dist = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            dist += 1
    diff = 0
    for i in range(len(s1)-1):
        if s1[i] != s2[i]:
            if (s1[i] == "-" and s1[i+1] == "-" and s2[i+1] != "-") \
                    or (s2[i] == "-" and s2[i+1] == "-" and s1[i+1] != "-"):
                diff += 1

    return (dist-diff)


def normcustomdist(s1, s2):

    if len(s1) != len(s2):
        print("sequences must be the same length")
        sys.exit()

    dist = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            dist += 1
    diff = 0
    for i in range(len(s1)-1):
        if s1[i] != s2[i]:
            if (s1[i] == "-" and s1[i+1] == "-" and s2[i+1] != "-") \
                    or (s2[i] == "-" and s2[i+1] == "-" and s1[i+1] != "-"):
                diff += 1

    dist = dist-diff
    normdist = dist / len(s1)

    return normdist


def main(ref, align_file, outpath, name):

    print(align_file)

    # set the outfile name
    name = name + "_divergence.csv"
    outfile = os.path.join(outpath, name)

    # import the reference
    ref = fasta_to_dct(ref)
    refseq = str(ref.seq).upper()

    # store all sequnces in a dict
    all_sequences = fasta_to_dct(align_file)

    # get hxb2 seq and remove it from the dict
    hxb2_name, hxb2_seq = gethxb2(all_sequences)
    del all_sequences[hxb2_name]

    # write the headings to the outfile
    with open(outfile, "w") as handle:
        handle.write("Time,Normalised_hamming_distance_adjusted_(changes_per_100_bases),sequence_id\n")

    # calculate the divergence from the reference for each sequnce
    for seq_name, seq in sorted(all_sequences.items()):
        if len(seq) != len(refseq):
            print("input sequence and reference sequence were not the same length.")
            sys.exit()
        else:
            time = seq_name.split("_")[2][:-3]

            # hamdist = distance.hamming(seq, refseq, normalized=True)
            normadjustdist_perc = round(normcustomdist(seq, refseq) * 100, 2)

            with open(outfile, "a") as handle:
                handle.write(",".join([str(x) for x in [time, normadjustdist_perc,  seq_name]]) + "\n")


    print("Divergence calculations are complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate genetic distance from a reference sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=str, required=True,
                        help='The path to the fasta file (ie: /path/to/5haplotype')
    parser.add_argument('-r', '--reference', type=str, required=True,
                        help='The fasta file with only the reference sequence. '
                             'Usually taken from the most abundant haplotype at the first time point')
    parser.add_argument('-o', '--outpath', type=str, required=True,
                        help='The path to the where the outfile should be written')
    parser.add_argument('-n', '--name', type=str, required=True,
                        help='The name of the outfile (no extension')
    args = parser.parse_args()

    ref = args.reference
    infile = args.infile
    outpath = args.outpath
    name = args.name

    main(ref, infile, outpath, name)
