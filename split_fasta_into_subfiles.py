#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from itertools import groupby


__author__ = 'colin.anthony'


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


def gethxb2(dictionary):
    """
    finds HXB2 in the dictionary keys and returns the full name and the sequence
    :param dictionary: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    """
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dictionary.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dictionary[k]
            print("Found hxb2 ref. seq. Its full name is: ", hxb2_key)
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")

    return hxb2_key, str(hxb2_seq)


def main(infile, outpath, field):

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    # import the data and remove HXB2 from alignment
    alignment_d = fasta_to_dct(infile)
    hxb2_key, hxb2_seq = gethxb2(alignment_d)
    print(hxb2_key)
    if hxb2_key is None:
        print("no hxb2 to remove")
    elif hxb2_key is not None:
        print("removing hxb2")
        del alignment_d[hxb2_key]

    # adjust the field value for python zero indexing
    adjust_field = field

    # sort time-point/sample sequence entries into dictionaries by the specified unique field
    split_d = collections.defaultdict(dict)
    for name, sequence in alignment_d.items():
        unique_field = name.split("_")[0:adjust_field]
        new_name = "_".join(unique_field)
        split_d[new_name][name] = sequence

    # write the grouped sequences to their separate outfiles
    for group, group_d in split_d.items():
        out_name = group + "_sep.fasta"
        outfile = os.path.join(outpath, out_name)
        if os.path.isfile(outfile):
            print("{0} file already exists, overwriting file".format(outfile))
            os.unlink(outfile)
        for seq_name, seq in group_d.items():
            with open(outfile, 'a') as handle:
                handle.write(">{0}\n{1}\n".format(seq_name, seq))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Splits a fasta file into multiple fasta files, '
                                                 'based on a unique field in the fasta headers)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=str, required=True,
                        help='The input fasta file')
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-f', '--field', type=int, default=4, required=False,
                        help="The field that differentiates your samples/time points."
                             "(ie: for the sequence name 'CAP177_2000_004wpi_V3C4_GGGACTCTAGTG_28_001_0.953,"
                             " 1 = CAP177, 4 = CAP177_2000_004wpi_V3C4, etc...")

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    field = args.field

    main(infile, outpath, field)
