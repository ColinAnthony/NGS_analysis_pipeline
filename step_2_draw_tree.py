#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from Bio import SeqIO

__author__ = 'Colin Anthony'


def fasta_to_dct(fn):
    '''
    :param fn: a fasta file
    :return: a dictionary
    '''
    dct = collections.OrderedDict()
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        dct[seq_record.description.replace(" ", "_")] = str(seq_record.seq).replace("~", "-").upper()
    return dct


def main(infile, outpath, name):
    print(infile)
    outfile = os.path.join(outpath, name)

    d = fasta_to_dct(infile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input fasta file, with all the time points in one file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the prefix name of your outfile', required=True)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name

    main(infile, outpath, name)
