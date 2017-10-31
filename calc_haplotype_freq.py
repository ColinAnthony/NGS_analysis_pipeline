#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import collections
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


def hapl_collapse(fastafile, seq_name, outfile):
    """
    Haplotype a fasta file by collapsing identical sequences into a single sequence with a frequency count.
    :param fastafile: input fasta file
    :param seq_name: (str) prefix for sequence name
    :param outfile: the path to where the outfile is created
    :return: dictionary of key = sequence, value  = count of that sequence
    """
    # import seqs to dict and get count of each unique sequence
    total = 0
    dct = collections.defaultdict(int)
    for seq_name, seq in py3_fasta_iter(fastafile):
        dct[str(seq).upper()] += 1
        total += 1

    # list of sequences, sorted by count, max to min
    sorted_sequences_by_count = sorted(dct, key=dct.get, reverse=True)

    # clear any existing file previously created
    if os.path.isfile(outfile):
        print("{0} already exists, overwriting file".format(outfile))
        os.unlink(outfile)

    # initialize a counter to add an incremental number to each sequence
    for counter, sequence in enumerate(sorted_sequences_by_count):
        frq = round(dct[sequence] / float(total), 3)
        num = str(counter).zfill(3)
        n = seq_name + "_" + str(num) + "_" + str(frq)

        # write the haplotyped sequence to file
        with open(outfile, "a") as handle:
            handle.write('>{0}\n{1}\n'.format(n, sequence))


def main(infile, outpath, field):

    print("Collapsing to haplotypes on file: \n\t {}".format(infile))

    # set the outfile name
    full_inpath = os.path.abspath(infile)
    full_outpath = os.path.abspath(outpath)
    name = os.path.split(full_inpath)[-1]
    name = name.replace("_sep.fasta", "_hap.fasta")
    out = os.path.join(full_outpath, name)

    # adjust field for python zero indexing
    adjust_field = field + 1

    # get prefix for sequence names from infile name
    seq_name = "_".join(name.split("_")[:adjust_field])

    # run the haplotyping function
    hapl_collapse(infile, seq_name, out)

    print("Done haplotyping")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='collapses sequences by unique, with freq added to the name',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=str, required=True,
                        help='The input fata file')
    parser.add_argument('-o', '--outpath', required=True, type=str,
                        help='path for output.')
    parser.add_argument('-f', '--field', type=int, default=4, required=False,
                        help="The field that differentiates your samples/time points."
                             "(ie: for the sequence name 'CAP177_2000_004wpi_V3C4_GGGACTCTAGTG_28_001_0.953,"
                             " 1 = CAP177, 4 = CAP177_2000_004wpi_V3C4, etc...")

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    field = args.field
    main(infile, outpath, field)
