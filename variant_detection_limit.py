#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import argparse
import collections
# import scipy.stats as ss
from glob import glob
from Bio import SeqIO
import os


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


def main(inpath, outpath, freq):

    infiles_search = os.path.join(inpath, "*sep.fasta")
    list_infiles = glob(infiles_search)
    outfile = os.path.join(outpath, "Variant_detection_probabilities.csv")
    print(outfile)
    if os.path.isfile(outfile):
        os.remove(outfile)

    with open(outfile, 'w') as handle:
        handle.write("name,sequencing_depth,prob_to_detect_{},variant_freq_with_95_percent_prob_to_detect\n"
                     .format(str(freq)))

    for infile in list_infiles:
        name = os.path.split(infile)[-1].replace("_sep.fasta", "")
        d = fasta_to_dct(infile)
        num_consensus_seqs = len(d)
        p = freq/100
        # k = range(1, num_consensus_seqs, 1)
        # binomial_prob = ss.binom.pmf(k, num_consensus_seqs, p)
        # a = round(sum(binomial_prob) * 100, 2)

        var_detection_limit = round((1 - ((1 - p) ** num_consensus_seqs))*100, 2)
        var_freq_with_95perc_prob = round((1 - (0.05 ** (1 / num_consensus_seqs))) * 100, 3)

        out_string = "{0},{1},{2},{3}\n".format(str(name), str(num_consensus_seqs),
                                                str(var_detection_limit), str(var_freq_with_95perc_prob))
        with open(outfile, 'a') as handle:
            handle.write(out_string)


    print("Variant detection limits have been calculated")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Returns the probability to detect a variant at a given frequency '
                                                 'and the frequency (as perc) of a variant that can be detected '
                                                 'with 95 % probability',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--inpath', type=str, required=True,
                        help='The path to your patient/time point files (should be "/path/to/5haplotype" '
                             'if you ran the NGS_processing_pipeline')
    parser.add_argument('-o', '--outpath', type=str, required=True,
                        help='The path to where the output file should go')
    parser.add_argument('-f', '--freq', type=float, default=1, required=True,
                        help='The percent frequency threshold of the variant to detect')

    args = parser.parse_args()
    inpath = args.inpath
    freq = args.freq
    outpath = args.outpath

    main(inpath, outpath, freq)
