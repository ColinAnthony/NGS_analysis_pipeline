#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import sys
import os
import argparse
import collections
from Bio import SeqIO
import subprocess


__author__ = 'Colin Anthony'


def main(alignment, parent_folder, start, script_folder, freq, reference, longitudinal, env, loops):

    print(alignment)
    name = "_".join(os.path.split(alignment)[-1].split("_")[:2])

    analysis_folder = os.path.join(parent_folder, "6analysis")
    haplotypes = os.path.join(parent_folder, "5haplotypes")
    analysis_sub_folders = os.listdir(analysis_folder)
    print(analysis_sub_folders)
    # get viral load csv

    check_folders = ['aa_freq', 'entropy', 'glycans', "divergence", "loops"]
    for folder in check_folders:
        if folder not in analysis_sub_folders:
            print("could not find analysis subfolder {0}, "
                  "check that you ran part 1 of the pipeline correctly".format(folder))
            sys.exit()

    # calc aa freq, glycs, entropy
    print("calculating aa freq, glycan freq and entropy")
    entropy_aa_calc_script = os.path.join(script_folder, "entropy_aa_freq_calculator.py")
    cmd1 = "python3 {0} -i {1} -o {2} -n {3} -s {4}".format(entropy_aa_calc_script, alignment, analysis_folder,
                                                            name, start)
    subprocess.call(cmd1, shell=True)

    # graph aa_freq, entropy, glycans

    # calc divergence
    if longitudinal:
        print("calculating divergence")
        divergence_script = os.path.join(script_folder, "divergence_calculator.py")
        divergence_folder = os.path.join(analysis_folder, "divergence")
        cmd2 = "python3 {0} -r {1} -i {2} -o {3} -n {4}".format(divergence_script, reference,
                                                                alignment, divergence_folder, name)

        subprocess.call(cmd2, shell=True)

        # graph divergence

    # calc loop stats
    if env:
        print("calculating v-loop statistics")
        loop_stats_script = os.path.join(script_folder, "characterise_v_loops.py")
        loops_folder = os.path.join(analysis_folder, "loops")
        loop_arg = " ".join(["-" + str(x) for x in loops])

        cmd3 = "python3 {0} -i {1} -o {2} -n {3} {4}".format(loop_stats_script, alignment, loops_folder, name, loop_arg)

        subprocess.call(cmd3, shell=True)

        # graph loopstats

    # calc detection level
    detection_limit_script = os.path.join(script_folder, "variant_detection_limit.py")

    cmd4 = "python3 {0} -i {1} -o {2} -f {3}".format(detection_limit_script, haplotypes, parent_folder, freq)

    subprocess.call(cmd4, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NGS data analysis pipeline',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str, required=True,
                        help='The name of the aligned protein fasta file, '
                             'with all the time points/samples in one file')
    parser.add_argument('-p', '--path', default=argparse.SUPPRESS, type=str, required=True,
                        help='The path to the gene region folder, created by running part1_ngs_processing_pipeline')
    parser.add_argument('-s', '--start', default=1, type=int, required=False,
                        help='the HXB2 start position of your alignment')
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str, required=True,
                        help='the path to the folder containing the pipeline scripts')
    parser.add_argument('-f', '--freq', type=float, default=1, required=False,
                        help='The percent frequency threshold of the variant to detect')
    parser.add_argument('-r', '--reference', type=str, required=False,
                        help='The fasta file with only the reference sequence. '
                             'Usually taken from the most abundant haplotype at the first time point')
    parser.add_argument('-l', '--longitudinal', default=False, action="store_true", required=False,
                        help='If you have longitudinal sequences and want to calculate divergence, use this flag')
    parser.add_argument('-e', '--env', default=False, action="store_true", required=False,
                        help='If you have envelope sequences and want to characterise the v-loop, use this flag')
    parser.add_argument('-vl', '--loops', nargs="+", type=str, required=False,
                        help='If you set -e, which v-loops are in your sequences? (eg: -vl c3 v4 c4 v5')

    args = parser.parse_args()
    infile = args.infile
    start = args.start
    path = args.path
    script_folder = args.script_folder
    freq = args.freq
    reference = args.reference
    longitudinal = args.longitudinal
    env = args.env
    loops = args.loops

    main(infile, path, start, script_folder, freq, reference, longitudinal, env, loops)
