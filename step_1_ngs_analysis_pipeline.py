#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import sys
import os
import argparse
import subprocess
from glob import glob


__author__ = 'Colin Anthony'


def aa_plotter(script_file, infile, outpath, ssnab, bnab):

    print("plotting amino acid frequencies")
    if bnab is None and ssnab is not None:
        ab_markers = '-t {0}'.format(ssnab)
    elif bnab is not None and ssnab is None:
        ab_markers = '-b {0}'.format(bnab)
    elif bnab is not None and ssnab is not None:
        ab_markers = '-t {0} -b {1}'.format(ssnab, bnab)
    else:
        ab_markers = ''

    cmd1 = "python3 {0} -i {1} -o {2} {3}".format(script_file, infile, outpath, ab_markers)

    subprocess.call(cmd1, shell=True)


def glycan_plotter(script_file, infile, outpath, ssnab, bnab):

    print("plotting glycan frequencies")
    if bnab is None and ssnab is not None:
        ab_markers = '-t {0}'.format(ssnab)
    elif bnab is not None and ssnab is None:
        ab_markers = '-b {0}'.format(bnab)
    elif bnab is not None and ssnab is not None:
        ab_markers = '-t {0} -b {1}'.format(ssnab, bnab)
    else:
        ab_markers = ''

    cmd1 = "python3 {0} -i {1} -o {2} {3}".format(script_file, infile, outpath, ab_markers)

    subprocess.call(cmd1, shell=True)


def entropy_plotter(script_file, infiles, outpath, ssnab, bnab):

    print("plotting entropy heatmaps")
    if bnab is None and ssnab is not None:
        ab_markers = '-t {0}'.format(ssnab)
    elif bnab is not None and ssnab is None:
        ab_markers = '-b {0}'.format(bnab)
    elif bnab is not None and ssnab is not None:
        ab_markers = '-t {0} -b {1}'.format(ssnab, bnab)
    else:
        ab_markers = ''

    for infile in infiles:
        name = os.path.split(infile)[-1].replace(".csv", "")

        cmd1 = "python3 {0} -i {1} -o {2} -n {3} {4}".format(script_file, infile, outpath, name, ab_markers)

        subprocess.call(cmd1, shell=True)


def divergence_plotter(script_file, infiles, outpath, ssnab, bnab, vl_file):
    print("plotting divergence")

    if vl_file is None:
        vl_inflie = ''
    else:
        vl_inflie = '-i2 {0}'.format(vl_file)

    if bnab is None and ssnab is not None:
        ab_markers = '-t {0}'.format(ssnab)
    elif bnab is not None and ssnab is None:
        ab_markers = '-b {0}'.format(bnab)
    elif bnab is not None and ssnab is not None:
        ab_markers = '-t {0} -b {1}'.format(ssnab, bnab)
    else:
        ab_markers = ''

    for infile in infiles:
        name = os.path.split(infile)[-1].replace("_divergence.csv", "")
        cmd1 = "python3 {0} -i1 {1} {2} -o {3} -n {4} {5}".format(script_file, infile, vl_inflie,
                                                                  outpath, name, ab_markers)

        subprocess.call(cmd1, shell=True)


def loop_stats_plotter(script_file, infiles, outpath, ssnab, bnab, vl_file):
    print("plotting loop_stats")
    if vl_file is None:
        vl_inflie = ''
    else:
        vl_inflie = '-i2 {0}'.format(vl_file)
    if bnab is None and ssnab is not None:
        ab_markers = '-t {0}'.format(ssnab)
    elif bnab is not None and ssnab is None:
        ab_markers = '-b {0}'.format(bnab)
    elif bnab is not None and ssnab is not None:
        ab_markers = '-t {0} -b {1}'.format(ssnab, bnab)
    else:
        ab_markers = ''

    for infile in infiles:
        name = os.path.split(infile)[-1].replace("_loop_stats.csv", "")
        cmd1 = "python3 {0} -i1 {1} {2} -o {3} -n {4} {5}".format(script_file, infile, vl_inflie,
                                                                  outpath, name, ab_markers)

        subprocess.call(cmd1, shell=True)


def main(alignment, viral_load_file, parent_folder, start, script_folder, freq, ab_time, bnab_time,
         reference, longitudinal, env, loops, run_step):

    print(alignment)
    name = "_".join(os.path.split(alignment)[-1].split("_")[:2])

    analysis_folder = os.path.join(parent_folder, "6analysis")
    haplotypes = os.path.join(parent_folder, "5haplotypes")
    analysis_sub_folders = os.listdir(analysis_folder)
    print(analysis_folder)
    # get viral load csv

    check_folders = ['aa_freq', 'entropy', 'glycans', "divergence", "loops"]
    for folder in check_folders:
        if folder not in analysis_sub_folders:
            print("could not find analysis subfolder {0}, "
                  "check that you ran part 1 of the pipeline correctly".format(folder))
            sys.exit()

    # step 1 do the calculations
    if run_step == 1:
        run_step += 1
        # calc detection level
        detection_limit_script = os.path.join(script_folder, "calc_variant_detection_limit.py")
        cmd4 = "python3 {0} -i {1} -o {2} -f {3}".format(detection_limit_script, haplotypes, parent_folder, freq)

        subprocess.call(cmd4, shell=True)

        # calc aa freq, glycs, entropy
        print("calculating aa freq, glycan freq and entropy")
        entropy_aa_calc_script = os.path.join(script_folder, "calc_entropy_aa_glycan_freq.py")
        cmd1 = "python3 {0} -i {1} -o {2} -n {3} -s {4}".format(entropy_aa_calc_script, alignment, analysis_folder,
                                                                name, start)
        subprocess.call(cmd1, shell=True)

        # calc divergence
        if longitudinal:
            print("calculating divergence")
            divergence_script = os.path.join(script_folder, "calc_divergence.py")
            divergence_folder = os.path.join(analysis_folder, "divergence")
            cmd2 = "python3 {0} -r {1} -i {2} -o {3} -n {4}".format(divergence_script, reference,
                                                                    alignment, divergence_folder, name)

            subprocess.call(cmd2, shell=True)

        # calc loop stats
        if env:
            print("calculating v-loop statistics")
            loop_stats_script = os.path.join(script_folder, "calc_v_loop_stats.py")
            loops_folder = os.path.join(analysis_folder, "loops")
            loop_arg = " ".join(["-" + str(x) for x in loops])

            cmd3 = "python3 {0} -i {1} -o {2} -n {3} {4}".format(loop_stats_script, alignment, loops_folder, name,
                                                                 loop_arg)

            subprocess.call(cmd3, shell=True)

    # graph aa_freq
    if run_step == 2:
        run_step += 1

        aa_frq_script_file = os.path.join(script_folder, "plot_aa_freq.py")
        aa_frq_search = os.path.join(analysis_folder, "aa_freq", "*aa_freq.csv")
        aa_frq_inpath = glob(aa_frq_search)
        aa_freq_infile = aa_frq_inpath[0]
        aa_frq_outpath = os.path.join(analysis_folder, "aa_freq")

        aa_plotter(aa_frq_script_file, aa_freq_infile, aa_frq_outpath, ab_time, bnab_time)

    # plot glycan freqs
    if run_step == 3:
        run_step += 1

        glyc_script_file = os.path.join(script_folder, "plot_glycan_freq.py")
        glyc_search = os.path.join(analysis_folder, "glycans", "*glycan_freq.csv")
        glyc_inpath = glob(glyc_search)
        glyc_infile = glyc_inpath[0]
        glyc_outpath = os.path.join(analysis_folder, "glycans")

        glycan_plotter(glyc_script_file, glyc_infile, glyc_outpath, ab_time, bnab_time)

    # plot entropy
    if run_step == 4:
        run_step += 1
        entropy_script_file = os.path.join(script_folder, "plot_entropy.py")
        entropy_search = os.path.join(analysis_folder, "entropy", "*.csv")
        entropy_inpath = glob(entropy_search)
        entropy_outpath = os.path.join(analysis_folder, "entropy")

        entropy_plotter(entropy_script_file, entropy_inpath, entropy_outpath, ab_time, bnab_time)

    # plot divergence
    if run_step == 5:
        run_step += 1
        if longitudinal:
            divergence_script_file = os.path.join(script_folder, "plot_divergence.py")
            divergence_search = os.path.join(analysis_folder, "divergence", "*divergence.csv")
            divergence_inpath = glob(divergence_search)
            divergence_outpath = os.path.join(analysis_folder, "divergence")

            divergence_plotter(divergence_script_file, divergence_inpath, divergence_outpath, ab_time, bnab_time,
                               viral_load_file)

    # plot loop stats
    if run_step == 6:
        run_step += 1
        if env:
            loop_stats_script_file = os.path.join(script_folder, "plot_loop_stats.py")
            loop_stats_search = os.path.join(analysis_folder, "loops", "*loop_stats.csv")
            loop_stats_inpath = glob(loop_stats_search)
            loop_stats_outpath = os.path.join(analysis_folder, "loops")

            loop_stats_plotter(loop_stats_script_file, loop_stats_inpath, loop_stats_outpath, ab_time, bnab_time,
                               viral_load_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NGS data analysis pipeline',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i1', '--infile1', default=argparse.SUPPRESS, type=str, required=True,
                        help='The name of the aligned protein fasta file, with all the time points/samples in one file')
    parser.add_argument('-i2', '--viral_load_file', default=None, type=str, required=False,
                        help='The csv file with column for: patient name, time points and viral load')
    parser.add_argument('-p', '--path', default=argparse.SUPPRESS, type=str, required=True,
                        help='The path to the gene region folder, created by running part1_ngs_processing_pipeline')
    parser.add_argument('-s', '--start', default=1, type=int, required=False,
                        help='the HXB2 start position of your alignment')
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str, required=True,
                        help='the path to the folder containing the pipeline scripts')
    parser.add_argument('-f', '--freq', type=float, default=1, required=False,
                        help='The percent frequency threshold of the variant to detect')
    parser.add_argument('-t', '--ab_time', type=int, required=False,
                        help='The time to mask, ie: start of nAb "-t 19"')
    parser.add_argument('-b', '--bnab_time', type=int, required=False,
                        help='The second time point to mask, ie: start of bnnAb "-b 50"')
    parser.add_argument('-r', '--reference', type=str, required=False,
                        help='The fasta file with only the reference sequence. '
                             'Usually taken from the most abundant haplotype at the first time point')
    parser.add_argument('-l', '--longitudinal', default=False, action="store_true", required=False,
                        help='If you have longitudinal sequences and want to calculate divergence, use this flag')
    parser.add_argument('-e', '--env', default=False, action="store_true", required=False,
                        help='If you have envelope sequences and want to characterise the v-loop, use this flag')
    parser.add_argument('-vl', '--loops', nargs="+", type=str, required=False,
                        help='If you set -e, which v-loops are in your sequences? (eg: -vl c3 v4 c4 v5')
    parser.add_argument('-rs', '--run_step', type=int, default=1,
                        help='rerun the pipeline from a given step:\n'
                             '1 = step 1: do the analysis calculations\n'
                             '2 = step 2: plot aa frequencies\n'
                             '3 = step 3: plot glycan frequencies\n'
                             '4 = step 4: plot entropy\n'
                             '5 = step 5: plot divergence (optional)\n'
                             '6 = step 6: plot loop statistics', required=False)

    args = parser.parse_args()
    infile1 = args.infile1
    viral_load_file = args.viral_load_file
    start = args.start
    path = args.path
    script_folder = args.script_folder
    freq = args.freq
    ab_time = args.ab_time
    bnab_time = args.bnab_time
    reference = args.reference
    longitudinal = args.longitudinal
    env = args.env
    loops = args.loops
    run_step = args.run_step

    main(infile1, viral_load_file, path, start, script_folder, freq, ab_time, bnab_time, reference, longitudinal, env,
         loops, run_step)
