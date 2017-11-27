#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
from itertools import groupby
from glob import glob
import argparse
import collections
import subprocess


__author__ = 'Colin Anthony'


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

    return str(hxb2_key), str(hxb2_seq)


def main(infile, outpath, name, limit, root, script_folder, sample_type):

    # get out_folder
    full_file_path = os.path.abspath(infile)
    path = os.path.split(full_file_path)[0]
    parent_path = os.path.split(path)[0]

    # make the output folder
    if outpath is None:
        tree_haplotype_folder = os.path.join(parent_path, "5haplotype", "tree_haplotype")
        full_outpath = tree_haplotype_folder
    else:
        full_outpath = os.path.abspath(outpath)
        tree_haplotype_folder = os.path.join(full_outpath, "tree_haplotype")

    if not os.path.exists(tree_haplotype_folder):
        print("out_folder does not exist, creating it")
        try:
            os.makedirs(tree_haplotype_folder, exist_ok=True)
        except PermissionError as e:
            print(e)

    # set the tree output folder
    if outpath is None:
        tree_analysis_folder = os.path.join(parent_path, "6analysis", "tree")
    else:
        tree_analysis_folder = os.path.join(outpath, "6analysis", "tree")
    if not os.path.exists(tree_analysis_folder):
        print("out_folder does not exist, creating it")
        try:
            os.makedirs(tree_analysis_folder, exist_ok=True)
        except PermissionError as e:
            print(e)
            sys.exit()
    try:
        os.makedirs(tree_haplotype_folder, exist_ok=True)
    except PermissionError as e:
        print(e)
        sys.exit()

    # split the aligned file into subfiles
    split_by_unique = os.path.join(script_folder, 'split_fasta_into_subfiles.py')
    cmd1 = 'python3 {0} -i {1} -o {2}'.format(split_by_unique, infile, tree_haplotype_folder)
    subprocess.call(cmd1, shell=True)

    # run haplotype script on subfiles
    split_fasta_files = os.path.join(tree_haplotype_folder, "*sep.fasta")
    haplotyper = os.path.join(script_folder, 'calc_haplotype_freq.py')
    for split_fasta in glob(split_fasta_files):
        cmd2 = 'python3 {0} -i {1} -o {2}'.format(haplotyper, split_fasta, tree_haplotype_folder)
        subprocess.call(cmd2, shell=True)

    # cat the top n haplotyped files into one file
    out_name = name + "_top_{}hap.fasta".format(str(limit))
    top_hap_outfile = os.path.join(tree_analysis_folder, out_name)
    if os.path.isfile(top_hap_outfile):
        print("{0} file already exists, overwriting file".format(top_hap_outfile))
        os.unlink(top_hap_outfile)
    search_haplotype = os.path.join(tree_haplotype_folder, "*hap.fasta")

    # if an external root is specified, get the name and sequence
    if root is not None:
        root_d = fasta_to_dct(root)
        if len(root_d) != 1:
            print("Incorrect number of sequences in {0} file, only one root allowed".format(root))
        for rname, rseq in root_d.items():
            root_name = rname
            root_seq = rseq
            with open(top_hap_outfile, 'a') as handle:
                handle.write(">{0}\n{1}\n".format(root_name, root_seq))

    rank = 0
    root_name = root
    for hap_file in glob(search_haplotype):
        top_rec_counter = 0
        records = fasta_to_dct(hap_file)

        for seq_name, seq in records.items():
            # get top record from first time point
            if root is None and rank == 0:
                root_name = seq_name
                rank += 1
                print("root", root_name)
            if top_rec_counter < limit:
                seq = seq.upper()
                print('name', seq_name)
                top_rec_counter += 1
                print(rank)
                with open(top_hap_outfile, 'a') as handle:
                    handle.write(">{0}\n{1}\n".format(seq_name, seq))



    # build the tree
    tree_infile = top_hap_outfile
    tree_outfile = tree_infile.replace(".fasta", ".nwk")
    if not sample_type:
        flag = '-nt'
    else:
        flag = ''

    cmd = "fasttree -gtr {0} {1} > {2}".format(flag, tree_infile, tree_outfile)
    subprocess.call(cmd, shell=True)

    # plot the tree figure
    tree_fig_name = name + "_tree"
    tree_script = os.path.join(script_folder, "bubble_tree.py")
    scale = 50000

    cmd3 = "python3 {0} -i {1} -o {2} -n {3} -r {4} -z -v {5}".format(tree_script, tree_outfile, tree_analysis_folder,
                                                                      tree_fig_name, root_name, scale)
    print("plot_script: ", cmd3)
    subprocess.call(cmd3, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Take a cleaned alignment fasta file, splits the file into subfiles '
                                                 'based on sample names/time-points'
                                                 'calculates the haplotype frequencies for each '
                                                 '(collapses by identical sequences)'
                                                 ' and produce a tree and image of the tree',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str, required=True,
                        help='The input fasta file, with all the time points in one file')
    parser.add_argument('-o', '--outpath', default=None, type=str, required=False,
                        help='the path to where the intermediate files will be created (the haplotyped files)'
                             'Only use this if you did not run the full pipeline or want a custom out_folder.'
                             'By default the out_folder "5haplotype" will be used')
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str, required=True,
                        help='the prefix name of your outfile')
    parser.add_argument('-l', '--limit', default=10, type=int, required=False,
                        help='the number of sequences to select from each file')
    parser.add_argument('-r', '--root', default=None, type=str, required=False,
                        help='fasta file containing the sequence to use if an external root sequence is required')
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str, required=True,
                        help='the path to the folder containing the pipeline scripts')
    parser.add_argument('-st', '--sample_type', default=False, action='store_true', required=False,
                        help='Use this flag if the file contains protein sequences. Default = DNA')
    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    limit = args.limit
    root = args.root
    script_folder = args.script_folder
    sample_type = args.sample_type
    main(infile, outpath, name, limit, root, script_folder, sample_type)
