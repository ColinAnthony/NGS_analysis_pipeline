#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
from glob import glob
import argparse
import collections
from Bio import SeqIO
import subprocess


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


def gethxb2(dict):
    '''
    :param dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    '''
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dict.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dict[k]
            print("Found hxb2 ref. seq. Its full name is: %s" %(hxb2_key))
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")
    return str(hxb2_key), str(hxb2_seq)


def main(infile, name, limit, root, script_folder):
    print(infile)



    # glob the 5haplotype folder for all hap.fasta files
    path_aln, name_aln = os.path.split(infile)[0]
    parent_path = os.path.split(path_aln)[0]
    tree_haplotype_folder = os.path.join(parent_path, "5haplotype", "tree_haplotype")
    os.makedirs(tree_haplotype_folder)
    top_hap_outfile = os.path.join(tree_haplotype_folder, "{0}_top_{1}.fasta".format(name, limit))
    tree_analysis_folder = os.path.join(parent_path, "6analysis", "tree")


    # split the aligned file into subfiles
    aln_d = fasta_to_dct(infile)
    hxb2_key, hxb2_seq = gethxb2(aln_d)
    tmp_file = "No_hxb2_align.fasta"
    with open(tmp_file, 'w') as handle:
        for name, seq in aln_d.items():
            handle.write(">{0}\n{1}\n".format(name, seq))

    split_by_unique = os.path.join(script_folder, 'split_fasta_into_subfiles.py')
    cmd1 = 'python3 {0} -i {1} -o {2}'.format(split_by_unique, tmp_file, tree_haplotype_folder)

    subprocess.call(cmd1, shell=True)

    # run haplotype script on subfiles
    split_fasta_files = os.path.join(tree_haplotype_folder, "*_sep.fasta")
    haplotyper = os.path.join(script_folder, 'calc_haplotype_freq.py')
    for split_fasta in glob(split_fasta_files):
        cmd2 = 'python3 {0} -i {1} -o {2}'.format(haplotyper, split_fasta, tree_haplotype_folder)

        subprocess.call(cmd2, shell=True)

    # cat the top n haplotyped files into one file
    search_haplotype = os.path.join(tree_haplotype_folder, "*hap.fasta")
    if root is not None:
        rank = 1
        root = SeqIO.parse(root, 'fasta')
        with open(top_hap_outfile, 'w') as handle:
            handle.write(">{0}\n{1}\n".format(root.name, str(root.seq)))
    else:
        rank = 0

    for hap in glob(search_haplotype):
        top_rec_counter = 0
        records = SeqIO.parse(hap, 'fasta')
        if limit is None:
            limit = len(records)

        for record in records:
            # get top record from first time point
            if root is None and rank == 0 and top_rec_counter == 0:
                root_name = record.name

            while top_rec_counter < limit:
                name = record.name
                if not "HXB2" in name.upper():
                    seq = str(record.seq).upper()
                    top_rec_counter += 1
                    with open(top_hap_outfile, 'a') as handle:
                        handle.write(">{0}\n{1}\n".format(name, seq))

        rank += 1

    # build the tree
    tree_infile = top_hap_outfile
    tree_outfile = tree_infile.replace(".fasta", "nwk")
    cmd = "fasttree -nt -gtr {0} > {1}".format(tree_infile, tree_outfile)
    subprocess.call(cmd, shell=True)
    tree_fig_name = name + "_tree"
    tree_script = os.path.join(script_folder, "bubble_tree.py")
    scale = 50000
    cmd3 = "python3 {0} -i {1} -o {2} -n {3} -r {4} -z -v {5}".format(tree_script, tree_outfile,
                                                                      tree_analysis_folder, tree_fig_name, root, scale)

    subprocess.call(cmd3, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='take a processed alignment file and produce a tree',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str, required=False,
                        help='The input fasta file, with all the time points in one file')
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str, required=False,
                        help='the prefix name of your outfile')
    parser.add_argument('-l', '--limit', default=None, type=str, required=False,
                        help='the number of sequences to select from each file')
    parser.add_argument('-r', '--root', default=None, type=str, required=False,
                        help='fasta file containing the sequence to use if an external root sequence is required')
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str, required=True,
                        help='the path to the folder containing the pipeline scripts')

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    limit = args.limit
    root = args.root
    script_folder = args.script_folder

    main(infile, name, limit, root, script_folder)
