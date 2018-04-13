from __future__ import print_function
from __future__ import division
import argparse
import collections
import os
import sys
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


def gethxb2(dict):
    """
    Fuction to find and return hxb2 seq and name from dict of sequences
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
            print("Found hxb2 ref. seq. Its full name is: %s" %(hxb2_key))
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")
        return None, None
    return str(hxb2_key), str(hxb2_seq)


def posnumcalc(hxb2seq, start):
    """
    gets positional hxb2 numbering
    :param hxb2seq: hxb2 sequence
    :param start: start position of alignment relative to hxb2 (if not at start of gene)
    :return: list of hxb2 numbering
    """
    pos_num = []
    n = start
    s = 0.01
    m = len(hxb2seq) - 1
    for i, resi in enumerate(hxb2seq):
        if i == 0 and resi == '-':
            print("Can't start numbering. HXB2 starts with a gap. Use a longer HXB2 sequence for the numbering")
        if i != m:
            if resi != '-' and hxb2seq[i+1] != '-':
                pos_num.append(n)
                n += 1
            elif resi != '-' and hxb2seq[i+1] == '-':
                g = n
                pos_num.append(g)
            elif resi == '-' and hxb2seq[i+1] == '-':
                g = n + s
                pos_num.append(g)
                s += 0.01
            elif resi == '-' and hxb2seq[i+1] != '-':
                g = n + s
                pos_num.append(g)
                s = 0.01
                n += 1
        else:
            if resi != '-':
                pos_num.append(n)
            elif resi == '-':
                g = n + s
                pos_num.append(g)
    return pos_num


def main(infile, outpath, posis, name, start):
    print(infile)

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    # construct outfile string
    nam = name + ".csv"
    outfile = os.path.join(outpath, nam)

    # initialise outfile with headers (overwrites pre-existing file)
    with open(outfile, 'w') as handle:
        handle.write('Time,Haplotype,Frequency\n')

    d = fasta_to_dct(infile)
    hxb2key, hxb2seq = gethxb2(d)
    try:
        del d[hxb2key]
    except:
        sys.exit("No HXB2 sequence in file, can't assign numbering, exiting")

    # get indexs for positions of interest
    numering = posnumcalc(hxb2seq, start)
    site_index_list = []
    for i, num in enumerate(numering):
        if num in posis:
            site_index_list.append(i)

    # convert dict format to key = time point, value = sequence
    time_d = collections.defaultdict(list)
    for names, seq in d.items():
        time = str(names.split("_")[2][:-3])
        time_d[time].append(seq)

    # dict to collect all the difference haplotypes
    master_hap_d = collections.defaultdict(int)

    # dict to collect everything: time point, haplotype and its frequency
    master_profile = collections.defaultdict(dict)

    for time, seqlist in time_d.items():
        total = len(seqlist)
        haps = collections.defaultdict(int)
        for seq in seqlist:
            sites = []
            for pos in site_index_list:
                sites.append(seq[pos])
            h = "".join(sites)
            haps[h] += 1
            master_hap_d[h] = 1
        # convert count to frequency and store vals in master dict
        for resis, count in haps.items():
            master_profile[time][resis] = round((count / total*100), 4)

    # add null values for time points where haplotype only emerges later
    for time, dct in master_profile.items():
        hap_list = list(dct.keys())
        full_hap_list = list(master_hap_d.keys())
        extra_h = list(set(full_hap_list).difference(hap_list))
        for item in extra_h:
            dct[item] = 0

    # write out the frequency .csv file
    print("writing file")
    with open(outfile, 'a') as handle:
        for time, dct in sorted(master_profile.items(), key=lambda x: x[0]):
            for hap, frq in dct.items():
                handle.write("{0},{1},{2}\n".format(time, hap, str(frq)))

    print("multi-site haplotype frequencies have been calculated")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate the frequency of several positions as a haplotype',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input fasta file, with all the time points in one file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-p', '--posis', type=int, nargs='+', default=argparse.SUPPRESS,
                        help='a list of positions you want to analyse (eg: 210 211 223 240)',
                        required=True)
    parser.add_argument('-n', '--name', default=1, type=str,
                        help='the name of your outfile (without the suffix (.csv)',
                        required=True)
    parser.add_argument('-s', '--start', default=1, type=int,
                        help='the HXB2 start position of your alignment', required=False)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    posis = args.posis
    name = args.name
    start = args.start

    main(infile, outpath, posis, name, start)
