#!/usr/bin/python2
from __future__ import print_function
from __future__ import division
import argparse
import numpy as np
import collections
from Bio import SeqIO
import os
import regex
import pandas as pd
import sys
__author__ = 'colin.anthony001@gmail.com'



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


def posnumcalc(hxb2seq, start):
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


def d_freq_lists(dna_list):
    n = len(dna_list[0])
    dist_dict = {'A': [0]*n, 'C': [0]*n, 'G': [0]*n, 'T': [0]*n, '-': [0]*n}
    total = 0
    for seq in dna_list:
        total += 1
        for index, dna in enumerate(seq):
            dist_dict[dna][index] += 1

    for base, freqlist in dist_dict.items():
        for i, cnt in enumerate(freqlist):
            frq = round((cnt/total*100), 4)
            freqlist[i] = frq
        dist_dict[base] = freqlist

    return dist_dict


def p_freq_lists(prot_list):
    n = len(prot_list[0])
    dist_dict = {'A': [0]*n, 'C': [0]*n, 'D': [0]*n, 'E': [0]*n, 'F': [0]*n, 'G': [0]*n, 'H': [0]*n, 'I': [0]*n,
                 'K': [0]*n, 'L': [0]*n, 'M': [0]*n, 'N': [0]*n, 'P': [0]*n, 'Q': [0]*n, 'R': [0]*n, 'S': [0]*n,
                 'T': [0]*n, 'V': [0]*n, 'W': [0]*n,  'Y': [0]*n, '-': [0]*n, 'X': [0]*n}
    total = 0
    for seq in prot_list:
        total += 1
        for index, aa in enumerate(seq):
            dist_dict[aa][index] += 1

    for resi, freqlist in dist_dict.items():
        for i, cnt in enumerate(freqlist):
            frq = round((cnt/total*100), 4)
            freqlist[i] = frq
        dist_dict[resi] = freqlist

    return dist_dict


def reorder_freq_dict(master_profile, pos_num):
    pos_d = collections.defaultdict(dict)
    for time, dct in master_profile.items():
        for indx, pos in enumerate(pos_num):
            new_freq_list = []
            for resi, old_freq_list in dct.items():
                new_freq_list.append((resi, old_freq_list[indx]))
            pos_d[pos][time] = new_freq_list

    return pos_d


def jsd(t1, t2, pos, time):  # Jensen-shannon divergence
    '''
    :param t1: list of AA/DNA frequencies for a given position and time point
    :param t2: list of AA/DNA frequencies for a given position and (time point +1)
    :param pos: The number of the sequence position
    :param time: the value for the time point (004)
    :return:
    '''
    # adapted from @author: jonathanfriedman
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    x = [i[1]/100 for i in t1]
    y = [i[1]/100 for i in t2]
    x = np.array(x)
    y = np.array(y)
    d1 = x * np.log2(2 * x / (x + y))
    d2 = y * np.log2(2 * y / (x + y))
    d1[np.isnan(d1)] = 0
    d2[np.isnan(d2)] = 0
    d = float(0.5 * np.sum(d1 + d2))
    d = round(d, 6)

    return (str(time), str(pos), str(d))


def shannon(t, pos, time):
    '''
    :param t: list of AA/DNA frequencies for a given position and time point
    :param pos: The number of the sequence position
    :param time: the value for the time point (004)
    :return:
    '''
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    entropy = 0.0
    x = [i[1] / 100 for i in t]

    for aa_freq in x:
        if aa_freq == 0:
            continue
        else:
            entropy += aa_freq * np.log2(aa_freq)
    entropy = entropy * -1
    entropy = round(entropy, 6)

    return (str(time), str(pos), str(entropy))


def glyc_finder(d, pos_num, outfile):
    '''
    :param d: dictionary of aligned protein sequences
    :param pos_num: list of HXB2 number codes for the alignment
    :param outfile: string for the outfile
    :return: None, writes csv file of glycan frequencies over time
    '''

    regex_pattern = 'N[\-]*[^P\-][\-]*[TS][^P]'
    master_profile = collections.OrderedDict()
    sub_dict = collections.defaultdict(list)
    for name, seq in d.items():
        time = name.split("_")[2][:-3]
        sub_dict[time].append(seq)
    for time, seq_list in sub_dict.items():
        total = len(seq_list)
        s = collections.defaultdict(int)
        s_frq = collections.OrderedDict()

        for i in seq_list:
            l = regex.finditer(regex_pattern, i, regex.BESTMATCH, overlapped=True)
            for j in l:
                start = j.span()[0]
                s[pos_num[start]] += 1

        for site, count in s.items():
            s_frq[site] = (float(count) / total) * 100
        master_profile[time] = s_frq

    all_sites_list = []
    for time, pos_frq_d in master_profile.items():
        for pos, frq in pos_frq_d.items():
            if pos not in all_sites_list:
                all_sites_list.append(pos)
    all_sites_list.sort(key=int)
    print(all_sites_list)

    with open(outfile, 'w') as handle:
        handle.write("Time,")
        for i in range(len(all_sites_list)):
            handle.write(str("N{0} glycan,".format(all_sites_list[i])))
        handle.write("\n")
    s_d = collections.OrderedDict(sorted(master_profile.items()))
    for time, pos_frq_d in s_d.items():
        with open(outfile, 'a') as handle:
            handle.write(str(time) + ",")
            for position in all_sites_list:
                try:
                    f = str(pos_frq_d[position])
                except:
                    f = '0'
                handle.write(str(f) + ',')
            handle.write('\n')


def main(infile, outpath, dna, start, name):
    print(infile)

    aa_out_name =  name + "_aa_freq.csv"
    aa_outfile = os.path.join(outpath, "aa_freq", aa_out_name)
    se_out_name = name + "_shannon_entropy.csv"
    se_outfile = os.path.join(outpath, "entropy", se_out_name)
    jsd_outfile = name + "_jsd.csv"
    jsd_outfile = os.path.join(outpath, "entropy", jsd_outfile)
    glyc_out_name = name + "_glycan_freq.csv"
    glycan_outfile = os.path.join(outpath, "glycans", glyc_out_name)

    d = fasta_to_dct(infile)
    hxb2key, hxb2seq = gethxb2(d)
    del d[hxb2key]

    pos_num = posnumcalc(hxb2seq, start)

    time_d = collections.defaultdict(list)
    master_profile = collections.OrderedDict()
    time_list = []
    for names, seq in d.items():
        time = str(names.split("_")[2][:-3])
        time_list.append(time)
        time_d[time].append(seq)

    if dna is True:
        func = d_freq_lists
    else:
        func = p_freq_lists

    # calculate AA/DNA frequencies per time point
    for time, seqlist in time_d.items():
        master_profile[time] = func(seqlist)

    pos_d = reorder_freq_dict(master_profile, pos_num)
    ksort = sorted(pos_d.keys())

    if dna is False:
        # calc glyc sites freq over time
        glyc_finder(d, pos_num, glycan_outfile)

        # write out aa frequencies
        with open(aa_outfile, 'w') as handle:
            handle.write("HXB2_position,Time,-,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,X,Y\n")
        for pos, time_d in sorted(pos_d.items()):
            for time, f_list in sorted(time_d.items()):
                frq_out = []
                for resi, frq in sorted(f_list):
                    frq_out.append(str(frq))
                with open(aa_outfile, 'a') as handle:
                    handle.write(str(pos) + "," + str(time) + "," + (",".join(frq_out)) + "\n")

    with open(jsd_outfile, 'w') as handle:
        handle.write('Time,Position,Jensen-Shannon_divergence' + "\n")

    with open(se_outfile, 'w') as handle1:
        handle1.write('Time,Position,Shannon_entropy' + "\n")

    for pos in ksort:
        nksort = sorted(pos_d[pos])
        for indx, time in enumerate(nksort):
            if indx == 0:
                # print('This must be the first time point:', time, "\nposition:", pos)
                t1 = nksort[indx]
                rsort_tn1 = sorted(pos_d[pos][t1])
            else:
                t1 = nksort[indx - 1]
                rsort_tn1 = sorted(pos_d[pos][t1])

            t2 = nksort[indx]
            rsort_tn2 = sorted(pos_d[pos][t2])

            js_dist = jsd(rsort_tn1, rsort_tn2, pos, time)
            s_entropy = shannon(rsort_tn2, pos, time)

            outwrite = ",".join(js_dist)
            outwrite1 = ",".join(s_entropy)

            with open(jsd_outfile, 'a') as handle:
                handle.write(outwrite + "\n")

            with open(se_outfile, 'a') as handle1:
                handle1.write(outwrite1 + "\n")

    print("entropy calculations complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate the Jensen-Shannon entropy for a '
                                                 'longitudinal sequence alignment',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input fasta file, with all the time points in one file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The parent path containing the aa_freq, entropy and glycans subfolders. '
                             '(usually "/path/to/6analysis")', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the name of your outfile', required=True)
    parser.add_argument('-d', '--dna', default=False, action='store_true',
                        help='is the fasta file a DNA sequence? (default = False (protein))', required=False)
    parser.add_argument('-s', '--start', default=1, type=int,
                        help='the HXB2 start position of your alignment', required=False)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    dna = args.dna
    start = args.start
    name = args.name

    main(infile, outpath, dna, start, name)
