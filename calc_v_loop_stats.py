import argparse
import collections
import os
from itertools import groupby
import regex


__author__ = 'colin anthony'


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
    :param dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence and name, as a strings
    """
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dict.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dict[k]
            print("Found hxb2 ref. seq. Its full name is: {0}".format(hxb2_key))
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")
    return str(hxb2_key), str(hxb2_seq)


def find_loop(hxb2, loopstart, loopend):
    """
    :param hxb2: protein sequence of HXB2 from an alignment to your sequences
    :param loopstart: region of interest (start)
    :param loopend: region of interest (end)
    :return: index where loop/region starts in alignment, index where loop/region ends
    """
    options = {
    "V1start" : "(L[-]*K[-]*C)",  #C[-]*V[-]*S[-]*
    "V1end" : "(C[-]*S[-]*F[-]*N[-]*I)",
    "V2start" : "(G[-]*E[-]*I[-]*K[-]*N[-]*C[-]*S)",
    "V2end" : "(C[-]*N[-]*T[-]*S[-]*V)",
    "C2start" : "(K[-]*L[-]*T[-]*S[-]*C[-]*N)",
    "C2end" : "(S[-]*V[-]*E[-]*I[-]*N)",
    "V3start" : "(S[-]*V[-]*E[-]*I[-]*N[-]*C)",
    "V3end" : "(C[-]*N[-]*I[-]*S[-]*R)",
    "C3start" : "(N[-]*M[-]*R[-]*Q[-]*A[-]*H[-]*C)",
    "C3end" : "(I[-]*I[-]*F[-]*K[-]*Q)", # should be "(Y[-]*C[-]*N[-]*S[-]*T[-]*Q)",
    "V4start" : "(G[-]*G[-]*E[-]*F[-]*F[-]*Y[-]*C)",
    "V4end" : "(C[-]*R[-]*I[-]*K[-]*Q)",
    "C4start" : "(T[-]*I[-]*T[-]*L[-]*P[-]*C[-]*R)",
    "C4end" : "(L[-]*T[-]*R[-]*D[-]*G[-]*G)",
    "V5start" : "(L[-]*T[-]*R[-]*D[-]*G[-]*G[-]*N)",
    "V5end" : "(R[-]*P[-]*G[-]*G[-]*G)",
    }

    a = regex.search(options[loopstart], hxb2, regex.BESTMATCH)
    matchstart = a.end()
    b = regex.search(options[loopend], hxb2, regex.BESTMATCH)
    matchend = b.start()

    return matchstart, matchend


def loopseq(v, start, end):
    """
    :param v: a protein sequence string
    :param start: the loop start index (from find_loop)
    :param end: the loop end index (from find_loop)
    :return: a degapped string of the sequence between the start and end index's
    """
    try:
        seq = v[start:end+1]
    except:
        seq = v[start:]

    s = (seq.replace("-", ""))
    return s


def looplen(seq):
    """
    :param seq: a protein sequence string
    :return: loop_length
    """
    l = len(seq)
    return int(l)


def loopcharge(seq):
    """
    :param seq: a protein sequence string
    :return: loop av charge
    """
    charge = 0
    c = {"R": 1, "K": 1, "H": 1, "E": -1, "D": -1}
    for item in seq:
        if item in c.keys():
            charge += c[item]
    return str(charge)


def loop_all_charge(seq):
    """
    :param seq: a protein sequence string
    :return: number of pos and neg charged residues in the loop
    """
    charge_pos = 0
    charge_neg = 0
    pos = {"R": 1, "K": 1, "H": 1}
    neg = {"E": -1, "D": -1}
    for item in seq:
        if item in pos.keys():
            charge_pos += pos[item]
        elif item in neg.keys():
            charge_neg += neg[item]

    return charge_pos, charge_neg


def loopglyc(seq): #NOTE: detects overlapping glycosylation sites, IE: "NNST", as two sites.
    """
    :param seq: a protein sequence string
    :return: number of glycan sites in the sequence
    """
    glyc = 0
    motif = "(N[-]*[^P-][-]*[TS])"
    c = regex.finditer(motif, seq, overlapped=True)
    for match in c:
        glyc += 1
    return str(glyc)


def main(infile, outpath, name, v1, v2, c2, v3, c3, v4, c4, v5):
    """
    :param infile: (str) inpath name to sequences aligned in coding space, with HXB2
    :param outpath: (str) path to the outfile
    :param name: (str) outfile name prefix
    :param v1 to v5: the v-loop(s) to compute stats for
    :param dna: (bool) is the alignment dna?
    :return: a csv file with all the loop summary values for all sequences
    """
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    print("file name:", infile)

    name = name + "_loop_stats.csv"
    outfile = os.path.join(outpath, name)

    loops = {
        "v1": ("V1start", "V1end"),
        "v2": ("V2start", "V2end"),
        "c2": ("C2start", "C2end"),
        "v3": ("V3start", "V3end"),
        "c3": ("C3start", "C3end"),
        "v4": ("V4start", "V4end"),
        "c4": ("C4start", "C4end"),
        "v5": ("V5start", "V5end"),
            }

    # read sequences to dict (translate if necessary and get HXB2 seq
    d = fasta_to_dct(infile)

    hxb2key, hxb2seq = gethxb2(d)
    del d[hxb2key]

    # find requested loops
    todo = [x.lower() for x in (v1, v2, c2, v3, c3, v4, c4, v5) if x is not None]

    # write csv headers
    outfile_string1 = []
    for i in range(len(todo)):
        outfile_string1.append(str(todo[i].upper()) + "_loop_length,")
        outfile_string1.append(str(todo[i].upper()) + "_loop_av_charge,")
        outfile_string1.append(str(todo[i].upper()) + "_total_pos_charge,")
        outfile_string1.append(str(todo[i].upper()) + "_total_neg_charge,")
        outfile_string1.append(str(todo[i].upper()) + "_glycan_sites,")
    outfile_string1.append("sequence_id")
    outfile_string1.append("\n")
    outfilename_all = "".join(outfile_string1)
    with open(outfile, "w") as handle:
        handle.write("Time," + outfilename_all)

    # get loop locations
    sites_loc = []
    for item in todo:
        print(item)
        start, end = find_loop(hxb2seq, loops[str(item)][0], loops[str(item)][1])
        reg = (item, start, end)
        sites_loc.append(reg)

    # calculate loop stats for each sequence
    for k, v in d.items():
        ids = k.upper().split("_")
        # get time point
        time = ids[2][:-3]
        with open(outfile, "a") as handle:
            handle.write(time + ",")
        for item in sites_loc:
            loop_seq = loopseq(v, item[1], item[2])
            loop_len = str(looplen(loop_seq))
            loop_charge = str(loopcharge(loop_seq))
            loop_charge_all_pos, loop_charge_all_neg = loop_all_charge(loop_seq)
            loop_glyc_no = str(loopglyc(loop_seq))

            # write stats to csv
            with open(outfile, "a") as handle:
                handle.write(f"{loop_len},{loop_charge},{loop_charge_all_pos},{loop_charge_all_neg},{loop_glyc_no},")
        with open(outfile, "a") as handle:
            handle.write(f"{str(k)}\n")

    print("done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate HIV-1 variable loop statistics',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--infile', type=str, required=True,
                        help='The input fasta file (all time points in one file')
    parser.add_argument('-o', '--outpath', type=str, required=True,
                        help='The path to where the output file should go')
    parser.add_argument('-n', '--name', type=str, required=True,
                        help='The name of the output file prefix: ie: "CAP177"')
    parser.add_argument('-v1', const="v1", action='store_const', required=False,
                        help='include v1 loop')
    parser.add_argument('-v2', const="v2", action='store_const', required=False,
                        help='include v2 loop')
    parser.add_argument('-c2', const="c2", action='store_const', required=False,
                        help='include c2 region')
    parser.add_argument('-v3', const="v3", action='store_const', required=False,
                        help='include v3 loop')
    parser.add_argument('-c3', const="c3", action='store_const', required=False,
                        help='include c3 region')
    parser.add_argument('-v4', const="v4", action='store_const', required=False,
                        help='include v4 loop')
    parser.add_argument('-c4', const="c4", action='store_const', required=False,
                        help='include c4 region')
    parser.add_argument('-v5', const="v5", action='store_const', required=False,
                        help='include v5 loop')

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    v1 = args.v1
    v2 = args.v2
    c2 = args.c2
    v3 = args.v3
    c3 = args.c3
    v4 = args.v4
    c4 = args.c4
    v5 = args.v5

    main(infile, outpath, name, v1, v2, c2, v3, c3, v4, c4, v5)
