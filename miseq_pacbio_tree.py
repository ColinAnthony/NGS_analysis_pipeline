from __future__ import print_function
from __future__ import division
import argparse
import sys, os
import collections
from itertools import groupby

try:
    from ete3 import Tree, faces, TreeStyle, NodeStyle, TreeNode, add_face_to_node, SequenceFace, SeqMotifFace
except ImportError:
    print("Ete3 is not installed correctly. For best results, install ete3 from the anaconda distribution, "
          "as per http://etetoolkit.org/download/")


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


def highlighter_dct(d, consens, colours):
    """
    :param d: a dictionary of names and sequences
    :return: a dictionary with sequences modified to show difference from a given template
    """
    ref = d[consens]
    del d[consens]
    if colours == 1:
        for name, seq in d.items():
            s = ''
            for idx, pos in enumerate(seq):
                if pos != '-':
                    s += 'X'
                else:
                    s += str('-')
            d[name] = s

    if colours == 2:
        for name, seq in d.items():
            s = ''
            for idx, pos in enumerate(seq):
                if pos == ref[idx]:
                    s += 'X'
                else:
                    s += str(pos)
            d[name] = s

    d[consens] = ref

    return d


def get_time(tree, field2):
    """
    :param tree: tree object parsed from ete module
    :param field2: the '_' delimited field containing the sample/time point identifier
    :return: non_redundant list of time points, number of time points
    """
    times = []
    for node in tree.traverse():
        if node.is_leaf() == True:
            try:
                name = node.name.split("_")
                t = name[field2]
            except:
                t = 'zero'
                print("incorrect name format: sequence assigned time point of 'zero")

            if t not in times:
                times.append(t)

    times.sort()

    return times


def col_map(times):
    """
    :param times: non_redundant list of time points from get_times function
    :param n_color: number of time points
    :return: dictionary mapping colour spectrum to time points, key = time point, value = colour code
    """

    colour_25 = ['#E50001', '#E32A00', '#E15700', '#E08300', '#DEAE00', '#DDD800', '#B4DB00', '#88D900',
                 '#5CD800', '#31D600', '#06D500', '#00D323', '#00D24C', '#00D075', '#00CE9D', '#00CDC4',
                 '#00ABCB', '#0082CA', '#0059C8', '#0031C6', '#000AC5', '#1C00C3', '#4200C2', '#6800C0', '#8D00BF']
    # colour = ['#9e0142', '#d53e4f', '#f46d43', '#abdda4', '#fdae61', '#66c2a5', '#3288bd', '#5e4fa2']
    # colour_25 = [1, 2, 5, 6, 4, 4, 1, 1]
    other = ['#1a1a1a', '#E50001']
    # colour_25 = [colour[0], colour[1], colour[1], colour[2], colour[2], colour[2], colour[2], colour[2],
    #             colour[3], colour[3], colour[3], colour[3], colour[3], colour[3], colour[4], colour[4],
    #             colour[4], colour[4], colour[5], colour[5], colour[5], colour[5], colour[6], colour[7]]

    n_colors = len(times)
    l = len(colour_25)
    assert n_colors <= l, "The colour gradient in this script can only handle 25 time points, add more colours to the colour_25 list"

    if 'zero' in times:
        n_colors = n_colors - 1

    step = l // n_colors

    j = 0
    cd = collections.OrderedDict()
    for i, v in enumerate(times):
        if str(v) == 'zero':
            cd[str(v)] = other[0]
        elif n_colors == 1:
            cd[str(v)] = other[0]
        elif n_colors == 2:
            cd[str(v)] = other[i]
        else:
            cd[str(v)] = colour_25[j]
            j += step

    return cd


def bub_tree(tree, fasta, outfile1, root, types, c_dict, show, size,
             colours, field1, field2, scale, multiplier, dna):
    """
    :param tree: tree object from ete
    :param fasta: the fasta file used to make the tree
    :param outfile1: outfile suffix
    :param root: sequence name to use as root
    :param types: tree type: circular (c) or rectangle (r)
    :param c_dict: dictionary mapping colour to time point (from col_map)
    :param show: show the tree in a gui (y/n)
    :param size: scale the terminal nodes by frequency information (y/n)
    :param colours: if using a matched fasta file, colour the sequence by charge/IUPAC
    :param field1: the field that contains the size/frequency value
    :param field2: the field that contains the size/frequency value
    :param scale: how much to scale the x axis
    :param multiplier
    :param dna true/false, is sequence a DNA sequence?
    :param t_list list of time points
    :return: None, outputs svg/pdf image of the tree
    """

    if multiplier is None:
        mult = 500
    else:
        mult = multiplier

    if dna:
        dna_prot = 'dna'
        bg_c = {'A': 'green',
                'C': 'blue',
                'G': 'black',
                'T': 'red',
                '-': 'grey',
                'X': 'white'}

        fg_c = {'A': 'black',
                'C': 'black',
                'G': 'black',
                'T': 'black',
                '-': 'black',
                'X': 'white'}
    else:
        dna_prot = 'aa'
        bg_c = {'K': '#145AFF',
                'R': '#145AFF',
                'H': '#8282D2',
                'E': '#E60A0A',
                'D': '#E60A0A',
                'N': '#00DCDC',
                'Q': '#00DCDC',
                'S': '#FA9600',
                'T': '#FA9600',
                'L': '#0F820F',
                'I': '#0F820F',
                'V': '#0F820F',
                'Y': '#3232AA',
                'F': '#3232AA',
                'W': '#B45AB4',
                'C': '#E6E600',
                'M': '#E6E600',
                'A': '#C8C8C8',
                'G': '#EBEBEB',
                'P': '#DC9682',
                '-': 'grey',
                'X': 'white'}

        fg_c = {'K': 'black',
                'R': 'black',
                'H': 'black',
                'E': 'black',
                'D': 'black',
                'N': 'black',
                'Q': 'black',
                'S': 'black',
                'T': 'black',
                'L': 'black',
                'I': 'black',
                'V': 'black',
                'Y': 'black',
                'F': 'black',
                'W': 'black',
                'C': 'black',
                'M': 'black',
                'A': 'black',
                'G': 'black',
                'P': 'black',
                '-': 'grey',
                'X': 'white'}

    if colours == 3:
        bg_c = None
        fg_c = None

    # outfile3 = str(outfile1.replace(".svg", ".nwk"))

    tstyle = TreeStyle()
    tstyle.force_topology = False
    tstyle.mode = types
    tstyle.scale = scale
    tstyle.min_leaf_separation = 0
    tstyle.optimal_scale_level = 'full'  # 'mid'
    # tstyle.complete_branch_lines_when_necessary = False
    if types == 'c':
        tstyle.root_opening_factor = 0.25

    tstyle.draw_guiding_lines = False
    tstyle.guiding_lines_color = 'slateblue'
    tstyle.show_leaf_name = False
    tstyle.allow_face_overlap = True
    tstyle.show_branch_length = False
    tstyle.show_branch_support = False
    TreeNode(format=0, support=True)
    # tnode = TreeNode()

    if root is not None:
        tree.set_outgroup(root)
    # else:
    #     r = tnode.get_midpoint_outgroup()
    #     print("r", r)
    #     tree.set_outgroup(r)
    time_col = []
    for node in tree.traverse():
        # node.ladderize()
        if node.is_leaf() is True:
            try:
                name = node.name.split("_")
                time = name[field2]
                kind = name[3]
                # print(name)
            except:
                time = 'zero'
                name = node.name
                print("Incorrect name format for ", node.name)

            if size is True:
                try:
                    s = 20 + float(name[field1]) * mult
                except:
                    s = 20
                    print("No frequency information for ", node.name)
            else:
                s = 20

            colour = c_dict[time]
            time_col.append((time, colour))
            nstyle = NodeStyle()
            nstyle["fgcolor"] = colour
            nstyle["size"] = s
            nstyle["hz_line_width"] = 10
            nstyle["vt_line_width"] = 10
            nstyle["hz_line_color"] = colour
            nstyle["vt_line_color"] = 'black'
            nstyle["hz_line_type"] = 0
            nstyle["vt_line_type"] = 0
            node.set_style(nstyle)

            if root is not None and node.name == root:  # place holder in case you want to do something with the root leaf
                print('root is ', node.name)
                # nstyle["shape"] = "square"
                # nstyle["fgcolor"] = "black"
                # nstyle["size"] = s
                # nstyle["shape"] = "circle"
                # node.set_style(nstyle)

            else:
                nstyle["shape"] = "circle"
                node.set_style(nstyle)

            if fasta is not None:
                seq = fasta[str(node.name)]
                seqFace = SequenceFace(seq, seqtype=dna_prot, fsize=10, fg_colors=fg_c, bg_colors=bg_c, codon=None,
                                       col_w=40, alt_col_w=3, special_col=None, interactive=True)
                # seqFace = SeqMotifFace(seq=seq, motifs=None, seqtype=dna_prot, gap_format=' ', seq_format='()', scale_factor=20,
                #              height=20, width=50, fgcolor='white', bgcolor='grey', gapcolor='white', )
                # seqFace = SeqMotifFace(seq, seq_format="seq", fgcolor=fg_c, bgcolor=bg_c) #interactive=True

                (tree & node.name).add_face(seqFace, 0, "aligned")

        else:
            nstyle = NodeStyle()
            nstyle["size"] = 0.1
            nstyle["hz_line_width"] = 10
            nstyle["vt_line_width"] = 10
            node.set_style(nstyle)
            continue
    tree.ladderize()
    # tnode.ladderize()
    legendkey = sorted(set(time_col))
    legendkey = [(tp, col) for tp, col in legendkey]
    # legendkey.insert(0, ('Root', 'black'))
    legendkey.append(('', 'white'))

    for tm, clr in legendkey:
        tstyle.legend.add_face(faces.CircleFace(30, clr), column=0)
        tstyle.legend.add_face(faces.TextFace('\t' + tm, ftype='Arial', fsize=60, fgcolor='black', tight_text=True),
                               column=1)
    if show is True:
        tree.show(tree_style=tstyle)

    tree.render(outfile1, dpi=600, tree_style=tstyle)
    # tree.write(format=1, outfile=outfile3) #nwk tree file


def main(infile, fasta, outpath, name, root, types, show, size, colours, field1, field2, scale, multiplier, dna,
         consens):
    print(infile)

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    name = name + '.svg'  # #svg, pdf or png possible
    outfile = os.path.join(outpath, name)
    if fasta is not None:
        dc = fasta_to_dct(fasta)
        if colours != 3 and consens is None:
            print("Must supply sequence name to use as reference for highlighter/indels with the -con flag")
            sys.exit()
        elif colours != 3 and consens is not None:
            d = highlighter_dct(dc, consens, colours)
        else:
            d = dc
    else:
        d = None

    if scale == 0:
        print("scale of zero is silly, everything will disappear! Don't do this")
        sys.exit()

    t = Tree(infile)
    t_list = get_time(t, field2)
    cdict = col_map(t_list)

    bub_tree(t, d, outfile, root, types, cdict, show, size, colours, field1, field2, scale, multiplier, dna)

    print("We are done here")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='render a phylogenetic tree from a newark file and fasta '
                                                 'sequence (optional',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--infile', type=str,
                        help='The input newick tree file', default=argparse.SUPPRESS, required=True)
    parser.add_argument('-f', '--fasta', type=str,
                        help='The fasta file used to make the tree', required=False, default=None)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path for the output tree and image', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='The filename for the output tree and image (with no suffix: ".png")', required=True)
    parser.add_argument('-r', '--root', default=None, type=str,
                        help='the sequence name to root on', required=False)
    parser.add_argument('-con', '--cons', default=None, type=str,
                        help='the sequence to use as the consensus for highlighter plot', required=False)
    parser.add_argument('-t', '--types', type=str, default='r',
                        help='default is regular tree, for circular tree use "-t c"', required=False)
    parser.add_argument('-s', '--show', default=False, action='store_true',
                        help='show tree', required=False)
    parser.add_argument('-z', '--size', default=False, action='store_true',
                        help='do the sequence names have frequency information at last "_" separated field',
                        required=False)
    parser.add_argument('-c', '--colours', type=int, default=2,
                        help='Requires the -f flag, Changes the format of the sequence images: show indel blocks (1), '
                             'use highlighter plot (2), colour by IUPAC chemistry (3)\n'
                             'Options 1 and 2 require the -con flag',
                        required=False)
    parser.add_argument('-x', '--field1', default=-1, type=int,
                        help='the field containing the size/frequency info (based on "_" delimeter), '
                             'uses zero based indexing (ie field 1 = 0) default is "-1" ie: the last position',
                        required=False)
    parser.add_argument('-y', '--field2', default=2, type=int,
                        help='the field to be used to colour the leaf nodes (time, sample...), based in "_" delimiter, '
                             'using zero based indexing (field 1=0)',
                        required=False)
    parser.add_argument('-v', '--scale', default=10000, type=int,
                        help='how much to scale the x-axis, often need to increase in 10 fold increments, ',
                        required=False)
    parser.add_argument('-m', '--multiplier', default=500, type=int,
                        help='how much multiply the size/frequency value by',
                        required=False)
    parser.add_argument('-d', '--dna', default=False, action='store_true',
                        help='Use if the fasta file is a DNA sequence',
                        required=False)

    args = parser.parse_args()
    infile = args.infile
    fasta = args.fasta
    outpath = args.outpath
    name = args.name
    root = args.root
    types = args.types
    show = args.show
    size = args.size
    colours = args.colours
    field1 = args.field1
    field2 = args.field2
    scale = args.scale
    multiplier = args.multiplier
    dna = args.dna
    cons = args.cons

    main(infile, fasta, outpath, name, root, types, show, size, colours, field1, field2, scale, multiplier, dna, cons)
