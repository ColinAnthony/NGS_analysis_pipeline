import argparse
import sys, os
import collections
import math
from itertools import groupby
try:
    from ete3 import Tree, faces, TreeStyle, NodeStyle, TreeNode, add_face_to_node,  SequenceFace, SeqMotifFace
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


def get_colour_maps(tree, field2):
    """
    :param tree: tree object parsed from ete module
    :param field2: the '_' delimited field containing the sample/time point identifier
    :return: non_redundant list of time points, number of time points
    """
    colours = []
    for node in tree.traverse():
        if node.is_leaf() == True:
            name = node.name.split("_")
            try:
                t = name[field2]
            except:
                t = 'black'
                print("incorrect name format: sequence colour code will be black", node.name)

            if t not in colours:
                colours.append(t)

    colours.sort()

    return colours


def col_map(colours):
    """
    :param colours: non_redundant list of time points from get_times function
    :param n_color: number of time points
    :return: dictionary mapping colour spectrum to time points, key = time point, value = colour code
    """

    full_colour_list = ['#E50001', '#E32A00', '#E15700', '#E08300', '#DEAE00', '#DDD800', '#B4DB00', '#88D900',
              '#5CD800', '#31D600', '#06D500', '#00D323', '#00D24C', '#00D075', '#00CE9D', '#00CDC4',
               '#00ABCB', '#0082CA', '#0059C8', '#0031C6', '#000AC5', '#1C00C3', '#4200C2', '#6800C0', '#8D00BF']

    other = ['#000000', '#42c2f4']

    n_colors = len(colours)
    l = len(full_colour_list)
    if n_colors <= l:
        print(f"The colour gradient in this script can only handle 25 different colours, you require {l} colours\n"
              "duplicating colours: different sequences will have the same colour")
        required_cols = math.ceil(n_colors / l)
        full_colour_list = full_colour_list * required_cols

    cd = collections.OrderedDict()
    for i, v in enumerate(colours):
        if str(v) == 'black' or str(v) == root.split("_")[field2].lower():
            cd[str(v)] = other[0]
        elif n_colors == 1:
            cd[str(v)] = other[1]
        elif n_colors == 2:
            cd[str(v)] = other[i]
        else:
            cd[str(v)] = full_colour_list[i]

    return cd


def bub_tree(tree, fasta, outfile, root, types, c_dict, show, size,
             colours, field1, field2, scale, multiplier, dna):
    """
    :param tree: tree object from ete
    :param fasta: the fasta file used to make the tree
    :param outfile: outfile suffix
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
        bg_c = {'A': '#1b7837',
                'C': '#053061',
                'G': 'black',
                'T': '#b2182b',
                '-': 'grey',
                'X': 'white',
                'Y': 'grey',
                'R': 'grey',
                'W': 'grey'}

        fg_c = {'A': 'black',
                'C': 'black',
                'G': 'black',
                'T': 'black',
                '-': 'black',
                'X': 'white',
                'Y': 'grey',
                'R': 'grey',
                'W': 'grey'}
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

    tstyle = TreeStyle()
    tstyle.force_topology = False
    tstyle.mode = types
    tstyle.scale = scale
    tstyle.min_leaf_separation = 0
    tstyle.optimal_scale_level = 'full' #'mid'
    # tstyle.complete_branch_lines_when_necessary = False
    if types == 'c':
        tstyle.root_opening_factor = 0.25

    tstyle.draw_guiding_lines = False
    tstyle.guiding_lines_color = 'slateblue'
    tstyle.show_leaf_name = True
    tstyle.allow_face_overlap = True
    tstyle.show_branch_length = True
    tstyle.show_branch_support = False
    TreeNode(format=0, support=True)
    # tnode = TreeNode()
    # edit leaf names for cleaner labels
    # leaves = tree.get_leaves()
    # for leaf in leaves:
    #     leaf.name = leaf.name.split("_")[0] + "_" + leaf.name.split("_")[-2] + "_" + leaf.name.split("_")[-1]
    if root is not None:
        tree.set_outgroup(root)
    else:
        r = tree.get_midpoint_outgroup()
        tree.set_outgroup(r)
    time_col = []
    tree.ladderize()

    tree_branch_width = 2
    for node in tree.traverse():
        node.ladderize()
        if node.is_leaf() is True:
            name = node.name.split("_")
            color_by = name[field2]

            if size is True:
                try:
                    s = 20 + float(name[field1]) * mult
                except:
                    s = 20
                    print("No frequency information for ", node.name)
            else:
                s = 20
            try:
                colour = c_dict[color_by]
            except KeyError as e:
                print("colour code not found, check maing format and colour field,\nUsing black as the colour")
                colour = "#000000"

            time_col.append((color_by, colour))
            nstyle = NodeStyle()
            nstyle["shape"] = "circle"
            nstyle["fgcolor"] = colour
            nstyle["size"] = s
            nstyle["hz_line_width"] = tree_branch_width
            nstyle["vt_line_width"] = tree_branch_width
            nstyle["hz_line_color"] = colour
            nstyle["vt_line_color"] = 'black'
            nstyle["hz_line_type"] = 0
            nstyle["vt_line_type"] = 0
            node.set_style(nstyle)

            if root is not None and node.name == root: #place holder in case you want to do something with the root leaf
                print('root is ', node.name)
                #nstyle["shape"] = "square"
                # nstyle["fgcolor"] = "black"
                # nstyle["size"] = s
                # nstyle["shape"] = "circle"
                # node.set_style(nstyle)

            else:
                node.set_style(nstyle)

            if fasta is not None:
                seq = fasta[str(node.name)]
                seqFace = SequenceFace(seq, seqtype=dna_prot, fsize=2, fg_colors=fg_c, bg_colors=bg_c, codon=None,
                                       col_w=2, alt_col_w=1, special_col=None, interactive=True)
                # seqFace = SeqMotifFace(seq=seq, motifs=None, seqtype=dna_prot, gap_format=' ', seq_format='()', scale_factor=20,
                #              height=20, width=50, fgcolor='white', bgcolor='grey', gapcolor='white', )
                # seqFace = SeqMotifFace(seq, seq_format="seq", fgcolor=fg_c, bgcolor=bg_c) #interactive=True

                (tree & node.name).add_face(seqFace, 0, "aligned")

        else:
            nstyle = NodeStyle()
            nstyle["size"] = 0.1
            nstyle["hz_line_width"] = tree_branch_width
            nstyle["vt_line_width"] = tree_branch_width
            node.set_style(nstyle)
            continue

    legendkey = sorted(set(time_col))
    legendkey = [(tp, col) for tp, col in legendkey]
    # legendkey.insert(0, ('Root', 'black'))
    legendkey.append(('', 'white'))

    for tm, clr in legendkey:
        tstyle.legend.add_face(faces.CircleFace(30, clr), column=0)
        tstyle.legend.add_face(faces.TextFace('\t' + tm, ftype='Arial', fsize=60, fgcolor='black', tight_text=True), column=1)
    if show is True:
        tree.show(tree_style=tstyle)

    tree.render(outfile + ".svg", dpi=600, tree_style=tstyle)
    tree.render(outfile + ".png", dpi=600, tree_style=tstyle)
    tree.render(outfile + ".pdf", dpi=600, tree_style=tstyle)
    #tree.write(format=1, outfile=outfile3) #nwk tree file


def main(infile, fasta, outpath, name, root, types, show, size, colours, field1, field2, scale, multiplier, dna, consens):
    print(infile)

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    outfile = os.path.join(outpath, name)

    root = root.lower()
    t = Tree(infile)
    leaves = t.get_leaves()
    for leaf in leaves:
        leaf.name = leaf.name.lower()
    leaf_names = [leaf.name for leaf in leaves]

    if fasta is not None:
        dc = fasta_to_dct(fasta)
        new_dict = {}
        for name, seq in dc.items():
            new_dict[name.lower()] = seq
            if name.lower() not in leaf_names:
                print("The sequence name in your fasta file is not in the tree file:", name)
                sys.exit("exiting")
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



    t_list = get_colour_maps(t, field2)
    cdict = col_map(t_list)

    bub_tree(t, d, outfile, root, types, cdict, show, size, colours, field1, field2, scale, multiplier, dna)

    print("We are done here")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='render a phylogenetic tree from a newick file and fasta '
                                        'sequence (optional)', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--infile', type=str,
                        help='The input newick tree file', default=argparse.SUPPRESS, required=True)
    parser.add_argument('-f', '--fasta', type=str,
                        help='The fasta file used to make the tree', required=False, default=None)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path for the output tree and image', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='The filename for the output tree and image (with no suffix: ".png")', required=True)
    parser.add_argument('-r', '--root', default=None, type=str,
                        help='the sequence name to root on. If using an external root, it must include the '
                             'field1 and field 2 tags, at the same index positions as the rest of the sequences',
                        required=False)
    parser.add_argument('-con', '--cons', default=None, type=str,
                        help='the sequence to use as the consensus for highlighter plot', required=False)
    parser.add_argument('-t', '--types', type=str, default='r',
                        help='default is regular tree, for circular tree use "-t c"', required=False)
    parser.add_argument('-s', '--show', default=False,  action='store_true',
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
