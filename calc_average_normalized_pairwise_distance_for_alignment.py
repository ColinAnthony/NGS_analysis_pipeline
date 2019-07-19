import os
import argparse
from smallBixTools import smallBixTools as st


__author__ = 'David Matten'


def main(in_fasta, outpath):
    dct = st.fasta_to_dct(in_fasta)
    all_distances = []
    outfile = os.path.join(outpath, "pairwise_distances.csv")
    with open(outfile, "w") as handle:
        handle.write("seq_id1,seq_id2,Normalised_hamming_distance_adjusted_(changes_per_100_bases)\n")

        for k1, v1 in dct.items():
            for k2, v2 in dct.items():
                if k1 != k2:
                    dist = st.customdist(v1.upper(), v2.upper())
                    norm_dist = dist / len(v1)
                    normadjustdist_perc = round(norm_dist * 100, 2)

                    handle.write(",".join([str(x) for x in [k1, k2, normadjustdist_perc]]) + "\n")
                    all_distances.append(normadjustdist_perc)
    print(sum(all_distances)/len(all_distances))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate genetic distance from a reference sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--in_fasta', type=str, required=True,
                        help='The path to the input fasta file (ie: /path/to/fasta_file.fasta')
    parser.add_argument('-out', '--outpath', type=str, required=True,
                        help='The path to the where the outfile should be written')

    args = parser.parse_args()

    inpath = args.in_fasta
    outpath = args.outpath

    main(inpath, outpath)
