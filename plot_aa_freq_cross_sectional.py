import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


__author__ = 'colin anthony and David Matten'



def main(infile, outpath, treatment_grp):

    df = pd.read_csv(infile, sep=',', header=0)

    if outpath == None:
        outpath = os.getcwd()

    selection_bound = 5
    x_header = "participant"
    group_col = "HXB2_position"
    y_lim_min = 0
    y_lim_max = 105
    # xmin = 0


    for pos, df_grp in df.groupby(group_col):
        print(f"pos: {pos}")
        df_grp.drop(group_col, axis=1, inplace=True)
        aminos = list(df_grp)[1:]
        for resi in aminos:
            if df_grp[resi].min() > (100 - selection_bound):
                df_grp.drop(resi, axis=1, inplace=True)
            elif df_grp[resi].max() < selection_bound:
                df_grp.drop(resi, axis=1, inplace=True)

        my_colors = ["#C68A3B", "#7394CA", "#579F6A", "#C47BB2", "#D36C6E", "#53A7A6", "#8D909B", "#D6C463", "#A98B6D",
                     "grey"]

        aa_dict = {'A': '#C8C8C8',
                   'R': '#145AFF',
                   'N': '#00DCDC',
                   'D': '#E60A0A',
                   'C': '#E6E600',
                   'Q': '#00DCDC',
                   'E': '#E60A0A',
                   'G': '#EBEBEB',
                   'H': '#8282D2',
                   'I': '#0F820F',
                   'L': '#0F820F',
                   'K': '#145AFF',
                   'M': '#E6E600',
                   'F': '#3232AA',
                   'P': '#DC9682',
                   'S': '#FA9600',
                   'T': '#FA9600',
                   'W': '#B45AB4',
                   'Y': '#3232AA',
                   'V': '#0F820F',
                   '-': '#707070',
                   'X': '#000000'}

        df_new = df_grp
        # aa_dict = {
        #         'A': '#99ff99',
        #         'R': '#0000ff',
        #         'N': '#ff6666',
        #         'D': '#660000',
        #         'C': '#ffff66',
        #         'Q': '#ff3333',
        #         'E': '#660000',
        #         'G': '#000000',
        #         'H': '#6666ff',
        #         'I': '#006600',
        #         'L': '#336633',
        #         'K': '#3333cc',
        #         'M': '#cc9933',
        #         'F': '#666633',
        #         'P': '#666666',
        #         'S': '#ff6633',
        #         'T': '#cc3300',
        #         'W': '#666600',
        #         'Y': '#996633',
        #         'V': '#ff99ff',
        #         '-': '#707070',
        #         'X': '#ff00ff'}

        # if the site is not conserved
        if len(df_new.columns) > 2:

            headers = list(df_new)
            melt_values = headers[1:]
            df_plot = pd.melt(df_new, id_vars=['participant'], value_vars=melt_values)
            df_plot.rename(columns={"variable": "residue"}, inplace=True)

            fig, ax = plt.subplots()

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()

            # set the axis limits and axis tick intervals
            plt.yticks(np.arange(y_lim_min, y_lim_max, 20), fontsize=10)
            plt.ylim(0, y_lim_max)
            # plt.xticks([])

            sns.set_style("white")
            sns.set_style("ticks")
            sns.axes_style("white")

            p = sns.stripplot(x="participant", y="value", data=df_plot, size=10, hue="residue",
                              palette=my_colors, alpha=0.80)

            sns.despine(left=False, bottom=False)

            ax.set_ylabel('Frequency (%)', fontsize=12)
            ax.set_xlabel('', fontsize=12)

            # p.set(xticklabels=[])

            plt.title("Position " + str(int(pos)), x=0.0, y=1.01, fontsize=16)

            # add figure legend
            ax.legend_.remove()

            handles2, labels2 = ax.get_legend_handles_labels()
            plt.legend(handles2, labels2, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.legend(frameon=False)

            # set figure size
            f = plt.gcf()
            w = 6.875
            h = 4
            f.set_size_inches(w, h)

            # write figure to file
            out_fn = os.path.join(outpath, "Site_" + str(pos) + "_rel.png")
            plt.savefig(out_fn, ext='png', dpi=600, format='png', facecolor='white', bbox_inches='tight')
            plt.close('all')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot amino acid frequencies for each position, from a .csv file'
                                                 'format output from calc_entropy_aa_glycan_freq.py.')

    parser.add_argument('-in', '--infile', type=str,
                        help='The input .csv file', required=True)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=True)
    parser.add_argument('-s', '--treatment_grp', type=str,
                        help='The csv file with the seq id and treatment assignment', required=True)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    treatment_grp = args.treatment_grp

    main(infile, outpath, treatment_grp)
