from pysam import VariantFile
from collections import namedtuple
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import argparse

def read_chrom_sizes(file):
    ret = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip().split()
            ret[line[0]] = int(line[1])
    return ret


def read_centromere(file):
    ret = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip().split()
            # "acen" refers to centromere
            # "stalk" refers to the short arm of acrocentric chromosomes, chr13,14,15,21,22;
            # "gvar" bands tend to be heterochomatin, either pericentric or telomeric.
            if line[-1] not in ["acen", "gvar", "stalk"]:
                continue
            chr, start, end = line[0], int(line[1]), int(line[2])
            if chr not in ret:
                ret[chr] = [start, end]
            else:
                if ret[chr][0] > start:
                    ret[chr][0] = start
                if ret[chr][1] < end:
                    ret[chr][1] = end
    return ret


class VCF(object):
    def __init__(self, file, chrom_sizes_file, centromere_file):
        self.file_name = file
        self.var_dict = self._read_variant_file(file)
        self.chrom_sizes = read_chrom_sizes(chrom_sizes_file)
        self.centromere = read_centromere(centromere_file)

    def _read_variant_file(self, file):
        var_dict = {}
        SNP = namedtuple("SNP", ['chrom', 'pos', 'gt'])
        handle = VariantFile(file, 'r')
        for rec in handle.fetch():
            snp = SNP(rec.chrom, rec.pos, rec.samples.values()[0]['GT'])
            if rec.chrom not in var_dict:
                var_dict[rec.chrom] = [snp]
            else:
                var_dict[rec.chrom].append(snp)
        return var_dict

    def draw_variants(self, save=None):
        chrom_names = [*self.chrom_sizes][:-1]
        #chrom_names = [*self.var_dict]

        # set up axes, legend
        fig, ax = plt.subplots(1, 1, figsize=(25, 25))
        ax.set_axis_off()
        lines = [Line2D([], [], color=c, marker='|', linestyle='None', markersize=20, markeredgewidth=2) for c in ["blue", "red"]]
        labels = ['Ref allele', 'Alt allele']
        ax.legend(lines, labels, loc="center right", fontsize=20)

        # set track positions
        x0, y0, width, height = 0, 0, 1, 1
        n_tracks = len(chrom_names)
        track_height = 0.9 * height / n_tracks * 0.35
        track_space = 0.9 * height / n_tracks * 0.05
        track_start_x = x0 + 0.1 * width
        text_start_x = x0 + 0.04*width
        text_a_start_x = x0 + 0.09*width
        track_start_y = y0 + 0.95 * height

        # draw each chromosome tracks
        i = 0
        for chr in chrom_names:
            track_width = 0.8 * width * self.chrom_sizes[chr] / 250000000   # all chrom normalized to length of 250 Mb
            centro_start = track_start_x + track_width * self.centromere[chr][0] / self.chrom_sizes[chr]
            centro_width = track_width * (self.centromere[chr][1] - self.centromere[chr][0]) / self.chrom_sizes[chr]
   
            # draw haplotype 1
            ## draw chromosome frames
            track_start_y_i = track_start_y - track_height*(2*i+1) - track_space*5*i
            text_start_y_i = track_start_y_i - 0.25*track_height
            text_a_start_y_i = track_start_y_i + 0.5*track_height
            rect_chr = Rectangle((track_start_x, track_start_y_i), track_width, track_height, linewidth=1, fill=None, alpha=1)
            rect_centro = Rectangle((centro_start, track_start_y_i), centro_width, track_height, linewidth=1, color='darkgrey', alpha=1)
            ax.add_patch(rect_chr)
            ax.add_patch(rect_centro)
            ax.text(text_start_x, text_start_y_i, str(chr), ha="left", va="center", fontsize=20)
            ax.text(text_a_start_x, text_a_start_y_i, 'a ', ha="left", va="center", fontsize=15)
            # draw each variant
            for j, var in enumerate(self.var_dict[chr]):
                if j % 100 != 0: continue  # draw the 1st SNP of every 100 SNPs, otherwise takes too long
                color = 'blue' if var.gt[0] == 0 else 'red'
                line_x = track_start_x + track_width * var.pos / self.chrom_sizes[chr]
                line = Line2D([line_x, line_x], [track_start_y_i, track_start_y_i + track_height*0.95], color=color, linewidth=0.2)
                ax.add_line(line)

            # draw haplotype 2
            track_start_y_i = track_start_y - track_height*(2*i+2) - track_space*(5*i+1)  # different from haplotype 1
            text_start_y_i = track_start_y_i + 0.45*track_height
            text_b_start_y_i = track_start_y_i + 0.5*track_height
            ax.text(text_a_start_x, text_b_start_y_i, 'b ', ha="left", va="center", fontsize=15)
            rect_chr = Rectangle((track_start_x, track_start_y_i), track_width, track_height, linewidth=1, fill=None, alpha=1)
            rect_centro = Rectangle((centro_start, track_start_y_i), centro_width, track_height, linewidth=1, color='darkgrey', alpha=1)
            ax.add_patch(rect_chr)
            ax.add_patch(rect_centro)
            ## draw each variant
            for j, var in enumerate(self.var_dict[chr]):
                if j % 100 != 0: continue
                color = 'blue' if var.gt[1] == 0 else 'red'        # different from haplotype 1
                line_x = track_start_x + track_width * var.pos / self.chrom_sizes[chr]
                line = Line2D([line_x, line_x], [track_start_y_i, track_start_y_i + track_height*0.95], color=color, linewidth=0.2)
                ax.add_line(line)
            i += 1

        if save is not None:
            plt.savefig(save, dpi=300, format="pdf")
        else:
            plt.show()


class CNV(object):
    def __init__(self, file, chrom_sizes_file, centromere_file):
        self.file_name = file
        self.var_dict = self._read_variant_file(file)
        self.chrom_sizes = read_chrom_sizes(chrom_sizes_file)
        self.centromere = read_centromere(centromere_file)

    def _read_variant_file(self, file):
        var_dict = {}
        CNV_T = namedtuple("CNV_T", ['chrom', 'start', 'end', 'cn'])
        with open(file, 'r') as f:
            for rec in f:
                rec = rec.strip().split()
                cnv = CNV_T(rec[0], int(rec[1]), int(rec[2]), rec[3])
                if cnv.chrom not in var_dict:
                    var_dict[cnv.chrom] = [cnv]
                else:
                    var_dict[cnv.chrom].append(cnv)
        return var_dict

    def draw_variants(self, save=None):
        chrom_names = [*self.chrom_sizes][:-1]
        color_dict = {'0': "slateblue", '1': "lightblue", '2': "white", '3': "lightpink", '4':"hotpink", "> 4":"red"}
        # set up axes, legend
        fig, ax = plt.subplots(1, 1, figsize=(25, 25))
        ax.set_axis_off()
        patches = [Patch(color=v, label="CN "+k) for k, v in color_dict.items()]
        ax.legend(handles=patches, loc="center right", fontsize=30)

        # set track positions
        x0, y0, width, height = 0, 0, 1, 1
        n_tracks = len(chrom_names)
        track_height = 0.9 * height / n_tracks * 0.7
        track_space = 0.9 * height / n_tracks * 0.3
        track_start_x = x0 + 0.1 * width
        text_start_x = x0 + 0.05*width
        track_start_y = y0 + 0.95 * height

        # draw each chromosome tracks
        for i, chr in enumerate(chrom_names):
            track_width = 0.8 * width * self.chrom_sizes[chr] / 250000000   # all chrom normalized to length of 250 Mb
            centro_start = track_start_x + track_width * self.centromere[chr][0] / self.chrom_sizes[chr]
            centro_width = track_width * (self.centromere[chr][1] - self.centromere[chr][0]) / self.chrom_sizes[chr]
   
            ## draw chromosome frames
            track_start_y_i = track_start_y - track_height*(i+1) - track_space*i
            text_start_y_i = track_start_y_i + 0.5*track_height
            rect_chr = Rectangle((track_start_x, track_start_y_i), track_width, track_height, linewidth=1, fill=False, edgecolor="grey", alpha=1)
            rect_centro = Rectangle((centro_start, track_start_y_i), centro_width, track_height, linewidth=1, facecolor='darkgrey', alpha=1)
            ax.add_patch(rect_chr)
            ax.add_patch(rect_centro)
            ax.text(text_start_x, text_start_y_i, str(chr), ha="left", va="center", fontsize=20)
            # draw each variant
            for var in self.var_dict[chr]:
                color = color_dict[var.cn] if var.cn in color_dict else color_dict["> 4"]
                rect_var_x = track_start_x + track_width * var.start / self.chrom_sizes[chr]
                rect_var_x_width = track_width * (var.end - var.start) / self.chrom_sizes[chr]
                rect_var = Rectangle((rect_var_x, track_start_y_i), rect_var_x_width, track_height, linewidth=1, edgecolor=None, facecolor=color, alpha=1)
                ax.add_patch(rect_var)

        if save is not None:
            plt.savefig(save, dpi=300, format="pdf")
        else:
            plt.show()


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--type", help="Types of variants to be drawn, currently supports snv or cnv")
    parser.add_argument("--chromsizes", help="A file with first column to be chromosome names, and second column to be chromosome sizes")
    parser.add_argument("--cytoband", help="A file that specifies the genomic coordinates of cytobands, can be downloaded from UCSC table browser")
    parser.add_argument("-o", help="Name for the output PDF file")
    parser.add_argument("input", help="A VCF file for SNV, or a BED file for CNV")
    return parser.parse_args()


def main():
    args = get_args()
    chrom_sizes = args.chromsizes
    centromere = args.cytoband
    input_file = args.input
    out_file = args.o

    if args.type == "snv":
        var = VCF(input_file, chrom_sizes, centromere)
    elif args.type == "cnv":
        var = CNV(input_file, chrom_sizes, centromere)
    else:
        raise ValueError("--type should be snv or cnv, currently supports plotting SNVs or CNVs")

    var.draw_variants(save=out_file)

if __name__ == "__main__":
    main()
