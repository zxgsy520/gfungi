#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from collections import OrderedDict


LOG = logging.getLogger(__name__)
__version__ = "1.1.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def stat_genome_size(file):

    base = 0

    for line in open(file):

        if line.startswith(">"):
            continue

        base += len(line.strip())

    return base


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def split_attr(attributes):

    r ={}

    for content in attributes.split(';'):
        if not content:
            continue
        if '=' not in content:
            print('%r is not a good formated attribute: no tag!')
            continue
        tag, value = content.split('=', 1)
        r[tag] = value

    return r


def stat_ncRNA(gff, genome):

    rna_dict = OrderedDict()
    rna_dict["regulatory"] = {"all": [0, 0]}
    rna_dict["tRNA"] = {"all": [0, 0]}
    rna_dict["rRNA"] = {"18S": [0, 0], "28S": [0, 0], "5.8S": [0, 0],
                        "5S": [0, 0], "all": [0, 0]}
    rna_dict["ncRNA"] = {"snRNA": [0, 0], "miRNA": [0, 0], "all": [0, 0],
                         "spliceosomal": [0, 0], "other": [0, 0]}

    genome_size = stat_genome_size(genome)

    for line in read_tsv(gff, "\t"):
        start, end = int(line[3]), int(line[4])
        if line[2] in rna_dict:
            rna_dict[line[2]]["all"][0] += 1
            rna_dict[line[2]]["all"][1] += end-start+1

        attr = split_attr(line[-1])
        if line[2] == "rRNA":
            rrnatype = attr["product"].split()[0]
            if rrnatype in rna_dict["rRNA"]:
                rna_dict["rRNA"][rrnatype][0] += 1
                rna_dict["rRNA"][rrnatype][1] += end-start+1
        elif line[2] == "ncRNA":
            if "microRNA" in attr["product"]:
                rna_dict["ncRNA"]['miRNA'][0] += 1
                rna_dict["ncRNA"]['miRNA'][1] += end-start+1
            elif 'mall nucleolar RNA' in attr["product"]:
                rna_dict["ncRNA"]['snRNA'][0] += 1
                rna_dict["ncRNA"]['snRNA'][1] += end-start+1
            elif 'spliceosomal' in attr["product"]:
                rna_dict["ncRNA"]['spliceosomal'][0] += 1
                rna_dict["ncRNA"]['spliceosomal'][1] += end-start+1
            else:
                rna_dict["ncRNA"]['other'][0] += 1
                rna_dict["ncRNA"]['other'][1] += end-start+1
        else:
            continue

    print("#Types\tType\tNumber\tAverage length(bp)\tTotal length(bp)\tPercentage(%)")
    for i in rna_dict:
        for j in rna_dict[i]:
            number, total = rna_dict[i][j]
            average = 0
            if number:
                average = total*1.0/number

            print("{0}\t{1}\t{2:,}\t{3:,.2f}\t{4:,}\t{5:.4f}".format(
                i, j, number, average, total, total*100.0/genome_size)
            )


    return 0


def add_help_args(parser):

    parser.add_argument('gff',  metavar='FILE', type=str,
        help='Input RNA annotation results (gff).')
    parser.add_argument('-g', '--genome',  metavar='FILE', type=str, required=True,
        help='Input genome file (fasta).')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    stat_ncRNA: stat ncRNA result.
attention:
    stat_ncRNA RNA.gff -g genome.fasta >stat_ncRNA.tsv
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_help_args(parser).parse_args()
    stat_ncRNA(args.gff, args.genome)


if __name__ == "__main__":
    main()
