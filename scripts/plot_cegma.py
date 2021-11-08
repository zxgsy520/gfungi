#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_cegma(file):

    for line in open(file, 'r'):
        line = line.strip()

        if not line or line.startswith("#"):
                continue
        if line.startswith("KO") or line.startswith("WARNING"):
                continue

        line = line.replace('Group ', 'Group').split()

        yield [line[0], int(line[1]), float(line[2])]


def truth_cegma(file):

    cegma_dict = {}

    for line in read_cegma(file):
        if line[0] == "Complete":
            name = "Complete"
            line[0] = "Total"
        if line[0] == "Partial":
            name = "Partial"
            line[0] = "Total"
        if name not in cegma_dict:
            cegma_dict[name] = {}
            n = 0
        cegma_dict[name][n] = [line[0],line[1], line[2]]
        n += 1

    return cegma_dict


def stat_cegma(cegma_dict, out_name):

    comp = cegma_dict["Complete"]
    part = cegma_dict["Partial"]

    with open("%s.cegma.tsv" % out_name, 'w') as fh:
        fh.write('#\tComplete\tComplete+Partial\n')
        fh.write('#Type\tProts\t%completeness\tProts\t%completeness\n')

        for i in range(5):
            fh.write('{}\t{}\t{}\t{}\t{}\n'.format(comp[i][0], comp[i][1], comp[i][2],part[i][1], part[i][2]))

    return 0


def plot_cegma(cegma_dict ,out_name):

    group = []
    partial = []
    complete = []
    part_list = cegma_dict['Partial']
    comp_list = cegma_dict['Complete']

    for i in range(1,5):
        group.append(part_list[i][0])
        partial.append(part_list[i][2])
        complete.append(comp_list[i][2])

    fig = plt.figure(figsize=[9, 5])
    ax = fig.add_subplot(1,1,1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.subplots_adjust(left=0.10, right=0.60, top=0.95, bottom=0.08)
    x = range(len(group))

    ax.bar(x, partial, width=0.5, label='Complete+Partial:{:.2f}%'.format(part_list[0][2]), color='#CCFF99')
    ax.bar(x, complete, width=0.5, label='Complete:{:.2f}%'.format(comp_list[0][2]), color='#99CC99')
    plt.xticks(x, group)
    plt.legend(loc='center right', bbox_to_anchor=(1.60, 0.95) ,frameon=False)
    plt.ylabel('%completeness' ,fontsize=12)
    plt.ylim((0, 100))
    plt.savefig("%s.cegma.png" % out_name, dpi=700)
    plt.savefig("%s.cegma.pdf" % out_name)


def stat_cegma_hlep(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Enter the   evaluation result.')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default='out',
        help='Output file prefix,default=out.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    plot_cegma.py  Statistics and plotting CEGMA results.

attention:
    plot_cegma.py -i p.completeness_report
''')
    args = stat_cegma_hlep(parser).parse_args()
    cegma_dict = truth_cegma(args.input)

    stat_cegma(cegma_dict, args.out)
    plot_cegma(cegma_dict, args.out)


if __name__ == "__main__":

    main()
