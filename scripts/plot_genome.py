#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def read_fasta(fasta):
    seq = ''
    seq_id = ''
    b = 0
    with open(fasta, 'r') as f_fa:
        for line in f_fa:
            if line.startswith('>'):
                if b == 1:
                    yield seq_id, seq
                    seq = ''
                seq_id = line.split()[0].replace('>', '')
                b = 1
            else:
                seq += line.strip()
        yield seq_id, seq


def cal_x_y(seq_len):
    x = []
    y = []
    ccu = 0
    seq_len.sort(reverse=True)
    for i in seq_len:
        ccu += float(i)/1000000
        y.append(float(i)/1000000)
        x.append(ccu)
    return x, y


def calc_n50(numlist):
    numlist.sort(reverse=True)
    s = sum(numlist)
    sub = 0
    limit = s * 0.5
    for l in numlist:
        sub = sub + l
        if sub >= limit:
            return float(l)/1000000, float(s)/1000000


def plot_geno(args):
    seq_len = []

    for seq_id, seq in read_fasta(args.genome):
        seq_len.append(len(seq))

    x, y = cal_x_y(seq_len)
    n50, gsize = calc_n50(seq_len)

    plt.switch_backend('agg')
    fig, ax = plt.subplots()
    ax.plot(x, y, c='r')
    plt.tick_params(labelsize=12)
    ax.set_xlabel('Accumulated contig span (Mb)', size=15)
    ax.set_ylabel('Maximum contig size (Mb)', size=15)
    ax.annotate('N50 : {:,.2f} Mb\n\nMax : {:,.2f} Mb\n\nSize : {:,.2f} Mb'.format(float(n50), float(max(seq_len))/1000000, float(gsize)),
                xy=(x[-1]-x[-1]/3.5, y[0]-y[0]/4),
                size=12)
    fig.tight_layout()
    plt.savefig('%s.genome.pdf' % args.name)
    plt.savefig('%s.genome.png' % args.name)


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-g', '--genome', required=True,
                        help='genome file, .fasta format')
    parser.add_argument('-n', '--name', default='result', type=str,
                        help='prefix of output')
    args = parser.parse_args()

    plot_geno(args)


if __name__ == "__main__":

    main()
