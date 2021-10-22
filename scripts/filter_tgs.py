#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import pysam
import logging
import argparse
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt


LOG = logging.getLogger(__name__)

__version__ = "1.2.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fa = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fa = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''
    for line in fa:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            seq = seq.split('\n')
            if len(seq)==2:
                yield seq[0], seq[1]
            seq = ''
            line = line.strip(">").split()[0]
            seq += "%s\n" % line
            continue
        seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]
    fa.close()


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file, 'r')
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith("@") and (len(seq)==0 or len(seq)>=5):
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith("@") and len(seq)==4:
            yield seq[0], seq[1]
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        seq.append(line)

    if len(seq)==4:
        yield seq[0], seq[1]
    fp.close()


def read_bam(file):

    fp = pysam.AlignmentFile(file, 'rb', check_sq=False)

    LOG.info('process %r' % file)
    for record in fp:
        yield record.qname, record.seq


def filter_reads(files, name, minlen=1000):

    raw_length = []
    filter_length = []
    output = open('%s.clean.fasta' % name, 'w')

    for file in files:
        if file.endswith('.fa') or file.endswith('.fasta') or file.endswith('.fa.gz') or file.endswith('.fasta.gz'):
            fh = read_fasta(file)
        elif file.endswith('.fq') or file.endswith('.fastq') or file.endswith('.fq.gz') or file.endswith('.fastq.gz'):
            fh = read_fastq(file)
        elif file.endswith('.bam'):
            fh = read_bam(file)
        else:
            raise Exception('Unknown format!')

        for seqid, seq in fh:
            seqlen = len(seq)

            raw_length.append(seqlen)
            if seqlen<=minlen:
                continue
            filter_length.append(seqlen)
            output.write('>%s\n%s\n' %(seqid, seq))
    output.close()

    return raw_length, filter_length


def n50(lengths):

    sumlen = sum(lengths)
    acculen = 0

    for i in sorted(lengths, reverse=True):
        acculen += i
        if acculen >= sumlen*0.5:
            break

    return i


def stat_len(lengths):

    k10 = 0
    k20 = 0
    k40 = 0

    for i in sorted(lengths, reverse=True):
        if i<10000:
            break
        k10 += 1
        if i >20000:
            k20 += 1
        if i >40000:
            k40 += 1
    return k10, k20, k40


def plot_read_length(lengths, name, xmin=0, xmax=35000, bins=100):

    #plt.style.use('ggplot')
    plt.switch_backend('agg')
    fig, ax = plt.subplots(figsize=(10, 6),)
#    ax.grid(color='w')
#    ax.patch.set_facecolor("#ffcfdc")
#    ax.spines['top'].set_visible(False) #去掉上边框
#    ax.spines['bottom'].set_visible(False) #去掉下边框
#    ax.spines['left'].set_visible(False) #去掉左边框
#    ax.spines['right'].set_visible(False) #去掉右边框

    ax.hist(lengths, bins=bins, range=(xmin, xmax), histtype='stepfilled', alpha=0.75)

    ax.set_xlim(xmin, xmax)
    font = {'weight': 'bold','size': 12,}
    ax.set_ylabel('Read Count', font)
    ax.set_xlabel('Read Length', font)
    plt.xticks()
    plt.savefig("%s.reads_length.pdf" % name)
    plt.savefig("%s.reads_length.png" % name, dpi=700)


def add_hlep(parser):
    parser.add_argument('-i', '--input', metavar='FILE', nargs='+', type=str, required=True,
        help='Input reads file(fasta, fatsq, bam).')
    parser.add_argument('--minlen', metavar='INT', type=int, default=1000,
        help='Set the minimum length of filtered reads, default=1000.')
    parser.add_argument('--xmin', metavar="INT", type=int, default=0,
        help='Minimum number of x axis, default=0.')
    parser.add_argument('--xmax', metavar="INT", type=int, default=80000,
        help='Maximum number of x axis, default=80000.')
    parser.add_argument('--bins', metavar="INT", type=int, default=200,
        help='Set the number of bins in the histogram, default=200.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default="out",
        help='The prefix name of the output file, default=out')

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
    filter_tgs.py  Quality control and statistics on three generations of data

attention:
    filter_tgs -i *.bam
    filter_tgs -i *.fa
    filter_tgs -i *.fq
version: %s
contact:  %s <%s>\
    ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep(parser).parse_args()
    raw_length, filter_length = filter_reads(args.input, args.name, args.minlen)

    plot_read_length(raw_length, args.name, args.xmin, args.xmax, args.bins)

    k10, k20, k40 = stat_len(raw_length)
    fk10, fk20, fk40 = stat_len(filter_length)
    raw_bases = sum(raw_length)
    filter_bases = sum(filter_length)
    raw_number = len(raw_length)
    filter_number = len(filter_length)
    output = open('%s.reads_stat.tsv' % args.name, 'w')

    output.write("""\
#Type\tBases(bp)\tReads number\tReads mean length(bp)\tReads max length(bp)\tReads N50 (bp)\tReads >10kb ratio(%)\tReads >20kb ratio(%)\tReads >40kb ratio(%)
Raw Reads\t{0:,}\t{1:,}\t{2:,.2f}\t{3:,}\t{4:,}\t{5:.2f}\t{6:.2f}\t{7:.2f}
Filtered Reads\t{8:,}\t{9:,}\t{10:,.2f}\t{11:,}\t{12:,}\t{13:.2f}\t{14:.2f}\t{15:.2f}
""".format(
        raw_bases, raw_number, raw_bases*1.0/raw_number, max(raw_length), n50(raw_length), k10*100.0/raw_number, k20*100.0/raw_number, k40*100.0/raw_number,
        filter_bases, filter_number, filter_bases*1.0/filter_number, max(filter_length), n50(filter_length), fk10*100.0/raw_number, fk20*100.0/raw_number, fk40*100.0/raw_number
    ))
    output.close()


if __name__ == "__main__":

    main()
