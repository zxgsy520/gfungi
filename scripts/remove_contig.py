#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)
    LOG.info("reading message from %r" % file)

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fp.close()


def read_gc_depth(file):

    lengths = 0
    total_gc = 0
    total_depth = 0
    r = {}

    for line in read_tsv(file, "\t"):
        lengths += float(line[1])
        total_gc += float(line[1])*float(line[2])
        total_depth += float(line[1])*float(line[3])
        r[line[0]] = [float(line[1]), float(line[2]), float(line[3])]

    return r, total_gc/lengths, total_depth/lengths


def pick_gc_depth(file, minlen=10000, diff_gc=0.5, diff_depth=0.5):

    data, aver_gc, aver_depth = read_gc_depth(file)
    maxgc = aver_gc*(1+diff_gc)
    mingc = aver_gc*(1-diff_gc)
    maxdepth = aver_depth*(1+diff_depth)
    mindepth = aver_depth*(1-diff_depth)
    LOG.info("Average GC content:%s\tAverage depth:%s" % (aver_gc, aver_depth))
    LOG.info("Remove ID\tLength\tGC(%)\tCoverage")

    r = []

    for i in data:
        line = data[i]
        if line[0] <= minlen:
            LOG.info("{0}\t{1:,}\t{2:.2f}\t{3:,.2f}".format(i, line[0], line[1], line[2]))
            continue
        if line[1] >= maxgc or line[1] <= mingc:
            LOG.info("{0}\t{1:,}\t{2:.2f}\t{3:,.2f}".format(i, line[0], line[1], line[2]))
            continue
        if line[2] >=maxdepth  or line[2] <= mindepth:
            LOG.info("{0}\t{1:,}\t{2:.2f}\t{3:,.2f}".format(i, line[0], line[1], line[2]))
            continue
        r.append(i)

    return r


def split_attr(attributes):

    r = {}

    for content in attributes.split(';'):
        if not content:
            continue
        if '=' not in content:
            LOG.info('%r is not a good formated attribute: no tag!')
            continue
        tag, value = content.split('=', 1)
        r[tag] = value

    return r


def remove_cont(file):

    cont_rna = ["16S_rRNA", "23S_rRNA", "12S_rRNA"]
    r = set()

    for line in read_tsv(file, "\t"):
        if len(line) <=8:
            continue
        attr = split_attr(line[-1])
        if ("note" in attr) or ("aligned only" in line[-1]):
            continue
        if attr["Name"] in cont_rna:
            LOG.info('Remove pollution sequence id:%s' % line[0])
            r.add(line[0])
    return r


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')

        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            if seq:
                yield seq.split('\n')
            seq = "%s\n" % line
        else:
            seq += line
    if seq:
        yield seq.split('\n')
    fp.close()


def read_fastq(file):
    '''Read fastq file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
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
        if not seq:
            seq.append(line.strip("@"))
            continue
        seq.append(line)

        if len(seq)==4:
            yield seq
            seq = []
    fp.close()


def remove_contig(file, gc_depth, rna, minlen=10000, diff_gc=0.5, diff_depth=0.5):

    seqids = pick_gc_depth(gc_depth, minlen=10000, diff_gc=0.5, diff_depth=0.5)
    rmids = remove_cont(rna)

    if file.endswith(".fastq.gz") or file.endswith(".fq.gz") or file.endswith(".fastq") or file.endswith(".fq"):
        fp = read_fastq(file)
    elif file.endswith(".fasta.gz") or file.endswith(".fa.gz") or file.endswith(".fasta") or file.endswith(".fa"):
        fp = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for record in fp:
        if record[0] in rmids:
            continue
        if record[0] in seqids:
            print(">%s\n%s" % (record[0], record[1]))
    return 0


def add_args(parser):

    parser.add_argument('genome',
        help='Input genome file(fasta, fastq)')
    parser.add_argument('-gd', '--gc_depth', metavar='FILE', type=str, required=True,
        help='Input statistics files of genome length, GC content and sequencing depth.')
    parser.add_argument('-r', '--rna', metavar='FILE', type=str, required=True,
        help='Input genomic RNA prediction result file(gff).')
    parser.add_argument('--minlen', metavar='INT', type=int, default=20000,
        help='Input the value to remove the minimum length, default=20000.')
    parser.add_argument('-gc', '--diff_gc', metavar='FLOAT', type=float, default=0.5,
        help='Input removed GC difference ratio, default=0.5.')
    parser.add_argument('-depth', '--diff_depth', metavar='FLOAT', type=float, default=0.5,
        help='Input removal depth difference ratio, default=0.5.')


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
    remove_contig.py -- Remove remove sequence

attention:
    remove_contig.py genome.fasta --gc_depth length_gc.tsv --rna rna.gff>genome_new.fasta
''')
    args = add_args(parser).parse_args()

    remove_contig(args.genome, args.gc_depth, args.rna, args.minlen, args.diff_gc, args.diff_depth)


if __name__ == "__main__":
    main()
