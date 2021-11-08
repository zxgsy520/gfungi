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


def read_paf(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        line = line.split('\t')
        match = int(line[9])
        blen = int(line[10])
        match = min(blen, match)

        if (match*100.0/blen)>50 or int(line[11])>0:
            yield line


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
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield [seq[0], seq[1]]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield [seq[0], seq[1]]
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
        line = line.strip()

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq[0], seq[1]
            seq = []
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

    if len(seq)==4:
        yield seq[0], seq[1]
    fp.close()


def choose_map_reads(file, pafile, name, base):

    read = set()
    n = 0
    output = open("%s.choose.fa" % name, "w")

    if base=="all":
        nbase = "all"
        base = 0
    else:
        nbase = ""
        base = int(base)

    for line in read_paf(pafile):
        read.add(line[0])

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for seqid,seq in fh:
        if seqid not in read:
            continue
        if n>=base and nbase!="all":
            break
        output.write('>%s\n%s\n' % (seqid, seq))
        n+=len(seq)
    output.close()


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, required=True,
        help="Input files in paf")
    parser.add_argument("-r", "--read", metavar='FILE', type=str, required=True,
        help="Input files in fastq or fasta format")
    parser.add_argument("-b", "--base", metavar='STR', type=str, default='all',
        help="Selected data volume, default=all")
    parser.add_argument("-n", "--name", metavar='STR', type=str, default='out',
        help="The name of the output file")
    
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
    choose_unmap_reads.py  Select no reads on the map

attention:
    choose_unmap_reads.py -i paf -r fastq -n unmap

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    choose_map_reads(args.read, args.input, args.name, args.base)


if __name__ == "__main__":

    main()
