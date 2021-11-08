#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys
import gzip
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


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
            line = line.decode("utf-8")
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
            line = line.decode("utf-8")
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


def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


def allot_base(file, base):

    name = file.split("/")[-1].split(".")[0]

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    n = 0
    lable = 1
    fname = "%s_%d.fa" % (name, lable)
    output = open(fname, 'w')

    for seqid, seq in fh:
        if n>=base:
            output.close()
            fname = check_path(fname)
            print(fname)
            n = len(seq)
            lable += 1
            fname = "%s_%d.fa" % (name, lable)
            output = open(fname, 'w')
            output.write(">%s\n%s\n" % (seqid, seq))
        else:
            output.write(">%s\n%s\n" % (seqid, seq))
            n += len(seq)

    output.close()
    fname = check_path(fname)
    print(fname)


def add_help(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input fasta file.')
    parser.add_argument('-b', '--base', metavar='INT', type=int, default=10000000,
        help='Split base size, default=10000000.')

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
    allot_base.py Split sequence file based on base amount

attention:
    allot_base.py fastq
    allot_base.py fasta -b 10000000 >split.list

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    allot_base(args.input, args.base)


if __name__ == "__main__":

    main()
