#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):
    """
    read tsv joined with sep
    :param file: file name
    :param sep: separator
    :return: list
    """
    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_paf(file, identities=50, mapq=0, depth='all'):

    data = {}
    read = set()

    for line in read_tsv(file, '\t'):
        match = int(line[9])
        blen = int(line[10])
        match = min(blen, match)

        if (match*100.0/blen)<=identities and int(line[11])<=mapq:
            continue
        if depth=='all':
            read.add(line[0])
            continue
        refid = line[5]
        reflen = int(line[6])

        if refid not in data:
            data[refid] = 0

        if data[refid]>int(depth)*reflen:
            continue
        data[refid] += int(line[1])
        read.add(line[0])

    return read


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
        if type(line) == type(b''):
            line = line.decode('utf-8')
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
        if type(line) == type(b''):
            line = line.decode('utf-8')
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


def choose_map_reads(file, pafile, name, identities, mapq, depth):

    output = open("%s.choose.fa" % name, "w")
    read = read_paf(pafile, identities, mapq, depth)

    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".fasta.gz") or file.endswith(".fa.gz"):
        fh = read_fasta(file)
    else:
        raise Exception("%r file format error" % file)

    for seqid,seq in fh:
        if seqid not in read:
            continue
        output.write('>%s\n%s\n' % (seqid, seq))
        read.remove(seqid)
        if len(read)==0:
            break

    output.close()


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, required=True,
        help="Input files in paf")
    parser.add_argument("-r", "--read", metavar='FILE', type=str, required=True,
        help="Input files in fastq or fasta format")
    parser.add_argument("-id", "--identities", metavar='FLOAT', type=float, default=50.0,
        help="Set the minimum identities of reads comparison, default=50.0")
    parser.add_argument("-mq", "--mapq", metavar='INT', type=int, default=2,
        help="Set the minimum quality value for reads alignment, default=2")
    parser.add_argument("-d", "--depth", metavar='STR', type=str, default=200,
        help="Set the reads depth selected by each contig(Select all=all), default=200")
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
    choose_map_reads.py  Choose appropriate reads for subsequent correction

attention:
    choose_map_reads.py -i paf -r fastq -n map

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    choose_map_reads(args.read, args.input, args.name, args.identities, args.mapq, args.depth)


if __name__ == "__main__":

    main()
