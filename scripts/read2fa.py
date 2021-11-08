#!/nextomics/Software/Base/miniconda3/bin/python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import pysam
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
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
        line = str(line).strip().strip("b'").strip("\n'")

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
        line = str(line).strip().strip("b'").strip("\n'")

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq[0], seq[1]
            seq = []
            seq.append(line.strip("@").split()[0])
            continue
        seq.append(line)

    if len(seq)==4:
        yield seq[0], seq[1]


def read_bam(file):

    fp = pysam.AlignmentFile(file, 'rb', check_sq=False)

    LOG.info('process %r' % file)
    for record in fp:
        yield record.qname, record.seq


def run_read_fa(files, number):

    n = 0

    if number!='all':
        number=int(number)

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
            if n==number:
                 break
            print('>%s\n%s' % (seqid, seq))
            n +=1
        if number!='all' and n>=number:
            break


def add_args(parser):

    parser.add_argument('-i', '--input', metavar='FILE', nargs='+', type=str, required=True,
        help='Input fastq file.')
    parser.add_argument('-n', '--number', default='all',
        help='Output reads number, default=all')
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

    read2fa.py -- Integrated database

attention:
    read2fa.py txt.fastq >txt.fasta
''')
    args = add_args(parser).parse_args()

    run_read_fa(args.input, args.number)


if __name__ == "__main__":
    main()
