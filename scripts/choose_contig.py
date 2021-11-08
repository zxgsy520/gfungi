#!/nextomics/Software/Base/miniconda3/bin/python3
# -*- coding: utf-8 -*-

import argparse

__version__ = "0.1.0"
__author__ = ("Dan Huang",)
__email__ = "2876129239@qq.com"
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
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield seq[0], seq[1]


def choose_contig(genome, shortest, longest, name):

    if genome.endswith(".fasta") or genome.endswith(".fa") or genome.endswith(".fasta.gz") or genome.endswith(".fa.gz"):
        fh = read_fasta(genome)
    else:
        print('%s file format err' % genome)

    output = open('%s.%s_%s.fasta' % (name, shortest, longest), 'w')
    for seqid, seq in fh:
        if shortest <= len(seq) and len(seq) <= longest :
            output.write('>{0}\n{1}\n'.format(seqid, seq))
    output.close()


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    choose_contig.py --genome --shortest --longest --name
''')

    parser.add_argument('--genome', metavar='FILE', type=str, required=True,
        help='The genome file')
    parser.add_argument('--shortest', type=int, default=0,
        help='Pick contigs longer than the length')
    parser.add_argument('--longest', type=int, default=500000,
        help='Pick contigs shorter than the length')
    parser.add_argument('--name', type=str, default='fungi',
        help='Output file name.')
    args = parser.parse_args()
    choose_contig(args.genome, args.shortest, args.longest, args.name)


if __name__ == "__main__":
    main()
