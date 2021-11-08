#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_fasta(fasta):
    '''Read fasta file'''

    seq = ''

    for line in open(fasta, 'r'):
        line = line.strip()

        if not line:
            continue
        if line.startswith('>'):
            if len(seq.split('\n'))==2:
                seq = seq.split('\n')
                yield seq[0], seq[1]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    if len(seq.split('\n'))==2:
        seq = seq.split('\n')
        yield [seq[0], seq[1]]


def stat_len(files):

    length = 0

    for line in read_fasta(files):
        length += len(line[1])

    return length


def read_vcf(files):

    for line in open(files, 'r'):
        line = line.strip()

        if not line or line.startswith('#'):
            continue
        yield line.split('\t')


def vcf_class(files):

    homo_snp_list = []
    heter_snp_list = []
    homo_indel_list = []
    heter_indel_list = []

    for line in read_vcf(files):
        if "COMPLEX" in line[7].upper():
            continue
        genotype = line[9].strip().split(':')
        if genotype[0]!='0/1' and genotype[0]!='1/1':
            continue

        type = genotype[0].strip().split('/')
        p = line[8].strip().split(':').index('AD')
        depth_list = genotype[p].strip().split(',')

        if genotype[0]=='0/1':
            depth = int(depth_list[0])+int(depth_list[1])
        else:
            depth = int(depth_list[1])

        if type[0]==type[1]:
            if "INDEL" in line[7].upper() or "INS" in line[7].upper():
                homo_indel_list.append(depth)
            else:
                homo_snp_list.append(depth)
        else:
            if "INDEL" in line[7].upper() or "INS" in line[7].upper():
                heter_indel_list.append(depth)
            else:
                heter_snp_list.append(depth)

    return homo_snp_list, heter_snp_list, homo_indel_list, heter_indel_list


def stat_number(vcf_list, depth):

    n = 0

    for i in vcf_list:
        if i>=depth:
            n += 1
        else:
            continue

    return n


def stat_vcf(files, genome, depth_list, out_file):

    out_dict = {}
    length = stat_len(genome)
    homo_snp_list, heter_snp_list, homo_indel_list, heter_indel_list = vcf_class(files)

    if len(depth_list)>=1:
        depth_list = depth_list.strip().split(',')
    else:
        print('Please enter the depth of coverage.')

    for i in depth_list:
        i = int(i)
        homo_snp_number = stat_number(homo_snp_list, i)
        heter_snp_number = stat_number(heter_snp_list, i)
        homo_indel_number = stat_number(homo_indel_list, i)
        heter_indel_number = stat_number(heter_indel_list, i)
        out_dict[i] = [heter_snp_number, heter_indel_number, homo_snp_number, homo_indel_number]

    output = open(out_file, 'w')
    output.write('#Depth\tHetero SNP\tHetero Indel\tHomo SNP\tError rate by Homo SNP(%)\tHomo Indel\tError rate by Homo Indel(%)\tError rate by homo variants(%)\tAccuracy genome(%)\n')

    for i in sorted(out_dict):
        line=out_dict[i]
        output.write('depth>={0}x\t{1:,}\t{2:,}\t{3:,}\t{4:6f}\t{5:,}\t{6:.6f}\t{7:.6f}\t{8:.6f}\n'.format(i, line[0], line[1], line[2], line[2]*100.0/length, line[3], line[3]*100.0/length, (line[2]+line[3])*100.0/length, 100-(line[2]+line[3])*100.0/length))
    output.close()


def stat_vcf_help(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Input the file in vcf format.')
    parser.add_argument('-g', '--genome', metavar='FILE', type=str, required=True,
        help='Input genome file.')
    parser.add_argument('-d', '--depth', metavar='STR', type=str, default='1,5,10',
        help='Input the statistical coverage depth,default=1,5,10.')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default='out.csv',
        help='The name of the output file.')

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
    stat_consistency.py  Sequence consistency statistics.

attention:
    stat_consistency.py -i var.vcf -g genome.fasta
    stat_consistency.py -i var.vcf  -g genome.fasta -d 1,5,10  -o stat_vcf.csv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = stat_vcf_help(parser).parse_args()

    stat_vcf(args.input, args.genome, args.depth, args.out)


if __name__ == "__main__":

    main()
