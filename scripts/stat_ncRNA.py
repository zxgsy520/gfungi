#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import logging

LOG = logging.getLogger(__name__)

__version__ = "20191223"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def stat_genome_length(file):
    total = 0
    for line in open(file, 'r'):
        if line.startswith(">"):
            continue
        total += len(line.strip())
    return total


def stat_gff(args):
    genome_length = stat_genome_length(args.genome)
    stat = {'ncRNA': [0, 0], 'regulatory': [0, 0], 'rRNA': [0, 0], 'tRNA': [0, 0]}
    rrna = {'18S': [0, 0], "28S": [0, 0], '5.8S': [0, 0], '5S': [0, 0]}
    srna = {'snRNA': [0, 0], 'miRNA': [0, 0], 'spliceosomal': [0, 0], 'other': [0, 0]}
    for line in open(args.gff, 'r'):
        if line.startswith('#'):
            continue
        sp1 = line.strip().split()
        sp2 = line.strip().split('product=')
        stat[sp1[2]][0] += 1
        stat[sp1[2]][1] += int(sp1[4]) - int(sp1[3]) + 1
        if sp1[2] == 'rRNA':
            rrna[sp2[-1].split()[0]][0] += 1
            rrna[sp2[-1].split()[0]][1] += int(sp1[4]) - int(sp1[3]) + 1
        if sp1[2] == 'ncRNA':
            if 'microRNA' in sp2[-1]:
                srna['miRNA'][0] += 1
                srna['miRNA'][1] += int(sp1[4]) - int(sp1[3]) + 1
            elif 'Small nucleolar RNA' in sp2[-1]:
                srna['snRNA'][0] += 1
                srna['snRNA'][1] += int(sp1[4]) - int(sp1[3]) + 1
            elif 'spliceosomal' in sp2[1]:
                srna['spliceosomal'][0] += 1
                srna['spliceosomal'][1] += int(sp1[4]) - int(sp1[3]) + 1
            else:
                srna['other'][0] += 1
                srna['other'][1] += int(sp1[4]) - int(sp1[3]) + 1

    f_out = open('%s.ncRNA.stat' % args.p, 'w')
    f_out.write('#Type\tNumber\tAverage_length(bp)\tTotal_length(bp)\tPercentage(%)\n')
    for k, v in stat.items():
        if int(v[0]) > 0:
            f_out.write('{rnatype}\t{number:,}\t{ave_len:,.2f}\t{total_len:,}\t{percent:.4f}\n'.format(
                rnatype=k, number=int(v[0]), ave_len=float(v[1])/float(v[0]), total_len=int(v[1]), percent=100 * float(v[1])/float(genome_length)
            ))
        else:
            f_out.write('%s\t0\t0\t0\t0\n' % k)
    f_out.write('\n#rRNA_stat:\n')
    for k, v in rrna.items():
        if int(v[0]) > 0:
            f_out.write('{rnatype}\t{number:,}\t{ave_len:,.2f}\t{total_len:,}\t{percent:.4f}\n'.format(
                rnatype=k, number=int(v[0]), ave_len=float(v[1])/float(v[0]), total_len=int(v[1]), percent=100 * float(v[1])/float(genome_length)
            ))
        else:
            f_out.write('%s\t0\t0\t0\t0\n' % k)
    f_out.write('\n#ncRNA_stat:\n')
    for k, v in srna.items():
        if int(v[0]) > 0:
            f_out.write('{rnatype}\t{number:,}\t{ave_len:,.2f}\t{total_len:,}\t{percent:.4f}\n'.format(
                rnatype=k, number=int(v[0]), ave_len=float(v[1])/float(v[0]), total_len=int(v[1]), percent=100 * float(v[1])/float(genome_length)
            ))
        else:
            f_out.write('%s\t0\t0\t0\t0\n' % k)
    f_out.close()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
        stat ncRNA result
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-gff', required=True,
                        help='gff, such as *ncRNA.gff3')
    parser.add_argument('-genome',required=True,
                        help='genome fasta')
    parser.add_argument('-p', default='result',
                        help='prefix of output')
    args = parser.parse_args()

    stat_gff(args)


if __name__ == "__main__":
    main()
