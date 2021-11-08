#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re
import sys
import argparse
import logging


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def stat_blastn(args):

    seq_id = dict()

    f_out = open( args.p + '_alignment.info', 'w')
    f_blast = open('%s.uniq' % args.i, 'w')

    f_out.write('seq_ID\tsubject_ID\tquery_len\tsubject_len\tidentity\talign_len\tcoverage\tstitle\n')
    for line in open(args.i, 'r'):
        tem = line.strip().split('\t')
        if tem[0] not in seq_id.keys():
            seq_id[tem[0]] = 1

            f_blast.write(line)
            cov = round(100 * float(tem[3])/float(tem[12]), 2)
            f_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (tem[0], tem[1], tem[12], tem[13], tem[2], tem[3], cov))

            if re.search('chloroplast', tem[14].lower()):
                f_out.write('chloroplast\n')
            elif re.search('mitochondrion', tem[14].lower()):
                f_out.write('mitochondrion\n')
            elif re.search('ribosomal', tem[14].lower()):
                f_out.write('ribosomal\n')
            else:
                f_out.write('NA\n')

    f_out.close()
    f_blast.close()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''

For exmple:
        python stat_blast_out.py -i *.fasta.blastn.out -p prefix
        
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-i', required=True,
                        help='a file, *.fasta.blastn.out')
    parser.add_argument('-p', default='result',
                        help='prefix of output')
    args = parser.parse_args()

    stat_blastn(args)


if __name__ == "__main__":

    main()

