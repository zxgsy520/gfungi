#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import logging


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def contig_tax(args):

    boundary = ['Archaea', 'Bacteria', 'Fungi', 'Viridiplantae', 'Alveolata', 'Metazoa', 'Viruses', 'Viroids', 'Unclassified', 'Other']
    cns_num = dict()
    cns_sp = dict()

    for line in open(args.i, 'r'):
        tem = line.strip().split('\t')
        if len(tem) != 4:
            pass
        else:
            cns_num[tem[0]] = tem[3]
            classify = tem[3].strip().split(';')
            for cls in classify:
                if cls.strip() in boundary:
                    cns_sp[tem[0]] = tem[2] + '\t' + cls.strip()

    f_out = open(args.p + '_contig_species.txt', 'w')
    #f_stat = open(args.p + '_stat_species.txt', 'w')
    #f_stat.write('#species\ttaxonomy\tnumber\n')
    stat_sp = dict()
    for k, v in cns_sp.items():
        if v in stat_sp:
            stat_sp[v] += 1
        else:
            stat_sp[v] = 1
        f_out.write('%s\t%s\t%s\n' % (k, v, cns_num[k]))

    # b = zip(stat_sp.values(), stat_sp.keys())
    # c = reversed(list(sorted(b)))
    #
    # for k, v in c:
    #     f_stat.write('%s\t%s\n' % (v, k))
    #
    # f_stat.close()
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

For exmple:
        python contig_taxonomy.py -i *.species_annotation.txt -p prefix
        
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-i', required=True,
                        help='a file, *.species_annotation.txt')
    parser.add_argument('-p', default='result',
                        help='prefix of output')
    args = parser.parse_args()

    contig_tax(args)


if __name__ == "__main__":

    main()

