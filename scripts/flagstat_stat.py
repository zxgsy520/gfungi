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


def read_map(files):
    '''Read statistics file'''

    for line in open(files, 'r'):
        line = line.strip()

        if not line or line.startswith('#'):
            continue
        yield line.split(' + 0 ')


def stat_map(files):
    '''Statistical map related indicators'''

    map_dict = {}

    for line in read_map(files):
        if 'in total' in line[1]:
            map_dict['Total Reads'] = int(line[0].strip())

        elif 'mapped (' in line[1]:
            map_dict['Map Reads'] = int(line[0].strip())

        elif 'paired in sequencing' in line[1]:
            map_dict['Paired Reads'] = int(line[0].strip())

        elif 'with itself and mate mapped' in line[1]:
            map_dict['Paired Map Reads'] = int(line[0].strip())

        elif 'properly paired' in line[1]:
            map_dict['Properly Paired Reads'] = int(line[0].strip())
        else:
            continue

    return map_dict


def out_map(files, out_files):
    '''Output file'''

    map_dict = stat_map(files)
    output = open(out_files, 'w')
    if map_dict['Paired Reads']==0:
        properly_map = 0
    else:
        properly_map = map_dict['Properly Paired Reads']*100.0/map_dict['Paired Reads']

    output.write('''#Total Reads\tMap Reads\tMap Rate(%)\tPaired Reads\tPaired Map Reads\tProperly Paired Reads\tProperly Map Rate(%)
{:,}\t{:,}\t{:.2f}\t{:,}\t{:,}\t{:,}\t{:.2f}
'''.format(map_dict['Total Reads'], map_dict['Map Reads'], map_dict['Map Reads']*100.0/map_dict['Total Reads'],
    map_dict['Paired Reads'], map_dict['Paired Map Reads'], map_dict['Properly Paired Reads'], properly_map))
    output.close()


def out_map_help(parser):

    parser.add_argument('-i', '--input', metavar='FILE', type=str, required=True,
        help='Input the map file of samtool statistics.')
    parser.add_argument('-o', '--out', metavar='STR', type=str, default='map.csv',
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
    map_stat.py  Statistical summary of map.

attention:
    map_stat.py -i map.log
    map_stat.py -i map.log -o map.csv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = out_map_help(parser).parse_args()

    out_map(args.input, args.out)


if __name__ == "__main__":

    main()
