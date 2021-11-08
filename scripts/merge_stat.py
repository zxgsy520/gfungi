#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_csv(file):
    
    for line in open(file):
        line = line.strip()

        if not line or line.startswith('#'):
            continue

        yield line.split('\t')


def obtain_title(file):
    
    for line in open(file):
        line = line.strip()
        
        if not line:
            continue
        
        if line.startswith('#'):
            line = line.strip('#').split('\t')
            break
        else:
            line = []

    return line


def merge_stat(files, name):

    output = open('%s.tsv' % name, 'w')
    title = obtain_title(files[0])
    
    output.write('#sample\t%s\n' % ('\t'.join(title)))
    for i in files:
        sample = i.strip().split('/')[-1].split('.')[0]
        for line in read_csv(i):
            output.write('%s\t%s\n' % (sample, '\t'.join(line)))

    output.close()
        

def add_help_parser(parser):

    parser.add_argument('-i', '--input', metavar='FILE', nargs='+', type=str, required=True,
        help='Input file.')
    parser.add_argument('-n', '--name', metavar='STR', type=str, default='out',
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
    merge_stat.py Combine statistical results from different samples

attention:
    merge_stat.py -i *.stat_map.csv 
    merge_stat.py -i *.stat_map.csv -n name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_parser(parser).parse_args()

    merge_stat(args.input, args.name)


if __name__ == "__main__":

    main()
