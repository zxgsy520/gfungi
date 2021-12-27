#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

from collections import OrderedDict

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("invicoun@foxmail.com",)
__email__ = "113178210@qq.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)

    fp.close()


def split_attr(attributes):

    r =  OrderedDict()
    for content in attributes.split(";"):
        content = content.strip()
        if not content:
            continue

        if ' "' not in content:
            print("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split(' "', 1)
        r[tag] = value.strip('"')

    return r


def stat_stringtie(files, fpkm=0.01):

    total_gene = set()
    express_gene = set()
    data = {}
    n = 0

    for file in files:
        n += 1
        sample = os.path.basename(file)
        sample = sample.split(".")[0]
        if sample in data:
            sample = "%s_%s" % (sample, n)
        data[sample] = 0

        for line in read_tsv(file, "\t"):
            if "transcript" not in line[2]:
                continue
            attr = split_attr(line[-1])
            total_gene.add(attr["gene_id"])

            if float(attr["FPKM"]) <= fpkm:
                continue
            express_gene.add(attr["gene_id"])
            data[sample] += 1
    LOG.info("The total number of genes is %s" % len(total_gene))
    print("#Sample\tAnnotated transcipts\tPercentage(%s)")
    for i in data:
        print("{0}\t{1:,}\t{2:.2f}".format(i, data[i], data[i]*100.0/len(total_gene)))
    print("Total\t{0:,}\t{1:.2f}".format(len(express_gene), len(express_gene)*100.0/len(total_gene)))

    return 0


def add_help_args(parser):

    parser.add_argument("input", nargs="+", metavar="FILE", type=str,
        help='Input file.')
    parser.add_argument("-f", '--fpkm', metavar="FLOAT", type=float, default=0,
        help='Input fpkm value, default=0.')

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
    stat_stringtie.py:Statistical genome expression

attention:
    stat_stringtie.py */*.transcripts.gtf >stringtie_stat.xls
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()

    stat_stringtie(args.input, args.fpkm)


if __name__ == "__main__":
    main()
