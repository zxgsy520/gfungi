#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import shutil
import logging
import os.path
import argparse

from datetime import datetime
from jinja2 import Template
from docxtpl import DocxTemplate
try:
    from ConfigParser import RawConfigParser
except:
    from configparser import RawConfigParser

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.2.1"


def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


def check_paths(obj):
    """
    check the existence of paths
    :param obj:
    :return: abs paths
    """

    if isinstance(obj, list):
        r = []
        for path in obj:
            r.append(check_path(path))

        return r
    else:
        return check_path(obj)


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def bp2unit(string):

    string = str(string).strip().replace(',', '')
    string = float(string)

    if string >= 1e9:
        string = "%.2fGb" % (string/1e9)
    elif string >= 1e6:
        string = "%.2fMb" % (string/1e6)
    elif string >= 1e3:
        string = "%.2fKb" % (string/1e3)
    else:
        string = "%.2fbp" % string

    return string


def bp2mbgb(string, types="gb"):

    string = str(string).strip().replace(',', '')

    if types == "gb":
        base = 1e9
    elif types == "mb":
        base = 1e6
    elif types == "kb":
        base = 1e3
    else:
        base = 1.0

    string = float(string)/base

    return "{0:,.2f}".format(string)


def read_config(cfg):
    """
    read config fron ini
    :param cfg:
    :return:
    """
    cfg = check_path(cfg)

    r = {}
    #config = ConfigParser()
    config = RawConfigParser()
    #config.read(cfg, encoding='utf-8')
    config.read(cfg)

    for section in config.sections():
        r[section] = {}

        for option in config.options(section):
            value = config.get(section, option)

            if isinstance(value, bytes):
                value = value.decode("utf-8")
            r[section][option] = value.strip()

    return r


def read_tsv(file, sep=None):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        yield line.split(sep)


def read_table_data(file):

    r = []
    for line in read_tsv(file, "\t"):
        r.append(line)

    return r


def read_table_species(file):

    r = []
    for line in read_tsv(file, "\t"):
        r.append(line)

    return r


def read_table_assembly(file):

    r = []
    for line in read_tsv(file, "\t"):
        r.append(line[0:3])

    return r


def read_table_busco(file):

    r = []
    for line in read_tsv(file, "\t"):
        if line[0] == "Type":
            continue
        r.append(line)

    return r


def read_table_cegma(file):

    r = []

    for line in read_tsv(file, "\t"):
        if line[0] == "Type" or not line[0]:
            continue
        r.append(line)

    return r


def read_table_ngs_map(file):

    r = []

    for line in read_tsv(file, "\t"):
        r = line

    return r


def read_table_coverage(file):

    r = []
    coverage = 0

    for line in read_tsv(file, "\t"):
        if line[0] == "Depth" or not line[0]:
            continue
        if "Coverage" in line[0]:
            coverage = line[1]
            continue
        r.append(line)

    return r, coverage


def read_table_snp_indel(file):

    r = []
    for line in read_tsv(file, "\t"):
        if line[0] == "Depth" or not line[0]:
            continue
        r.append(line)

    return r


def read_table_tgs_map(file):

    r = []

    for line in read_tsv(file, "\t"):
        r = line[0:3]

    return r


def run_report(cfg, assemble_json, table_reads, table_species, table_assembly,
               table_busco, table_cegma, table_ngs_map, table_ngs_coverage,
               table_snp_indel, table_tgs_map, table_tgs_coverage, figures,
               others, tpl_html, out_dir):

    out_dir = mkdir(out_dir)
    now = datetime.now()
    config = read_config(cfg)
    table_reads = read_table_data(table_reads)
    table_species = read_table_species(table_species)
    table_assembly = read_table_assembly(table_assembly)
    table_busco= read_table_busco(table_busco)
    table_cegma = read_table_cegma(table_cegma)
    table_ngs_map = read_table_ngs_map(table_ngs_map)
    table_ngs_coverage, ngs_coverage = read_table_coverage(table_ngs_coverage)
    table_snp_indel = read_table_snp_indel(table_snp_indel)
    table_tgs_map = read_table_tgs_map(table_tgs_map)
    table_tgs_coverage, tgs_coverage = read_table_coverage(table_tgs_coverage)

    r = {
        "project": "",
        "id": "",
        "sample": "",
        "sequencer": "SequelII",
        "author": "",
        "reviewer": "",
        "year": now.year,
        "month": now.month,
        "day": now.day,
        "cell_num": 1,
        "dtype": "",
        "raw_bases": table_reads[0][1],
        "pass_bases": table_reads[1][1],
        "genome_size": bp2unit(table_assembly[6][1]),
        "genome_n50": bp2unit(table_assembly[0][1]),
        "busco_db": "",
        "busco_ratio": table_busco[0][2],
        "assembly_core": table_cegma[0][3],
        "cegma_ratio": table_cegma[0][4],
        "complete_core": table_cegma[0][1],
        "cegma_core_ratio": table_cegma[0][2],
        "ngs_coverage": ngs_coverage,
        "ngs_map_rate": table_ngs_map[-1],
        "homo_snps": table_snp_indel[1][3],
        "homo_snp_ratio": table_snp_indel[1][4],
        "homo_indels": table_snp_indel[1][5],
        "homo_indel_ratio": table_snp_indel[1][6],
        "single_base_accuracy": table_snp_indel[1][8],
        "tgs_coverage": tgs_coverage,
        "tgs_map_rate": table_tgs_coverage[0][2],
        "ngs_cov1": table_ngs_coverage[0][2],
        "table_reads": table_reads,
        "table_species": table_species,
        "table_assembly": table_assembly,
        "table_busco": table_busco,
        "table_cegma": table_cegma,
        "table_ngs_map": table_ngs_map,
        "table_ngs_coverage": table_ngs_coverage,
        "table_snp_indel": table_snp_indel,
        "table_tgs_map": table_tgs_map,
        "table_tgs_coverage": table_tgs_coverage,
        "gc_depth_description": "",
        "software": {},
        "database": {}
    }

    r.update(config["general"])

    with open(assemble_json) as fh:
        js = json.load(fh)
        for k in js:
            r[k].update(js[k])

#    html_report
    for i in ["images", "static"]:
        temp = os.path.join(out_dir, i)
        if os.path.exists(temp):
            shutil.rmtree(temp)
        shutil.copytree(os.path.join(tpl_html, i), temp)
    
    for i in figures:
        shutil.copy(i, os.path.join(out_dir, "images/"))

    temp = mkdir(os.path.join(out_dir, "tables/"))
    for i in others:
        shutil.copy(i, temp)

    for j in ["main.html", "index.html"]:
        temp = open(os.path.join(tpl_html, j)).read()
        if isinstance(temp, bytes):
            temp = temp.decode('utf-8')
        tpl = Template(temp)

        with open(os.path.join(out_dir, j), "w") as fh:
            line = tpl.render(r)
            if isinstance(line, bytes):
                line = line.decode('utf-8')
            try:
                fh.write(line)
            except:
                reload(sys)
                sys.setdefaultencoding("utf-8")
                fh.write(line)


    return r


def report(args):

    run_report(
        cfg=args.cfg,
        assemble_json=args.assemble,
        table_reads=args.data,
        table_species=args.species,
        table_assembly=args.stat_asm,
        table_busco=args.busco,
        table_cegma=args.cegma,
        table_ngs_map=args.ngs_map,
        table_ngs_coverage=args.ngs_coverage,
        table_snp_indel=args.snp_indel,
        table_tgs_map=args.tgs_map,
        table_tgs_coverage=args.tgs_coverage,
        figures=args.figures,
        others=args.others, 
        tpl_html=args.html,
        out_dir=args.out
    )

    return 0


def add_report_args(parser):

    parser.add_argument("cfg", help="Input configuration file.")
    parser.add_argument("-o", "--out", metavar='FILE', type=str, required=True,
        help="Output result path.")
    parser.add_argument("--others",  nargs='+', metavar='FILE', type=str,
        help="Input statistics table.")

    json_group = parser.add_argument_group(title="Json", )
    json_group.add_argument("--assemble", metavar='FILE', type=str, required=True,
        help="Input the jison file for assemble analysis: assemble.json.")

    table_group = parser.add_argument_group(title="Tables", )
    table_group.add_argument("--data", metavar='FILE', type=str, required=True,
        help="Input data to count the quality control results: reads_stat.tsv.")
    table_group.add_argument("--species", metavar='FILE', type=str, required=True,
        help="Input species distribution statistics.")
    table_group.add_argument("--stat_asm", metavar='FILE', type=str,  required=True,
        help="Input assembly statistics.")
    table_group.add_argument("--busco", metavar='FILE', type=str, required=True,
        help="Input busco evaluation results.")
    table_group.add_argument("--cegma", metavar='FILE', type=str, required=True,
        help="Input cegma evaluation results.")
    table_group.add_argument("--ngs_map", metavar='FILE', type=str,  required=True,
        help="Input the result of the second-generation data mapping.")
    table_group.add_argument("--ngs_coverage", metavar='FILE', type=str, required=True,
        help="Input the coverage statistics of the second-generation data.")
    table_group.add_argument("--snp_indel", metavar='FILE', type=str, required=True,
        help="Input snp statistics.")
    table_group.add_argument("--tgs_map", metavar='FILE', type=str, required=True,
        help="Input the mapping results of the three generations of data.")
    table_group.add_argument("--tgs_coverage", metavar='FILE', type=str, required=True,
        help="Input the coverage statistics of the three-generation data.")

    figure_group = parser.add_argument_group(title="Figure", )
    figure_group.add_argument("--figures", nargs='+', metavar='FILE', type=str,
        help="Input the path of all pictures.")

    template_group = parser.add_argument_group(title="Template", )
    template_group.add_argument("--html", required=True,
        help="Input the htlm template file.")

    return parser


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Configuration file
________________________________________________________
[general]
project=1个真菌近完成图测序
id=BYMG2021081601
sample=Pj
sequencer=SequelII
author=张兴国
reviewer=张兴国
cell_num=1
dtype=hifi
busco_db=fungi
________________________________________________________
version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_report_args(parser)
    args = parser.parse_args()

    report(args)


if __name__ == "__main__":
    main()
