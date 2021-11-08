#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import sys
import json
import shutil
import logging
import argparse
import os.path

from gfungi.config import *
from gfungi.busco_cegma import run_busco_cegma
from gfungi.genome_mapping import run_ngs_mapping, run_gc_depth

from gfungi.parser import add_assemble_qc_args
from thirdparty.dagflow import DAG, Task, do_dag
from gfungi.common import check_path, check_paths, mkdir, read_tsv, cat, get_version

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.1.0"


ASMQC_VERSION = {
    "quast": {
        "GETVER": "%s/quast -v 2>&1" % QUAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "5.0.2"
    },
}


def stat_genome_task(genome, prefix, thread, job_type, work_dir, out_dir,
                     concurrent=2, refresh=10):
    option = {}
    option["splitfp"] = {
        "quast": get_version(ASMQC_VERSION["quast"]),
        "option": "--min-contig 0 "
    }

    dag = DAG("stat_genome")
    task = Task(
        id="stat_genome",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
{script}/sort_genome.py {genome} >{prefix}.genome.fasta
{script}/stat_genome.py -s {prefix}.genome.fasta -r {prefix}.genome.tsv
{script}/plot_genome.py -g {prefix}.genome.fasta -n {prefix}

export PATH={quast}:$PATH
quast --threads {thread} --min-contig 0 {prefix}.genome.fasta --output-dir quast
cp {prefix}.genome.fasta {prefix}.genome.pdf {prefix}.genome.png {prefix}.genome.tsv {out_dir}
cp quast/report.html {out_dir}/{prefix}.report.html

#python {script}/html2png.py --url quast/report.html --paths {phantomjs} -o {prefix}.report.png
#cp {prefix}.report.png {out_dir}
#rm -rf quast
""".format(quast=QUAST_BIN,
            phantomjs=PHANTOMJS,
            script=SCRIPTS,
            genome=genome,
            prefix=prefix,
            thread=thread,
            out_dir=out_dir
        )
    )
    dag.add_task(task)
    do_dag(dag, concurrent, refresh)

    return option, os.path.join(work_dir, "%s.genome.fasta" % prefix)


def run_assemble_qc(genome, prefix, read1, read2, reads, mode, busco_database,
                    thread, job_type, concurrent, refresh, work_dir, out_dir,
                    dtype="ont", window=10000, interval=10000):

    genome = check_paths(genome)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    work_dict = {
        "genome": "00_Genome",
        "complete": "01_Complete",
        "ngs_mapping": "02_NGS_MAPPING",
        "tgs_mapping": "03_TGS_MAPPING"
    }

    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    option, genome = stat_genome_task(
        genome=genome,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["genome"]),
        out_dir=os.path.join(out_dir, work_dict["genome"])
    )

    options = run_busco_cegma(
        genome=genome,
        prefix=prefix,
        mode=mode,
        busco_database=busco_database,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["complete"]),
        out_dir=os.path.join(out_dir, work_dict["complete"])
    )
    options["software"].update(option)

    noptions, stat_vcf, stat_map, stat_coverage = run_ngs_mapping(
        genome=genome,
        read1=read1,
        read2=read2,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["ngs_mapping"]),
        out_dir=os.path.join(out_dir, work_dict["ngs_mapping"])
    )
    options["software"].update(noptions["software"])

    noptions, map_stat, gc_depth = run_gc_depth(
        genome=genome,
        reads=reads,
        prefix=prefix,
        dtype=dtype,
        window=window,
        interval=interval,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=os.path.join(work_dir, work_dict["tgs_mapping"]),
        out_dir=os.path.join(out_dir, work_dict["tgs_mapping"])
    )
    options["software"].update(noptions["software"])

    return options, gc_depth


def assemble_qc(args):

    options, gc_depth = run_assemble_qc(
        genome=args.genome,
        prefix=args.prefix,
        reads=args.reads,
        dtype=args.dtype,
        read1=args.read1,
        read2=args.read2,
        mode=args.mode,
        busco_database=args.busco_database,
        window=args.window,
        interval=args.interval,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )
    with open(os.path.join(args.out_dir, "assemble_qc.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
assemble_qc:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_assemble_qc_args(parser)
    args = parser.parse_args()
    assemble_qc(args)


if __name__ == "__main__":
    main()
