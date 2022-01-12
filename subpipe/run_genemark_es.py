#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
import logging
import argparse

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from gfungi.config import *
from gfungi.common import check_paths, mkdir, check_path, get_version, read_files
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag
from thirdparty.seqkit.split import seq_split

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.0"


BRAKER_BIN = "/Work/pipeline/software/Base/braker2/v2.1.6/bin/"
GMES_BIN = "/Work/pipeline/software/meta/gmes/v4.68/"
GFFREAD_BIN = "/Work/pipeline/software/Base/gffread/v0.12.7/"


def split_genome(genome, prefix, job_type, work_dir, size="10m"):

    dag = DAG("split_genome")
    task = Task(
        id="split_genome",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
{bin}/biotool seqsplit {genome} --size {size} \\
  --prefix {prefix} --workdir {work_dir}
""".format(bin=BIN,
            genome=genome,
            prefix=prefix,
            size=size,
            work_dir=work_dir,
        )
    )
    dag.add_task(task)
    do_dag(dag, 8, 10)

    genomes = read_files(work_dir, "%s.part*.fa" % prefix)

    return genomes


def create_gmes_tasks(genomes, prefix, threads, job_type,
                      work_dir="", out_dir=""):

    id = "gmes"

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={braker}:$PATH
#cp {gmes}/gm_key ~/.gm_key
perl {gmes}/gmes_petap.pl  --ES --cores {threads} --sequence {{genomes}}
{gffread}/gffread genemark.gtf -o- >genemark.gff
rm -rf data info output run
""".format(braker=BRAKER_BIN,
            gmes=GMES_BIN,
            gffread=GFFREAD_BIN,
            threads=threads),
        genomes=genomes,
    )

    join_task = Task(
        id="merge_%s" % id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
cat {id}*/genemark.gff >{prefix}.genemark.gff
cp {prefix}.genemark.gff {out_dir}
""".format(id=id,
           prefix=prefix,
           out_dir=out_dir)
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task


def run_genemark_es(genomes, prefix, job_type, work_dir, out_dir,
                    threads=10, concurrent=10, refresh=30):

    genomes = check_paths(genomes)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    if len(genomes) <= 1:
        data = mkdir(os.path.join(work_dir, "00_genome"))
        genomes = split_genome(
            genome=genomes[0],
            prefix=prefix,
            job_type=job_type,
            work_dir=data,
            size="10m"
        )

    dag = DAG("run_genemark_es")
    tasks, join_task = create_gmes_tasks(
        genomes=genomes,
        prefix=prefix,
        threads=threads,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(*tasks)
    dag.add_task(join_task)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return 0


def add_hlep_args(parser):

    parser.add_argument("-g", "--genomes", nargs='+', metavar="STR", type=str,
        help="Input genome file(fasta).")
    parser.add_argument("--prefix", metavar="STR", type=str, default="ZXG",
        help="Input sample name.")
    parser.add_argument("--threads", metavar="INT", type=int, default=4,
        help="Threads used to run blastp (default: 4)")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", type=str, default="work",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", type=str, default="./",
        help="Output directory (default: current directory)")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
URL: https://github.com/zxgsy520/gfungi
name:
    run_genemark_es.py Genome De novo Prediction Module.

attention:
    run_genemark_es.py -p genome.fasta --prefix name
    run_genemark_es.py -p *genome.fasta --prefix name

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    run_genemark_es(genomes=args.genomes, prefix=args.prefix, threads=args.threads,
                    job_type=args.job_type, work_dir=args.work_dir, out_dir=args.out_dir,
                    concurrent=args.concurrent, refresh=args.refresh)


if __name__ == "__main__":

    main()
