#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../"))
from gfungi.config import *
from gfungi.common import check_paths, mkdir, check_path, get_version
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.0"

SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../stru_script/")
SOFTWARE_VERSION = {
    "fastp":{
        "GETVER": "%s/fastp --version 2>&1 |grep 'fastp'" % FASTP_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.20.0",
    },
    "STAR": {
        "GETVER": "%s/STAR 2>&1 |grep 'STAR version'" % STAR_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.7.9a",
    },
}


def create_star_index_task(genome, job_type, work_dir, threads=4):

    task = Task(
        id="star_index",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={star}:$PATH
STAR --runMode genomeGenerate --genomeDir {work_dir}/index \\
  --genomeFastaFiles {genome} --runThreadN {threads}
""".format(star=STAR_BIN,
            genome=genome,
            threads=threads,
            work_dir=work_dir
        )
    )

    return task, os.path.join(work_dir, "index")


def create_star_tasks(genome, reads1, reads2, index, work_dir,
                      threads=8, trim=5, qvalue=20):

    if len(prefixs) != len(reads1):
        prefixs = []
        n = 0
        for i in reads1:
            n += 1
            i = os.path.basename(i)
            i = i.split(".")[0]
            if i in prefixs:
                i = "%s_%s" % (i, n)
            prefixs.append(i)

    id = "star"
    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={fastp}:$PATH
fastp -i {{reads1}} -I {{reads2}} \\
-o {{prefixs}}.clean.r1.fq.gz -O {{prefixs}}.clean.r2.fq.gz \\
--thread {threads} -n 0 -f {trim} -F {trim} -t {trim} -T {trim} \\
-q {qvalue} --json fastp.json --html fastp.html
export PATH={star}:$PATH
STAR --runThreadN {threads} --genomeDir {index} \\
  -readFilesIn {{prefixs}}.clean.r1.fq.gz {{prefix}}.clean.r2.fq.gz\\
  --readFilesCommand zcat  --outWigType bedGraph  --outSAMtype BAM  SortedByCoordinate  --outSAMstrandField intronMotif

#export PATH={samtools}:$PATH
#samtools view --threads {threads} -Sb Aligned.sortedByCoord.out.bam | samtools sort --threads {threads} -m 4G -o {{prefix}}.sorted.bam
#rm {{prefix}}.clean.r*.fq.gz
""".format(fastp=FASTP_BIN,
           star=STAR_BIN,
           samtools=SAMTOOLS_BIN,
           index=index,
           threads=threads,
           trim=trim,
           qvalue=qvalue),
        reads1=reads1,
        reads2=reads2,
        prefixs=prefixs
    )

    join_task = Task(
        id="merge_star" % id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={samtools}:$PATH
samtools merge -f -c --threads {threads} Aligned.bam */Aligned.sortedByCoord.out.bam
samtools --threads {threads} -m 4G -o Aligned.sorted.bam Aligned.bam
export PATH={augustus}:$PATH
filterBam --uniq --in Aligned.sorted.bam  --out Aligned.clean.bam
bam2hints --intronsonly --in=Aligned.clean.bam --out=introns.gff
perl {script}/filterIntronsFindStrand.pl {genome} introns.gff --score > introns.f.gff
#rm */Aligned.sortedByCoord.out.bam Aligned.bam
""".format(samtools=SAMTOOLS_BIN,
           augustus=AUGUSTUS_BIN,
           script=SCRIPTS,
           genome=genome,
           )
    )

    join_task.set_upstream(*tasks)

    return tasks, join_task, os.path.join(work_dir, "Aligned.sorted.bam")


def run_star(genome, reads1, reads2, job_type, work_dir,
             threads=4, concurrent=10, refresh=15):

    genome = check_path(genome)
    reads1 = check_paths(reads1)
    reads2 = check_paths(reads2)
    work_dir = mkdir(os.path.abspath(work_dir))
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    options["software"]["STAR"] = {
        "version": get_version(SOFTWARE_VERSION["STAR"]),
        "option": "default"
    }
    options["software"]["fastp"] = {
        "version": get_version(SOFTWARE_VERSION["fastp"]),
        "option": "-n 0 -f 5 -F 5 -t 5 -T 5 -q 20"
    }

    dag = DAG("run_star")
    index_task, index = create_star_index_task(
        genome=genome,
        job_type=job_type,
        work_dir=work_dir,
        threads=threads
    )
    dag.add_task(index_task)

    star_tasks, join_task, bam = create_star_tasks(
        genome=genome,
        reads1=reads1,
        reads2=reads2,
        index=index,
        work_dir=work_dir,
        threads=threads,
        trim=5,
        qvalue=20
    )
    index_task.set_downstream(*star_tasks)
    dag.add_task(*star_tasks)
    dag.add_task(join_task)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return options


def star(args):

    options = run_star(
        genome=args.genome,
        reads1=args.reads1,
        reads2=args.reads2,
        threads=args.threads,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir
    )
    with open(os.path.join(args.out_dir, "star.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def add_hlep_args(parser):

    parser.add_argument("genome", metavar="STR", type=str,
        help="Input genome file(fata).")
    parser.add_argument("-r1", "--read1", metavar="FILE", nargs='+', type=str, required=True,
        help="Input the second-generation transcriptome data R1.")
    parser.add_argument("-r2", "--read2", metavar="FILE", nargs='+', type=str, required=True,
        help="Input the second-generation transcriptome data R2.")
    parser.add_argument("-t", "--threads", metavar="INT", type=int, default=4,
        help="Threads used to run stra (default: 4)")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("-w", "--work_dir", metavar="DIR", type=str, default="work",
        help="Work directory (default: current directory)")

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
    star.py Star comparison of transcriptome data.

attention:
    star.py genome.fasta --read1 *.R1.fq.gz --read2 *.R2.fq.gz

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    star(args)


if __name__ == "__main__":

    main()
