#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import logging
import argparse

from collections import OrderedDict

from gfungi.config import *
from gfungi.common import check_path, check_paths, mkdir, read_files, get_version
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.1.0"



MAP_VERSION = {
    "splitfp": {
        "GETVER": "%s/splitfp.py -h 2>&1|grep 'version'" % SCRIPTS,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.1.0"
    },
    "bwa": {
        "GETVER": "%s/bwa 2>&1|grep 'Version'" % BWA_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "0.7.17"
    },
    "samtools": {
        "GETVER": "%s/samtools 2>&1|grep -i '^Version:'" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.9"
    },
}

def split_ngs(read1, read2, prefix, number, job_type, work_dir, out_dir="./"):

    option = {}
    option["splitfp"] = {
        "version": get_version(MAP_VERSION["splitfp"]),
        "option": "-r1 -r2"
    }

    dag = DAG("split_ngs")
    task = Task(
        id="split_ngs",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
{script}/splitfp.py -r1 {read1} -r2 {read2} -o {prefix} -n {number}
""".format(script=SCRIPTS,
            read1=read1,
            read2=read2,
            prefix=prefix,
            number=number
        )
    )

    dag.add_task(task)
    do_dag(dag, 8, 10)

    temp = read_files(work_dir, "%s.r1.part_*.fastq" % prefix)
    reads1 = []
    reads2 = []

    for i in temp:
        reads1.append(i)
        reads2.append(i.replace('.r1.part_','.r2.part_'))

    return option, reads1, reads2


def bwa_index_task(genome, prefix, job_type, work_dir, out_dir="./"):

    task = Task(
        id="bwa_index",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={bwa}:$PATH
ln -sf {genome} {prefix}.fasta
bwa index {prefix}.fasta
""".format(
            bwa=BWA_BIN,
            genome=genome,
            prefix=prefix
        )
    )

    return task, os.path.join( work_dir, "%s.fasta" % prefix)


def bwa_mem_tasks(reads1, reads2, genome, thread, job_type,
                  work_dir, out_dir="./"):

    option = {}
    option["bwa"] = {
        "version": get_version(MAP_VERSION["bwa"]),
        "option": "default"
    }
    option["samtools"] = {
        "version": get_version(MAP_VERSION["samtools"]),
        "option": "default"
    }
    number = list(range(len(reads1)))

    id = "bam"
    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={bwa}:{samtools}:$PATH
bwa mem -t {thread} {genome} {{reads1}} {{reads2}}\\
  |samtools view --threads {thread} -bS\\
  |samtools sort --threads {thread} -m 4G -o {{number}}.sort.bam
""".format(bwa=BWA_BIN,
           samtools=SAMTOOLS_BIN,
           genome=genome,
           thread=thread,
           ),
        reads1=reads1,
        reads2=reads2,
        number=number)

    return tasks, option, os.path.join(work_dir, "*/*.sort.bam")


def bwa_merge_bam_task(bams, prefix, thread, job_type, work_dir, out_dir="./"):

    task = Task(
        id="merge_bam",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={samtools}:$PATH
samtools merge -f -c --threads {thread} {prefix}.sorted.bam {bams}
samtools index {prefix}.sorted.bam
rm -rf {bams}
""".format(samtools=SAMTOOLS_BIN,
            bams=bams,
            prefix=prefix,
            thread=int(thread*1.5),
        )
    )

    return task, os.path.join(work_dir, "%s.sorted.bam" % prefix)


def run_bwa_mem(reads1, reads2, genome, prefix, thread, job_type,
                concurrent, refresh, work_dir, out_dir):

    dag = DAG("bwa_mem")
    index_task, genome = bwa_index_task(
        genome=genome,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(index_task)

    bwa_tasks, option, bams = bwa_mem_tasks(
        reads1=reads1,
        reads2=reads2,
        genome=genome,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(*bwa_tasks)

    merge_task, bam = bwa_merge_bam_task(
        bams=bams,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(merge_task)
    index_task.set_downstream(*bwa_tasks)
    merge_task.set_upstream(*bwa_tasks)
    do_dag(dag, concurrent, refresh)

    return option, bam


def bwa_mem(read1, read2, genome, prefix, number, thread, job_type,
            concurrent, refresh, work_dir, out_dir="./"):

    genome = check_path(genome)
    if isinstance(read1, list):
        read1 = " ".join(check_paths(read1))
        read2 = " ".join(check_paths(read2))
    else:
        read1 = check_path(read1)
        read2 = check_path(read2)

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    work_dict = {
        "split": "01_split",
        "bwa": "02_bwa",
    }
    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))

    option, reads1, reads2 =split_ngs(
        read1=read1,
        read2=read2,
        prefix=prefix,
        number=number,
        job_type=job_type,
        work_dir=work_dict["split"],
        out_dir=work_dict["split"],
    )
    options["software"].update(option)

    option, bam = run_bwa_mem(
        reads1=reads1,
        reads2=reads2,
        genome=genome,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dict["bwa"],
        out_dir=work_dict["bwa"]
    )
    options["software"].update(option)

    return bam, options
