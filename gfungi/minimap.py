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
__version__ = "v1.0.1"


MAPX = {
    "ont": "map-ont",
    "clr": "map-pb",
    "clr": "pb",
    "ccs": "splice:hq",
    "hifi": "splice:hq"
}

MAP_VERSION = {
    "splitfp": {
        "GETVER": "%s/splitfp.py -h 2>&1|grep 'version'" % SCRIPTS,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "2.1.0"
    },
    "minimap2": {
        "GETVER": "%s/minimap2 --version" % MINIMAP_BIN,
        "REGEXP": "\d+\.\d+\-r\d+",
        "MINVER": "2.11-r797"
    },
    "samtools": {
        "GETVER": "%s/samtools 2>&1|grep -i '^Version:'" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.9"
    },
}

def split_reads_task(reads, prefix, number, job_type, work_dir, out_dir="./"):

    option = {}
    option["splitfp"] = {
        "version": get_version(MAP_VERSION["splitfp"]),
        "option": "default"
    }

    dag = DAG("split_tgs")
    task = Task(
        id="split_tgs",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
{script}/splitfp.py -r1 {reads} --number {number} -o {prefix}
""".format(script=SCRIPTS,
            reads=reads,
            prefix=prefix,
            number=number,
        )
    )
    dag.add_task(task)
    do_dag(dag, 8, 10)

    reads = read_files(work_dir, '%s.part_*.fast*' % prefix)

    return option, reads


def minimap_align_tasks(reads, genome, dtype, thread,
                        job_type, work_dir, out_dir="./"):

    option = {}
    option["minimap2"] = {
        "version": get_version(MAP_VERSION["minimap2"]),
        "option": "--secondary=no %s" % MAPX[dtype]
    }
    option["samtools"] = {
        "version": get_version(MAP_VERSION["samtools"]),
        "option": "default"
    }

    id = "minimap"
    number = list(range(len(reads)))

    tasks = ParallelTask(
        id=id,
        work_dir="%s/{id}" % work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PATH={samtools}:{minimap}:$PATH
minimap2 --secondary=no -t {thread} -ax {x} {genome} {{reads}} | \\
  samtools view --threads {thread} -bS | samtools sort --threads {thread} -m 4G \\
  -o {{number}}.sort.bam
""".format(minimap=MINIMAP_BIN,
            samtools=SAMTOOLS_BIN,
            genome=genome,
            x=MAPX[dtype],
            thread=thread,
        ),
        reads=reads,
        number=number
    )

    return tasks, option, os.path.join(work_dir, "*/*.sort.bam")


def merge_bam_task(bams, prefix, thread, job_type, work_dir, out_dir="./"):

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


def run_minimap(reads, genome, dtype, prefix, thread, job_type,
                concurrent, refresh, work_dir, out_dir):

    dag = DAG("minimap")
    minimap_tasks, option, bams = minimap_align_tasks(
        reads=reads,
        genome=genome,
        dtype=dtype,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(*minimap_tasks)

    merge_task, bam = merge_bam_task(
        bams=bams,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(merge_task)
    merge_task.set_upstream(*minimap_tasks)
    do_dag(dag, concurrent, refresh)

    return option, bam


def minimap(reads, genome, prefix, dtype, number, thread, job_type,
            concurrent, refresh, work_dir, out_dir):

    genome = check_path(genome)
    if isinstance(reads, list):
        reads = " ".join(check_paths(reads))
    else:
        reads = check_path(reads)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    work_dict = {
        "split": "01_split",
        "minimap": "02_minimap",
    }
    for k, v in work_dict.items():
        work_dict[k] = mkdir(os.path.join(work_dir, v))


    option, reads = split_reads_task(
        reads=reads,
        prefix=prefix,
        number=number,
        job_type=job_type,
        work_dir=work_dict["split"],
        out_dir=work_dict["split"],
    )
    options["software"].update(option)

    option, bam = run_minimap(
        reads=reads,
        genome=genome,
        dtype=dtype,
        prefix=prefix,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dict["minimap"],
        out_dir=work_dict["minimap"]
    )
    options["software"].update(option)

    return bam, options
