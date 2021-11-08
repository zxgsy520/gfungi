#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import logging
import argparse

from gfungi.config import *
from gfungi.bwa_mem import bwa_mem
from gfungi.minimap import minimap
from gfungi.common import check_path, mkdir, get_version
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.1.0"


MAP_VERSION = {
    "freebayes": {
        "GETVER": "%s/freebayes --version 2>&1" % FREEBAYES_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "1.3.5"
    },
    "bcftools": {
        "GETVER": "%s/bcftools --version 2>&1|grep -i 'bcftools'" % BCFTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.14"
    },
    "samtools": {
        "GETVER": "%s/samtools 2>&1|grep -i '^Version:'" % SAMTOOLS_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.9"
    },
}


def run_freebayes_task(genome, bam, prefix, job_type, work_dir, out_dir):

    option = {}
    option["freebayes"] = {
        "version": get_version(MAP_VERSION["freebayes"]),
        "option": "--min-coverage 1 --max-coverage 300"
    }

    task = Task(
        id="freebayes",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={freebayes}:$PATH
freebayes -f {genome} {bam} --min-coverage 1 --max-coverage 300 \\
  --min-base-quality 1 > {prefix}.vcf
{script}/stat_snp_indel.py -i {prefix}.vcf -d 1,5,10 -g {genome} \\
  -o {prefix}.stat_snp_indel.tsv
cp {prefix}.stat_snp_indel.tsv {prefix}.vcf {out_dir}
""".format(freebayes=FREEBAYES_BIN,
            script=SCRIPTS,
            genome=genome,
            bam=bam,
            prefix=prefix,
            out_dir=out_dir
        )
    )

    return task, option, os.path.join(out_dir, "%s.stat_snp_indel.tsv" % prefix)


def run_mpileup_task(genome, bam, prefix, job_type, work_dir, out_dir):

    option = {}
    option["bcftools"] = {
        "version": get_version(MAP_VERSION["bcftools"]),
        "option": "call -v -m "
    }
    option["samtools"] = {
        "version": get_version(MAP_VERSION["samtools"]),
        "option": "mpileup --max-idepth 100 -t AD,ADF,ADR,DP,SP -u -g"
    }

    task = Task(
        id="mpileup",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={samtools}:{bcftools}:$PATH
samtools mpileup --max-idepth 100 -t AD,ADF,ADR,DP,SP -u -g -f {genome} \\
  {bam}|bcftools call -v -m > {prefix}.vcf
{script}/stat_snp_indel.py -i {prefix}.vcf -d 1,5,10 -g {genome} \\
  -o {prefix}.stat_snp_indel.tsv
cp {prefix}.stat_snp_indel.tsv {prefix}.vcf {out_dir}
""".format(samtools=SAMTOOLS_BIN,
            bcftools=BCFTOOLS_BIN,
            script=SCRIPTS,
            genome=genome,
            bam=bam,
            prefix=prefix,
            out_dir=out_dir
        )
    )

    return task, option, os.path.join(out_dir, "%s.stat_snp_indel.tsv" % prefix)


def run_flagstat_task(bam, prefix, job_type, work_dir, out_dir):

    task = Task(
        id="flagstat",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={samtools}:$PATH
samtools flagstat {bam} >{prefix}.map
{script}/flagstat_stat.py -i {prefix}.map -o {prefix}.map.tsv
cp {prefix}.map.tsv {out_dir}
""".format(samtools=SAMTOOLS_BIN,
            script=SCRIPTS,
            bam=bam,
            prefix=prefix,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.map.tsv" % prefix)


def stat_coverage_task(genome, bam, prefix, job_type, work_dir, out_dir):

    task = Task(
        id="stat_coverage",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1 %s" % QUEUE,
        script="""
export PATH={samtools}:$PATH
samtools depth -aa {bam} > {prefix}.depth
{script}/stat_coverage.py -i {prefix}.depth -d 1,5,10,20 -o {prefix}.coverage.tsv
{script}/stat_length_gc.py -d {prefix}.depth -g {genome} -n {prefix}
cp {prefix}.coverage.tsv {prefix}.length_gc.tsv {out_dir}
""".format(samtools=SAMTOOLS_BIN,
            script=SCRIPTS,
            genome=genome,
            bam=bam,
            prefix=prefix,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.coverage.tsv" % prefix)


def stat_gc_depth_task(genome, bam, prefix, job_type, work_dir, out_dir,
                       window=10000, interval=10000):

    task = Task(
        id="stat_coverage",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2",
        script="""
export PATH={samtools}:$PATH
samtools depth -aa {bam} > {prefix}.depth
{script}/stat_coverage.py -i {prefix}.depth -d 1,5,10,20 -o {prefix}.coverage.tsv
{script}/stat_gc_depth.py -d {prefix}.depth -g {genome} \\
   -b {interval} -w {window} -e 100 -n {prefix}
{script}/draw_depth_gc.py -gcd {prefix}.stat_gc_depth.tsv -n {prefix}
{script}/stat_length_gc.py -g {genome} -d {prefix}.depth -n {prefix}
cp {prefix}.coverage.tsv {prefix}.gc_depth.pdf {prefix}.length_gc.tsv {out_dir}
cp {prefix}.gc_depth.png {prefix}.map.tsv {prefix}.stat_gc_depth.tsv {out_dir}
""".format(samtools=SAMTOOLS_BIN,
            script=SCRIPTS,
            genome=genome,
            bam=bam,
            prefix=prefix,
            window=window,
            interval=interval,
            out_dir=out_dir
        )
    )

    return task, os.path.join(out_dir, "%s.stat_gc_depth.tsv" % prefix) 


def run_ngs_mapping(genome, read1, read2, prefix, thread, job_type,
                       concurrent, refresh, work_dir, out_dir):

    genome = check_path(genome)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    bam, options = bwa_mem(
        read1=read1,
        read2=read2,
        genome=genome,
        prefix=prefix,
        number=5000000,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=work_dir
    )

    dag = DAG("genome_mapping")

    mpileup_task, option, stat_vcf = run_mpileup_task(
        genome=genome,
        bam=bam,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(mpileup_task)
    options["software"].update(option)

    flagstat_task, stat_map = run_flagstat_task(
        bam=bam,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(flagstat_task)

    coverage_task, stat_coverage = stat_coverage_task(
        genome=genome,
        bam=bam,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(coverage_task)

    do_dag(dag, concurrent, refresh)

    return options, stat_vcf, stat_map, stat_coverage


def run_gc_depth(genome, reads, prefix, dtype, window, interval, thread,
                 job_type, concurrent, refresh, work_dir, out_dir):

    genome = check_path(genome)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    bam, options = minimap(
        reads=reads,
        genome=genome,
        prefix=prefix,
        dtype=dtype,
        number=80000,
        thread=thread,
        job_type=job_type,
        concurrent=concurrent,
        refresh=refresh,
        work_dir=work_dir,
        out_dir=work_dir
    )

    dag = DAG("gc_depth")
    gc_depth_task, gc_depth = stat_gc_depth_task(
        genome=genome,
        bam=bam,
        prefix=prefix,
        window=window,
        interval=interval,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    dag.add_task(gc_depth_task)

    flagstat_task, stat_map = run_flagstat_task(
        bam=bam,
        prefix=prefix,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    dag.add_task(flagstat_task)

    do_dag(dag, concurrent, refresh)

    return options, stat_map, gc_depth
