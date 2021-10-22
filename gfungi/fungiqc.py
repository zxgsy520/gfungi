#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse

from gfungi.config import *
from gfungi.common import check_paths, mkdir, check_path, get_version
from thirdparty.dagflow import DAG, Task, ParallelTask, do_dag
from gfungi.parser import add_fungi_qc_args

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.0.0"


def filter_tgs_task(reads, name, minlen, xmax, bins, job_type, work_dir, out_dir):

    task = Task(
        id="tgs_qc",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={python}:$PATH
{script}/filter_tgs.py -i {reads} -n {name} --minlen {minlen} --xmax {xmax} --bins {bins}
cp {name}.reads_length.pdf {name}.reads_length.png {name}.reads_stat.tsv {out_dir}
""".format(python=PYTHON_BIN,
            script=SCRIPTS,
            reads=reads,
            name=name,
            minlen=minlen,
            xmax=xmax,
            bins=bins,
            out_dir=out_dir
        )
    )

    return task, os.path.join(work_dir, '%s.clean.fasta' % name)


def contamination_eval_task(reads, name, number, thread, job_type ,work_dir, out_dir):

    option = {}
    option["blastn"] = {
        "version": get_version(SOFTWARE_VERSION["blastn"]),
        "option": "-outfmt 6 -max_target_seqs 5"
    }

    task = Task(
        id="cont_eval",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s" % thread,
        script="""
export PATH={blast}:{python}:$PATH
python {scripts}/fq2fa.py {reads} -n {number} --cut >{name}.choose.fa
blastn -query {name}.choose.fa -db {nt} \
-outfmt "6 std staxid sskingdom staxids" -max_target_seqs 5 -num_threads {thread} \
-out {name}.m6
python {scripts}/obtain_taxonomy.py -i {name}.m6 -t {taxonomy} -n {name}
python {scripts}/stat_taxonomy.py -i {name}.species_annotation.txt -rn {number} -n {name}
python {scripts}/plot_stat_species.py -i {name}.stat_species.tsv -n {name}
cp {name}.top10_species.png {name}.top10_species.pdf {name}.top10_species.tsv {out_dir}
cp {name}.species_classify.tsv {name}.species_annotation.txt {out_dir}
""".format(scripts=SCRIPTS,
            blast=BLAST_BIN,
            nt=NT_DATABASE,
            python=PYTHON_BIN,
            taxonomy=TAXONOMY,
            name=name,
            reads=reads,
            number=number,
            thread=thread,
            out_dir=out_dir
        )
    )

    return task, option


def run_fungi_qc(reads, name, minlen, job_type, concurrent, refresh, work_dir, out_dir, thread=4):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    dag = DAG("fungi_qc")
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    qc_task, clean_reads = filter_tgs_task(
        reads=reads,
        name=name,
        minlen=minlen,
        xmax=80000,
        bins=300,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)

    eval_task, option = contamination_eval_task(
        reads=reads,
        name=name,
        number=10000,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir)
    options["software"].update(option)

    dag.add_task(qc_task)
    dag.add_task(eval_task)
    do_dag(dag, concurrent, refresh)

    return clean_reads, options


def fungi_qc(args):

    reads, options =run_fungi_qc(
        reads=' '.join(check_paths(args.reads)),
        name=args.name,
        minlen=args.minlen,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir)

    with open(os.path.join(args.out_dir, "qc.json"), "w") as fh:
        json.dump(options, fh, indent=2)


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_fungi_qc_args(parser)
    args = parser.parse_args()
    fungi_qc(args)


if __name__ == "__main__":
    main()
