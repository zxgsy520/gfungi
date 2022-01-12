#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import json
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

UNIVEC = "/Work/database/UniVec/UniVec_Core.fasta,/Work/database/UniVec/UniVec.fasta"
SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../stru_script/")
SOFTWARE_VERSION = {
    "stringtie":{
        "GETVER": "%s/stringtie -h 2>&1 |grep 'StringTie v'" % STRINGTIE_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "v2.2.0",
    },
}

PASA_CONFIG = """
# PASA admin settings
#emails sent to admin on job launch, success, and failure
PASA_ADMIN_EMAIL=invicoun@foxmail.com
#database to manage pasa jobs; required for daemon-based processing.
PASA_ADMIN_DB=PASA2_admin
DATABASE={work_dir}/pasa.sqlite
#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=85
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
"""


def create_stringtie_task(genome, bam, job_type, work_dir, threads=4):

    task = Task(
        id="stringtie",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={stringtie}:{gffread}:$PATH
stringtie -o stringtie.gtf -p {threads} {bam}
ls {work_dir}/stringtie.gtf > gtf.list
stringtie --merge -p 8 -o merged.gtf gtf.list
gffread -w rnaseq.transcript.fasta -g {genome} merged.gtf
""".format(stringtie=STRINGTIE_BIN,
            gffread=GFFREAD_BIN,
            bam=bam,
            genome=genome,
            threads=threads,
            work_dir=work_dir
        )
    )

    return task, os.path.join(work_dir, "rnaseq.transcript.fasta")


def create_seqclean_task(rna_ngs, rna_tgs, database, job_type, work_dir,
                          threads=4):
    if rna_tgs:
        temp = """
python {script}/rename_id.py {rna_tgs}  -p transngs >trans.rename_tgs.fasta
seqclean trans.rename_tgs.fasta  -c {threads} -v {database}
less trans.rename_tgs.fasta.clean | grep \'>\'|sed \'s/>//g\'|awk \'{{print $1}}\' >tgs.acc
""".format(script=SCRIPTS, rna_tgs=rna_tgs, threads=threads, database=UNIVEC)
        rnatgs = os.path.join(work_dir, "tgs.acc")
    else:
        temp = ""
        rnatgs = ""

    task = Task(
        id="seqclean",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={seqclean}:$PATH
python {script}/rename_id.py {rna_ngs} -p transngs >trans.rename.fasta
seqclean trans.rename.fasta  -c {threads} -v {database}
{temp}
""".format(seqclean=SEQCLEAN_BIN,
            database=UNIVEC,
            script=SCRIPTS,
            rna_ngs=rna_ngs,
            rna_tgs=rna_tgs,
            threads=threads,
            temp=temp
        )
    )
    rnaseq = os.path.join(work_dir, "trans.rename.fasta")
    clean_rnaseq = os.path.join(work_dir, "trans.rename.fasta.clean")

    return task, rnaseq, clean_rnaseq, rnatgs


def create_pasa_task(prefix, genome, rnaseq, clean_rnaseq, rnatgs, job_type, work_dir, out_dir,
                     threads=4):

    if rnatgs:
        temp = "-f %s" % rnatgs
    else:
        temp = ""
    fo = open(os.path.join(work_dir, "pasa.config"), 'w')
    fo.write(PASA_CONFIG.format(work_dir=work_dir))
    fo.close()
    task = Task(
        id="pasa",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (threads, QUEUE),
        script="""
export PATH={pasa}:{gmap}:$PATH
cut -f 1 {genome} >genome.fasta
{pasa}/Launch_PASA_pipeline.pl -c pasa.config \\
  -C -R -g genome.fasta -T\\
  -u {rnaseq} \\
  -t {clean_rnaseq} \\
  {temp} --CPU {threads} --ALIGNERS gmap
python {script}/rename_pasa_gtf.py pasa.sqlite.pasa_assemblies.gtf >pasa_assemblies.rename.gtf
cp pasa.sqlite.pasa_assemblies.gff3 {out_dir}/{prefix}.pasa.gff3
""".format(pasa=PASA_BIN,
            gmap=GMAP_BIN,
            script=SCRIPTS,
            prefix=prefix,
            genome=genome,
            rnaseq=rnaseq,
            clean_rnaseq=clean_rnaseq,
            temp=temp,
            threads=threads,
            out_dir=out_dir
        )
    )
    gff = os.path.join(work_dir, "pasa.sqlite.pasa_assemblies.gff3")
    gtf = os.path.join(work_dir, "pasa_assemblies.rename.gtf")

    return task, gff, gtf


def run_pasa(prefix, genome, bam, rna_tgs, job_type, work_dir, out_dir,
             threads=4, concurrent=10, refresh=15):

    genome = check_path(genome)
    bam = check_path(bam)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    if rna_tgs:
        rna_tgs = check_path(rna_tgs)

    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }
    options["software"]["stringtie"] = {
        "version": get_version(SOFTWARE_VERSION["stringtie"]),
        "option": "default"
    }
    options["database"]["UniVec"] = {
        "version": "latest",
    }

    dag = DAG("run_pasa")
    stringtie_task, rna_ngs = create_stringtie_task(
        genome=genome,
        bam=bam,
        job_type=job_type,
        work_dir=os.path.join(work_dir, "01_stringtie"),
        threads=threads
    )
    dag.add_task(stringtie_task)

    seqclean_task, rnaseq, clean_rnaseq, rnatgs = create_seqclean_task(
        rna_ngs=rna_ngs,
        rna_tgs=rna_tgs,
        database=UNIVEC,
        job_type=job_type,
        work_dir=os.path.join(work_dir, "02_seqclean"),
        threads=threads
    )
    dag.add_task(seqclean_task)
    seqclean_task.set_upstream(stringtie_task)

    pasa_task, gff, gtf = create_pasa_task(
        prefix=prefix,
        genome=genome,
        rnaseq=rnaseq,
        clean_rnaseq=clean_rnaseq,
        rnatgs=rnatgs,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir,
        threads=threads
    )
    dag.add_task(pasa_task)
    pasa_task.set_upstream(seqclean_task)

    do_dag(dag, concurrent_tasks=concurrent, refresh_time=refresh)

    return options, gff, gtf

def pasa(args):

    options, gff, gtf = run_pasa(
        prefix=args.prefix,
        genome=args.genome,
        bam=args.bam,
        rna_tgs=args.rna_tgs,
        job_type=args.job_type,
        work_dir=args.work_dir,
        out_dir=args.out_dir,
        threads=args.threads,
        concurrent=args.concurrent,
        refresh=args.refresh
    )

    with open(os.path.join(args.out_dir, "pasa.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def add_hlep_args(parser):

    parser.add_argument("genome", metavar="FILE", type=str,
        help="Input genome file(fata).")
    parser.add_argument("-p", "--prefix", metavar="FILE", type=str, default="D1",
        help="Input sample name, default=D1.")
    parser.add_argument("-b", "--bam", metavar="FILE", type=str, required=True,
        help="Input the second-generation transcriptome data and genome data map results(bam).")
    parser.add_argument("-rt", "--rna_tgs", metavar="FILE", type=str, default="",
        help="Input the transcript data processed by the three generations(fasta).")
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
    parser.add_argument("-o", "--out_dir", metavar="DIR", type=str, default="out",
        help="Out directory (default: current directory)")

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
    pasa.py: Using transcriptome data to predict genes.

attention:
    pasa.py genome.fasta --bam rna.bam 

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    pasa(args)


if __name__ == "__main__":

    main()
