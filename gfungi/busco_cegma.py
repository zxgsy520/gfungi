#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
import argparse

from gfungi.config import *
from gfungi.common import check_path, mkdir, get_version
from thirdparty.dagflow import DAG, Task, do_dag

LOG = logging.getLogger(__name__)
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__version__ = "v1.1.0"


AUGUSTUS_CONFIG = "/Work/pipeline/software/Base/augustus/lastest/config"
AUGUSTUS_BIN = "/Work/pipeline/software/Base/augustus/lastest/bin/"
BUSCO_CONFIG_FILE = "/Work/pipeline/software/Base/busco/lastest/bin/config/myconfig.ini"
BUSCO_BIN = "/Work/pipeline/software/Base/busco/lastest/bin/"
PERL_LIB = "/Work/pipeline/software/Base/cegma/lastest/lib"
CEGMA_BIN = "/Work/pipeline/software/Base/cegma/lastest/bin/"
RMBLAST_BIN = "/Work/pipeline/software/Base/rmblast/lastest/bin/"
GENEWISE_BIN = "/Work/pipeline/software/Base/genewise/lastest/bin/"
WISECONFIGDIR = "/Work/pipeline/software/Base/genewise/lastest/wisecfg"
GENEID_BIN = "/Work/pipeline/software/Base/geneid/lastest/bin/"
HMMER_BIN = "/Work/pipeline/software/Base/hmmer/lastest/bin/"

SOFTWARE_VERSION = {
    "busco": {
        "GETVER": "%s/busco --version 2>&1" % BUSCO_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "5.2.2"
    },
    "augustus": {
        "GETVER": "%s/augustus 2>&1|grep AUGUSTUS |grep prediction" % BUSCO_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "3.4.0"
    },
    "cegma": {
        "GETVER": "%s/cegma -h 2>&1|grep 'cegma -'" % CEGMA_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "2.5"
    },
    "blastn": {
        "GETVER": "%s/blastn -version 2>&1|grep 'blastn:'" % RMBLAST_BIN,
        "REGEXP": "\d+\.\d+\.\d+\+",
        "MINVER": "2.11.0+"
    },
    "genewise": {
        "GETVER": "%s/genewise -help 2>&1|grep 'wise2'" % GENEWISE_BIN,
        "REGEXP": "\d+\-\d+\-\d+",
        "MINVER": "2-4-1"
    },
    "geneid": {
        "GETVER": "%s/geneid -v 2>&1|grep 'geneid v'" % GENEID_BIN,
        "REGEXP": "\d+\.\d+",
        "MINVER": "1.4"
    },
    "hmmer": {
        "GETVER": "%s/hmmscan -h 2>&1|grep 'HMMER'" % HMMER_BIN,
        "REGEXP": "\d+\.\d+\.\d+",
        "MINVER": "3.3.2"
    },
}

def run_cegma_task(genome, prefix, thread, job_type, work_dir, out_dir):

    option = {}
    option["cegma"] = {
        "version": get_version(SOFTWARE_VERSION["cegma"]),
        "option": "default"
    }
    option["blastn"] = {
        "version": get_version(SOFTWARE_VERSION["blastn"]),
        "option": "default"
    }
    option["genewise"] = {
        "version": get_version(SOFTWARE_VERSION["genewise"]),
        "option": "default"
    }
    option["geneid"] = {
        "version": get_version(SOFTWARE_VERSION["geneid"]),
        "option": "default"
    }
    option["hmmer"] = {
        "version": get_version(SOFTWARE_VERSION["hmmer"]),
        "option": "default"
    }

    task = Task(
        id="cegma",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export PERL5LIB="{perl_lib}"
export PATH={rmblast}:{genewise}:{geneid}:{hmmer}:{cegma}:$PATH
export CEGMA={cegma}/../
export CEGMATMP=/tmp
export WISECONFIGDIR={wise_configdir}
cegma -g {genome} -T {thread} -o {prefix}
cp {prefix}.completeness_report {out_dir}
rm {prefix}.cegma.*
""".format(perl_lib=PERL_LIB,
            rmblast=RMBLAST_BIN,
            genewise=GENEWISE_BIN,
            geneid=GENEID_BIN,
            hmmer=HMMER_BIN,
            cegma=CEGMA_BIN,
            work_dir=work_dir,
            blast=BLAST_BIN,
            wise_configdir=WISECONFIGDIR,
            genome=genome,
            prefix=prefix,
            thread=thread,
            out_dir=out_dir
        )
    )

    stat_task = Task(
        id="stat_cegma",
        work_dir=work_dir,
        type=job_type,
        script="""
{script}/plot_cegma.py -i {prefix}.completeness_report -o {prefix}
cp {prefix}.cegma.pdf {prefix}.cegma.png {prefix}.cegma.tsv {out_dir}
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            prefix=prefix,
            out_dir=out_dir
        )
    )

    stat_task.set_upstream(task)

    return task, stat_task, option


def run_busco_task(genome, prefix, mode, busco_database, thread,
                   job_type, work_dir, out_dir):

    dbname = busco_database.split("/")[-1]
    option = {}
    option["busco"] = {
        "version": get_version(SOFTWARE_VERSION["busco"]),
        "option": "--mode %s --lineage_dataset %s" % (mode, dbname)
    }
    option["augustus"] = {
        "version": get_version(SOFTWARE_VERSION["augustus"]),
        "option": "default"
    }

    task = Task(
        id="busco",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp %s %s" % (thread, QUEUE),
        script="""
export BUSCO_CONFIG_FILE="{busco_config}"
export AUGUSTUS_CONFIG_PATH="{augustus_config}"
export PATH={busco}:{augustus}:$PATH
busco --cpu {thread} --mode {mode} --force --lineage_dataset {busco_database} \\
--offline --in {genome} --out {prefix}.BUSCO --out_path {work_dir}
rm -rf busco_downloads
""".format(augustus=AUGUSTUS_BIN,
            busco_config=BUSCO_CONFIG_FILE,
            augustus_config=AUGUSTUS_CONFIG,
            busco=BUSCO_BIN,
            genome=genome,
            prefix=prefix,
            busco_database=busco_database,
            mode=mode,
            thread=thread,
            work_dir=work_dir
        )
    )

    stat_task = Task(
        id="stat_busco",
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 1",
        script="""
export PATH={python}:$PATH
python {script}/plot_busco.py -i {prefix}.BUSCO/short_summary.*.txt -o {prefix}
cp {prefix}.busco.pdf {prefix}.busco.png {prefix}.busco.tsv {out_dir}
rm -rf {prefix}.BUSCO
""".format(script=SCRIPTS,
            python=PYTHON_BIN,
            prefix=prefix,
            out_dir=out_dir
        )
    )

    stat_task.set_upstream(task)

    return task, stat_task, option


def run_busco_cegma(genome, prefix, mode, busco_database, thread, job_type,
                    concurrent, refresh, work_dir, out_dir):

    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)
    genome = check_path(genome)

    check_path("%s_odb10" % busco_database)
    options = {
        "software": OrderedDict(),
        "database": OrderedDict()
    }

    dag = DAG("run_busco_cegma")
    cegma_task, stat_cegma_task, option = run_cegma_task(
        genome=genome,
        prefix=prefix,
        thread=int(thread*1.5),
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(cegma_task)
    dag.add_task(stat_cegma_task)

    busco_task, stat_busco_task, option = run_busco_task(
        genome=genome,
        prefix=prefix,
        mode=mode,
        busco_database=busco_database,
        thread=thread,
        job_type=job_type,
        work_dir=work_dir,
        out_dir=out_dir
    )
    options["software"].update(option)
    dag.add_task(busco_task)
    dag.add_task(stat_busco_task)
    do_dag(dag, concurrent, refresh)

    return options


def busco_cegma(args):

    options = run_busco_cegma(
        genome=args.genome,
        prefix=args.prefix,
        mode=args.mode,
        busco_database=args.busco_database,
        thread=args.thread,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )
    with open(os.path.join(args.out_dir, "busco_cegma.json"), "w") as fh:
        json.dump(options, fh, indent=2)

    return 0


def add_busco_cegma_args(parser):

    parser.add_argument("genome", metavar="FLIE", type=str,
        help="Input genome file(fasta)")
    parser.add_argument("-n", "--prefix", metavar="STR", type=str, default="fungi",
        help="Sample name, default=fungi.")
    parser.add_argument("-m", "--mode", type=str, default="genome",
        choices=["genome", "transcriptome", "proteins"],
        help="Specify which BUSCO analysis mode to run, default=genome")
    parser.add_argument("-db", "--busco_database", metavar="STR", type=str,
        default="/Work/database/busco_db/odb10/fungi",
        help="Specify the name of the BUSCO lineage to be used. default=fungi")
    parser.add_argument("-t", "--thread", metavar="INT", type=int, default=4,
        help="Number of threads, default=8")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", default=".",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", default=".",
        help="Output directory (default: current directory)")

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
busco_cegma:

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_busco_cegma_args(parser)
    args = parser.parse_args()
    busco_cegma(args)


if __name__ == "__main__":
    main()
