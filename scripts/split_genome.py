#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import sys
import logging


LOG = logging.getLogger(__name__)

__version__ = "0.1.0"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def read_fasta(fasta):
    seq = ''
    seq_id = ''
    b = 0
    with open(fasta,'r') as f_fa:
        for line in f_fa:
            if line.startswith('>'):
                if b == 1:
                    yield seq_id, seq
                    seq = ''
                seq_id = line.split()[0].replace('>', '')
                b = 1
            else:
                seq += line.strip()
        yield seq_id, seq


def sp_genome(genome, bin_size, prefix, work_dir):

    all_seq = []
    for seq_id, seq in read_fasta(genome):
        seq_length = len(seq)
        if seq_length > 1000000:
            seq_start = list(range(0, seq_length, 50000))
            seq_end = list(range(50000, seq_length, 50000))
            seq_end.append(seq_length)
            for i in range(0, len(seq_start)):
                split_seq = seq[seq_start[i]:seq_end[i]]
                all_seq.append('>%s_sub%s' % (seq_id, i+1))
                all_seq.append(split_seq)
        else:
            all_seq.append('>%s' % seq_id)
            all_seq.append(seq)

    num = 0
    file_num = 1
    f_out = open('%s/%s_part_%s.fasta' % (work_dir, prefix, file_num), 'w')
    for line in all_seq:
        if line.startswith('>'):
            num += 1
            if num >= bin_size:
                f_out.close()
                num = 0
                file_num += 1
                f_out = open('%s/%s_part_%s.fasta' % (work_dir, prefix, file_num), 'w')
                f_out.write('%s\n' % line)
            else:
                f_out.write('%s\n' % line)
        else:
            f_out.write('%s\n' % line)

    f_out.close()


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''split genome

For exmple:
        python split_genome.py -g genome.fa -bin 100
        
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-g', required=True,
                        help='genome file, .fasta format')
    parser.add_argument('-bin', default=100, type=int,
                        help='contig number of split genome file,default=100')
    parser.add_argument('-p', default='result',
                        help='prefix of output')
    parser.add_argument('-w', default='./',
                        help='work dir')
    args = parser.parse_args()

    sp_genome(args.g, args.bin, args.p, args.w)


if __name__ == "__main__":

    main()

