#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import logging
import collections
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch


LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Liu huifang",)
__email__ = "liuhuifang@grandomics.com"
__all__ = []


def read_fasta(fasta):
    seq = ''
    seq_id = ''
    b = 0
    with open(fasta, 'r') as f_fa:
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


def get_tax_percent(array):
    sub_tax = dict()
    sub_title = dict()
    sub_sp = dict()
    for i in array:
        tem = i.split('\t')
        if tem[0] not in sub_tax:
            sub_tax[tem[0]] = 0
        sub_tax[tem[0]] += 1
        if tem[1] not in sub_title:
            sub_title[tem[1]] = 0
        sub_title[tem[1]] += 1
        if tem[2] not in sub_sp:
            sub_sp[tem[2]] = 0
        sub_sp[tem[2]] += 1
    info = ''
    for k, v in sorted(sub_tax.items(), key=lambda x: x[1], reverse=True):
        per = round(v/len(array), 2)
        info += '%s=%s;' % (k, per)
    info1 = ''
    for k, v in sorted(sub_title.items(), key=lambda x: x[1], reverse=True):
        per = round(v/len(array), 2)
        info1 += '%s=%s;' % (k, per)
    info2 = ''
    for k, v in sorted(sub_sp.items(), key=lambda x: x[1], reverse=True):
        per = round(v/len(array), 2)
        info2 += '%s=%s;' % (k, per)
    return sorted(sub_sp.items(), key=lambda x: x[1], reverse=True)[0][0], sorted(sub_tax.items(), key=lambda x: x[1], reverse=True)[0][0], \
           sorted(sub_title.items(), key=lambda x: x[1], reverse=True)[0][0], info, info1, info2


def get_order(order_align_len, contig):

    subdict = dict()
    for k in order_align_len.keys():
        if contig in k:
            subdict[k] = order_align_len[k]
    for k, v in sorted(subdict.items(), key=lambda x: x[1], reverse=True):
        yield k


def plot_stat_bar(file, prefix):
   
    labels = []
    contig_num = []
    for line in open(file, 'r'):
        if line.startswith('#') or line.startswith('Total'):
            pass
        else:
            tem = line.strip().split()
            labels.append(tem[0])
            contig_num.append(float(tem[4]))

    plt.switch_backend('agg')
    fig, ax = plt.subplots()
    colors = plt.get_cmap('tab10')(np.linspace(0, 1, len(labels)))
    ax.bar(labels, contig_num, color=colors, width=0.5)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylim(0, max(contig_num) + max(contig_num)/10)
    ax.set_ylabel("Contig Length ratio(%)")
    for i in range(len(labels)):
        ax.text(labels[i], contig_num[i] + max(contig_num)/20, str(contig_num[i]) + '%', horizontalalignment='center', verticalalignment='center')
    fig.tight_layout()
    plt.savefig('%s.pdf' % prefix)
    plt.savefig('%s.png' % prefix)


def plot_stat(file, max_sp, prefix, index):
    
    labels = []
    contig_num = []
    other_labels = []
    other_num = []

    for line in open(file, 'r'):
        if line.startswith('#') or line.startswith('Total'):
            pass
        elif line.startswith(max_sp) or line.startswith('Nohit'):
            tem = line.strip().split()
            labels.append(tem[0])
            contig_num.append(float(tem[index]))
        else:
            tem = line.strip().split()
            other_labels.append('{}  {:.2f}%'.format(tem[0], float(tem[index])))
            other_num.append(float(tem[index]))
    labels.append('Others')
    contig_num.append(sum(other_num))

    # make figure and assign axis objects
    plt.switch_backend('agg')
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fig.subplots_adjust(wspace=-0.6, right=0.4)

    # pie chart parameters
    ratios = contig_num
    labels = labels
    explode = []
    for i in ratios:
        explode.append(0)
    colors = plt.get_cmap('Set3')(np.linspace(0, 1, len(labels)))
    # rotate so that first wedge is split by the x-axis
    ax1.pie(ratios, autopct='%1.2f%%', labels=labels, explode=explode, radius=1.25, colors=colors)

    # bar chart parameters
    xpos = 0
    bottom = 0
    ratios = other_num
    width = 1
    colors = [plt.get_cmap('Pastel1')(np.linspace(0, 1, 9))[i] for i in range(len(other_labels))]

    for j in range(len(ratios)):
        height = ratios[j]
        ax2.bar(xpos, height, width, bottom=bottom, color=colors[j])
        #ypos = bottom + ax2.patches[j].get_height() / 2
        bottom += height
        #ax2.text(xpos, ypos, "%0.2f%%" % (other_num[j]), ha='center')

    ax2.legend(other_labels, bbox_to_anchor=(0.65, 0), loc=3)
    ax2.axis('off')
    ax2.set_xlim(- 2.5 * width, 2.5 * width)

    # use ConnectionPatch to draw lines between the two plots
    # get the wedge data
    if len(ax1.patches) >= 3:
        theta1, theta2 = ax1.patches[2].theta1, ax1.patches[2].theta2
        center, r = ax1.patches[2].center, ax1.patches[2].r
        bar_height = sum([item.get_height() for item in ax2.patches])

        # draw top connecting line
        x = r * np.cos(np.pi / 180 * theta2) + center[0]
        y = r * np.sin(np.pi / 180 * theta2) + center[1]

        con = ConnectionPatch(xyA=(- width / 2, bar_height), xyB=(x, y),
            coordsA="data", coordsB="data", axesA=ax2, axesB=ax1)
        con.set_color("gray")
        con.set_linewidth(0.8)
        con.set_linestyle((0, (5, 10)))
        ax2.add_artist(con)

        # draw bottom connecting line
        x = r * np.cos(np.pi / 180 * theta1) + center[0]
        y = r * np.sin(np.pi / 180 * theta1) + center[1]
        con = ConnectionPatch(xyA=(- width / 2, 0), xyB=(x, y), coordsA="data",
            coordsB="data", axesA=ax2, axesB=ax1)
        con.set_color("gray")
        con.set_linewidth(0.8)
        con.set_linestyle((0, (5, 10)))
        ax2.add_artist(con)

    plt.tight_layout()
    plt.savefig('%s.pdf' % prefix)
    plt.savefig('%s.png' % prefix)


def generate_fasta(genome, stat_sp, genome_size, kingdom, prefix):
   
    for k in genome_size.keys():
        if kingdom not in stat_sp:
            stat_sp[kingdom] = []
        stat_sp[kingdom].append(k)

    for seqid, seq in read_fasta(genome):
        for k, v in stat_sp.items():
            if k == 'Nohit' or k == kingdom:
                f_out = open('%s_decontaminated_genome.fasta' % prefix, 'a')
            else:
                k = k.replace('/', '_')
                f_out = open('%s_%s.fasta' % (prefix, k), 'a')
            if seqid in v:
                f_out.write('>%s\n%s\n' % (seqid, seq))


def merge_sub(args):
    align_info = dict()
    seq_info = dict()
    order_align_len = dict()
    contig_id = dict()

    for line in open(args.align, 'r'):
        tem = line.strip().split('\t')
        if line.startswith('seq_ID'):
            pass
        else:
            if '_sub' in tem[0]:
                contig_id[tem[0].split('_sub')[0]] = 0
            else:
                contig_id[tem[0]] = 0
            order_align_len[tem[0]] = int(tem[5])
            align_info[tem[0]] = '%s\t%s\t%s\t%s\t%s\t%s\t%s' % (tem[1], tem[2], tem[3], tem[4], tem[5], tem[6], tem[7])

    for line in open(args.subcontig, 'r'):
        tem = line.strip().split('\t')
        if tem[0] in align_info.keys():
            seq_info[tem[0]] = '%s\t%s\t%s\t%s\t%s' % (tem[0], tem[1], tem[2], align_info[tem[0]], tem[3])

    f_out = open(args.p + '_blast_out.summary', 'w')
    f_out.write('seq_ID\tspecies\ttaxonomy\tncbi_id\tquery_len\tsubject_len\tidentity\talign_len\tcoverage\tstitle\ttaxonomy\n')

    seq_info_uniq = collections.OrderedDict()
    for k in sorted(contig_id.keys()):
        for subk in get_order(order_align_len, k):
            if subk in seq_info.keys():
                f_out.write('%s\n' % seq_info[subk])
                tem = seq_info[subk].split('\t')
                if float(tem[6]) > args.identity and float(tem[8]) > args.coverage:
                    tem[0] = tem[0].split('_sub')[0]
                    if tem[0] not in seq_info_uniq.keys():
                        seq_info_uniq[tem[0]] = []
                    seq_info_uniq[tem[0]].append('%s\t%s\t%s' % (tem[2], tem[9], tem[1]))
    f_out.close()

    stat_sp = dict()
    f_out1 = open(args.p + '_blast_out.summary.merge', 'w')
    f_out1.write('seq_ID\tmax_species\tmax_taxonomy\tmax_stitle\ttaxonomy\tstitle\tspecies\n')
    f_com = open(args.p + '_contig_contamination.txt', 'w')
    f_com.write('seq_ID\tmax_species\tmax_taxonomy\tmax_stitle\ttaxonomy\tstitle\tspecies\n')
    for k in seq_info_uniq.keys():
        if len(seq_info_uniq[k]) == 1:
            tem = seq_info_uniq[k][0].split('\t')
            f_out1.write('%s\t%s\t%s\t%s\tNA\tNA\tNA\n' % (k, tem[2], tem[0], tem[1]))

            if tem[1] != 'NA' and tem[1] != 'ribosomal':
                if 'Mitochondrion/Chloroplast' not in stat_sp.keys():
                    stat_sp['Mitochondrion/Chloroplast'] = []
                stat_sp['Mitochondrion/Chloroplast'].append(k)
            else:
                if tem[0] not in stat_sp.keys():
                    stat_sp[tem[0]] = []
                stat_sp[tem[0]].append(k)

            if tem[0] != args.kingdom:
                f_com.write('%s\t%s\t%s\t%s\tNA\tNA\tNA\n' % (k, tem[2], tem[0], tem[1]))
            elif tem[1] != 'NA' and tem[1] != 'ribosomal':
                f_com.write('%s\t%s\t%s\t%s\tNA\tNA\tNA\n' % (k, tem[2], tem[0], tem[1]))
        else:
            max_sp, max_tax, max_title, info, info1, info2 = get_tax_percent(seq_info_uniq[k])
            f_out1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (k, max_sp, max_tax, max_title, info, info1, info2))

            if max_title != 'NA' and max_title != 'ribosomal':
                if 'Mitochondrion/Chloroplast' not in stat_sp.keys():
                    stat_sp['Mitochondrion/Chloroplast'] = []
                stat_sp['Mitochondrion/Chloroplast'].append(k)
            else:
                if max_tax not in stat_sp.keys():
                    stat_sp[max_tax] = []
                stat_sp[max_tax].append(k)

            if max_tax != args.kingdom:
                f_com.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (k, max_sp, max_tax, max_title, info, info1, info2))
            elif max_title != 'NA' and max_title != 'ribosomal':
                f_com.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (k, max_sp, max_tax, max_title, info, info1, info2))

    f_out1.close()
    f_com.close()

    genome_size = dict()
    genome_sum = 0
    for seqid, seq in read_fasta(args.genome):
        seq_len = len(seq)
        genome_size[seqid] = seq_len
        genome_sum += int(seq_len)
    contig_sum = len(genome_size.keys())

    f_sp = open(args.p + '_stat_species.txt', 'w')
    f_sp.write('#sample\tContig Number\tContig Number ratio(%)\tContig Length(bp)\tContig Length ratio(%)\n')
    for k, v in sorted(stat_sp.items(), key=lambda x: len(x[1]), reverse=True):
        sub_sum = sum([genome_size[i] for i in v])
        f_sp.write('{}\t{:,}\t{:,.2f}\t{:,}\t{:,.2f}\n'.format(k, len(v), 100 * float(len(v)) / contig_sum, sub_sum, 100 * float(sub_sum)/genome_sum))
        for i in v:
            del genome_size[i]
    f_sp.write('Nohit\t{:,}\t{:,.2f}\t{:,}\t{:,.2f}\n'.format(len(genome_size.keys()), 100 * float(len(genome_size.keys())) / contig_sum, sum(genome_size.values()), 100 * float(sum(genome_size.values())) / genome_sum))
    f_sp.write('Total\t{:,}\t100.00\t{:,}\t100.00\n'.format(contig_sum, genome_sum))
    f_sp.close()

    generate_fasta(args.genome, stat_sp, genome_size, args.kingdom, args.p)
    plot_stat(args.p + '_stat_species.txt', args.kingdom, args.p + '_Contig_Number', 2)
    plot_stat_bar(args.p + '_stat_species.txt', args.p + '_Contig_Length')


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    parser.add_argument('-align',required=True,
                        help='*_alignment.info')
    parser.add_argument('-subcontig',required=True,
                        help='*_contig_species.txt')
    parser.add_argument('-g', '--genome', required=True,
                        help='genome sequence, *.fasta')
    parser.add_argument('-kingdom', choices=["Viridiplantae", "Metazoa", "Fungi"], default='Viridiplantae', type=str,
                        help='Choose kingdom,default=Viridiplantae')
    parser.add_argument('-identity', default=70, type=float,
                        help='blast identity, default=70')
    parser.add_argument('-coverage', default=1, type=float,
                        help='blast coverage, default=1')
    parser.add_argument('-p', default='result',
                        help='prefix of output')
    args = parser.parse_args()

    merge_sub(args)


if __name__ == "__main__":
    main()
