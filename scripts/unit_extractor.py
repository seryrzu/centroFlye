# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
from collections import defaultdict
import math
import os
import statistics
import subprocess
from bisect import bisect_left, bisect_right

# from timeit import default_timer as timer

from utils.bio import read_bio_seqs, write_bio_seqs
from utils.os_utils import smart_makedirs

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def get_repetitive_kmers(seq, k):
    kmers = defaultdict(list)
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        kmers[kmer].append(i)

    kmers = {kmer: pos for kmer, pos in kmers.items() if len(pos) > 1}
    return kmers


def get_convolution(rep_kmers):
    conv, union_conv = {}, []
    for kmer in rep_kmers:
        pos = rep_kmers[kmer]
        conv[kmer] = sorted(y - x for x, y in zip(pos[:-1], pos[1:]))
        union_conv += conv[kmer]
    union_conv.sort()
    return conv, union_conv


def get_period_info(conv, bin_size):
    # assert conv == sorted(conv)
    if len(conv) == 0:
        return [], [], None, None
    periods2bin_convs, bin_convs2periods = {}, {}

    periods, bin_convs = [], []
    l, r = 0, 0
    best_l, best_r = 0, 0
    while r < len(conv):
        while r < len(conv) and conv[r] - conv[l] <= 2 * bin_size:
            r += 1

        # period = statistics.median(conv[l:r])
        mid = l + (r-l) // 2
        if (r-l) % 2 == 0:
            period = (conv[mid] + conv[mid-1]) // 2
        else:
            period = conv[mid]

        # print(period, period_fast)
        # assert period == period_fast
        curr_bin_conv = r - l
        if period not in periods2bin_convs or \
                curr_bin_conv > periods2bin_convs[period]:
            bin_convs2periods[curr_bin_conv] = period
            if period in periods2bin_convs and \
                    curr_bin_conv > periods2bin_convs[period]:
                bin_convs2periods.pop(periods2bin_convs[period], None)
            periods2bin_convs[period] = curr_bin_conv
        if curr_bin_conv > best_r - best_l:
            best_l, best_r = l, r
        l += 1

    bin_convs, periods = zip(*sorted(bin_convs2periods.items(), reverse=True))
    return periods, bin_convs, conv[best_l], conv[best_r-1]


def get_hook_kmer(conv, bin_left, bin_right):
    hook_kmer, max_tandem_index = None, 0
    for kmer, dist in conv.items():
        tandem_index = \
            bisect_right(dist, bin_right) - bisect_left(dist, bin_left)
        if tandem_index > max_tandem_index:
            hook_kmer = kmer
            max_tandem_index = tandem_index
    return hook_kmer


def split_by_hook(seq, hook):
    hook_pos = []
    for i in range(len(seq) - len(hook) + 1):
        kmer = seq[i:i+len(hook)]
        if kmer == hook:
            hook_pos.append(i)
    splits = {}
    for start, end in zip(hook_pos[:-1], hook_pos[1:]):
        split = seq[start:end]
        split_id = f'split_{start}_{end}'
        splits[split_id] = split
    return splits


def run_on_read(seq, seq_id, k, bin_size, outdir):
    print("Getting repetitive kmers")
    rep_kmers = get_repetitive_kmers(seq, k)
    print("Getting union convolution")
    conv, union_conv = get_convolution(rep_kmers)
    print("Getting periods")
    periods, bin_convs, bin_left, bin_right = \
        get_period_info(union_conv, bin_size=bin_size)
    print(f"Selected period = {periods[0]}")
    print("Getting hook")
    hook = get_hook_kmer(conv, bin_left, bin_right)
    if hook is None:
        return
    print("Splitting by hook")
    splits = split_by_hook(seq, hook)
    med_len = \
        statistics.median_high([len(x) for x in splits.values()])

    for r_id in sorted(splits.keys()):
        r_al = splits[r_id]
        if len(r_al) == med_len:
            median_read_unit = r_al
            template_read = r_id
            break
    read_outdir = os.path.join(outdir, seq_id[:8])
    smart_makedirs(read_outdir)
    splits_outfile = os.path.join(read_outdir, 'splits.fasta')
    median_read_unit_fn = os.path.join(read_outdir, 'median_read_unit.fasta')
    write_bio_seqs(splits_outfile, splits)
    write_bio_seqs(median_read_unit_fn,
                   {template_read: median_read_unit})

    print("Running Flye")
    cmd = ['flye',
           f'--nano-raw', splits_outfile,
           '--polish-target', median_read_unit_fn,
           '-i', 2,
           '-t', 50,
           '-o', read_outdir]
    cmd = [str(x) for x in cmd]
    subprocess.check_call(cmd)

    plt.hist(union_conv, bins=100)
    plt.title(f'Tandem read convolution, {seq_id[:8]}, period={periods[0]}')
    plt.savefig(os.path.join(read_outdir, f'{seq_id[:8]}.pdf'), format='pdf')
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input reads", required=True)
    parser.add_argument("-o", "--outdir",
                        help="Output directory",
                        required=True)
    parser.add_argument("-k", help="kmer len", type=int, default=15)
    parser.add_argument("-b", "--bin-size",
                        help="bin size",
                        type=int,
                        default=10)
    params = parser.parse_args()
    smart_makedirs(params.outdir)

    reads = read_bio_seqs(params.input)
    for r_id, seq in reads.items():
        run_on_read(seq, seq_id=r_id, k=params.k, bin_size=params.bin_size,
                    outdir=params.outdir)


if __name__ == "__main__":
    main()
