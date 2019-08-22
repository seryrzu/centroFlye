import argparse
from collections import defaultdict
import os
import statistics
from utils.bio import read_bio_seqs
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
    if len(conv) == 0:
        return [], [], None, None
    periods, bin_convs = [], []
    l, r = 0, 0
    best_l, best_r = 0, 0
    while r < len(conv):
        while r < len(conv) and conv[r] - conv[l] <= 2 * bin_size:
            r += 1
        period = int(statistics.median(conv[l:r]))
        periods.append(period)
        bin_convs.append(r - l)
        if r - l > best_r - best_l:
            best_l, best_r = l, r
        l += 1

    sorted_bc_p = sorted(zip(bin_convs, periods), reverse=True)
    bin_convs = [x[0] for x in sorted_bc_p]
    periods = [x[1] for x in sorted_bc_p]

    uniq_periods = set()
    filt_periods, filt_bin_convs = [], []
    for period, bin_conv in zip(periods, bin_convs):
        if period not in uniq_periods:
            uniq_periods.add(period)
            filt_periods.append(period)
            filt_bin_convs.append(bin_conv)
    return filt_periods, filt_bin_convs, conv[best_l], conv[best_r-1]


def run_on_read(seq, seq_id, k, bin_size, outdir):
    rep_kmers = get_repetitive_kmers(seq, k)
    conv, union_conv = get_convolution(rep_kmers)
    periods, bin_convs, bin_left, bin_right = \
        get_period_info(union_conv, bin_size=bin_size)
    plt.hist(union_conv, bins=100)
    plt.title(f'Tandem read convolution, {seq_id[:8]}, period={periods[0]}')
    plt.savefig(os.path.join(outdir, f'{seq_id[:8]}.pdf'), format='pdf')
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input reads", required=True)
    parser.add_argument("-o", "--outdir",
                        help="Output directory",
                        required=True)
    parser.add_argument("-k", help="kmer len", type=int, default=19)
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
