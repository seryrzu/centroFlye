# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import re
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align.AlignInfo import SummaryInfo
from collections import defaultdict, Counter
import statistics
import os
from bisect import bisect_left, bisect_right
from subprocess import call

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans

from read import Read
from trf import TRF
from utils.various import take_closest
from utils.os_utils import smart_makedirs
from debruijn_graph import LongestPathInDeBruijn, DeBruijnGraph
from fitting_alignment import FittingAlignment
from hamming_distance import HammingDistance

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


class TandemReadAnalyser:
    def __init__(self, params, max_tandem_mult_index_diff=0.25):
        self.params = params
        self.inp_form = self.params.input.split('.')[-1]
        self.max_tandem_mult_index_diff = max_tandem_mult_index_diff

    def __get_repetitive_kmers(self, read):
        kmers = defaultdict(list)
        k = self.params.k
        for i in range(len(read)-k+1):
            kmer = read[i:i+k]
            kmers[kmer].append(i)

        kmers = {kmer: pos for kmer, pos in kmers.items() if len(pos) > 1}
        return kmers

    def __get_convolution(self, rep_kmers):
        conv, union_conv = {}, []
        for kmer in rep_kmers:
            pos = rep_kmers[kmer]
            conv[kmer] = sorted(y - x for x, y in zip(pos[:-1], pos[1:]))
            union_conv += conv[kmer]
        union_conv.sort()
        return conv, union_conv

    def __get_period_info(self, conv):
        if len(conv) == 0:
            return [], [], None, None
        periods, bin_convs = [], []
        l, r = 0, 0
        best_l, best_r = 0, 0
        while r < len(conv):
            while r < len(conv) and conv[r] - conv[l] <= 2 * self.params.bin_size:
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

    def __get_max_tandem_mult_index(self, periods):
        periods = sorted(list(periods))
        if len(periods) == 0:
            return None

        max_tandem_mult_index = 1
        for period in periods:
            curr_tandem_mult_index = 1
            used = set([period])
            while True:
                mult_period = (curr_tandem_mult_index + 1) * period
                i, closest = take_closest(periods, mult_period)
                if closest in used:
                    break
                used.add(closest)
                diff = abs(mult_period - closest) / max(mult_period, closest)
                if diff > self.max_tandem_mult_index_diff:
                    break
                curr_tandem_mult_index += 1
            if curr_tandem_mult_index > max_tandem_mult_index:
                max_tandem_mult_index = curr_tandem_mult_index
        return max_tandem_mult_index

    def __get_classification_info(self, rep_kmers,
                                  periods, bin_convs,
                                  max_tandem_mult_index,
                                  trf_tuples,
                                  n_repeats):
        classif_info = {}
        classif_info['NRepKmers'] = len(rep_kmers)
        classif_info['NPeriods'] = len(periods)
        classif_info['MaxTandemMultIndex'] = max_tandem_mult_index
        classif_info['NRepeats'] = n_repeats
        classif_info['TRF_Repeats'] = trf_tuples
        bin_convs += [None] * max(0, self.params.M - len(bin_convs))
        periods += [None] * max(0, self.params.M - len(periods))
        for i in range(self.params.M):
            classif_info['BinConv_%d' % (i+1)] = bin_convs[i]
            classif_info['Periods_%d' % (i+1)] = periods[i]
        return classif_info

    def __sort_kmers_tandem_index(self, conv, left, right):
        max_tandem_index = 0
        anchor_kmer = None
        ti_kmer = []
        for kmer, dist in conv.items():
            tandem_index = bisect_right(dist, right) - bisect_left(dist, left)
            # tandem_index = sum(x >= left and x < right for x in dist)
            ti_kmer.append((tandem_index, kmer))
        ti_kmer.sort(reverse=True)
        return ti_kmer

    def __get_strings_surrounding_kmer(self, read, pos):
        delta = self.params.delta
        return [read.seq[p-delta:p+delta] for p in pos]

    def __get_strings_between_pos(self, read, pos):
        return [read.seq[x:y] for x, y in zip(pos[:-1], pos[1:])]

    def __get_top_alignment_scores(self, last_column, buf_len, n=50):
        last_column = [(x, i) for i, x in enumerate(last_column)]
        last_column.sort(reverse=True)
        top = []
        i = 0
        while i < len(last_column) and len(top) < n:
            cand_score, cand_pos = last_column[i]
            far = True
            for pos, score in top:
                if abs(pos - cand_pos) < buf_len:
                    far=False
                    break
            if far:
                top.append((cand_pos, cand_score))
            i += 1
        top.sort()
        # TODO: 1D cluster the scores and select the corresponding substrings.
        # This is our Repeat.
        return top

    def __align(self, s1, s2):
        last_column, _, _, _ = FittingAlignment(s1, s2, sigma=3)
        top_sc_pos = self.__get_top_alignment_scores(last_column, len(s2))

        # Select good alignments
        top_sc = np.array([x[1] for x in top_sc_pos]).reshape(-1, 1)
        mixture = KMeans(n_clusters=2)
        mixture.fit(top_sc)
        pred_good_alignments = mixture.predict(top_sc)
        top_sc_pos = [x for i, x in enumerate(top_sc_pos) if pred_good_alignments[i]]

        # top positions of alignment + heuristic for starting position
        top_pos = [x[0] for x in top_sc_pos]
        top_pos.insert(0, max(0, top_pos[0] - (len(s2) + self.params.delta)))
        return top_pos

    def __get_multiple_alignment(self, seqs):
        seqs_fn = os.path.join(self.params.outdir, 'seqs.fasta')
        mult_al_fn = os.path.join(self.params.outdir, 'seqs_aligned.fasta')
        with open(os.path.join(self.params.outdir, 'seqs.fasta'), 'w') as f:
            for seq in seqs:
                print('>', file=f)
                print(seq, file=f)
        clustalomega_cline = ClustalOmegaCommandline(infile=seqs_fn,
                                                     outfile=mult_al_fn,
                                                     verbose=True, auto=True, force=True)
        clustalomega_cline()
        seqs_aligned = AlignIO.read(mult_al_fn, 'fasta')
        return seqs_aligned

    def __get_consensus(self, seqs_aligned):
        seqs_consensus = SummaryInfo(seqs_aligned).gap_consensus(require_multiple=1,threshold=0)
        seqs_consensus = str(seqs_consensus)
        seqs_consensus = seqs_consensus.replace('-','')
        seqs_consensus = seqs_consensus.replace('X','')
        return seqs_consensus

    def __align_and_get_consensus(self, read, longest_path):
        top_pos = self.__align(read.seq, longest_path)
        seqs = self.__get_strings_between_pos(read, top_pos)
        seqs_aligned = self.__get_multiple_alignment(seqs)
        seqs_consensus = self.__get_consensus(seqs_aligned)
        return seqs_consensus


    def __get_unit_and_positions(self, read, conv, bin_left, bin_right):
        # Use De Bruijn graph to construct a longer string from k-mers
        ti_kmer = self.__sort_kmers_tandem_index(conv, bin_left, bin_right)
        hook_kmers = [p[1] for p in ti_kmer if p[0] >= self.params.hook]
        gr = DeBruijnGraph.DeBruijnGraphFromKMers(hook_kmers)
        longest_path = LongestPathInDeBruijn(gr)

        unit = longest_path
        for i in range(self.params.n_consensus_iter):
            unit = self.__align_and_get_consensus(read, unit)
        return unit


    def __run_on_read(self, read):
        rep_kmers = self.__get_repetitive_kmers(read)
        conv, union_conv = self.__get_convolution(rep_kmers)
        print(union_conv)
        plt.hist(union_conv, bins=30)
        plt.title('Tandem read convolution')
        plt.savefig('test_tandem.pdf', format='pdf')
        periods, bin_convs, bin_left, bin_right = self.__get_period_info(union_conv)
        max_tandem_mult_index = self.__get_max_tandem_mult_index(periods)
        # unit = self.__get_unit_and_positions(read, conv, bin_left, bin_right)
        trf = TRF.run_trf_on_read(read, self.params.outdir)
        trf_tuples = trf.get_tuples_cons_begin_end()
        n_repeats = len(trf_tuples)
        classif_info = self.__get_classification_info(rep_kmers,
                                                      periods, bin_convs,
                                                      max_tandem_mult_index,
                                                      trf_tuples,
                                                      n_repeats)

        # if self.params.simulated:
        #     consensus = self.trf.df['Consensus'][0]
        #     print(len(consensus))
        #     print(consensus)
        #     consensus += consensus
        #     _, score, consensus_al, unit_al = FittingAlignment(consensus, unit, sigma=3)
        #     hd = HammingDistance(consensus_al, unit_al)
        #     len_cons, len_unit = len(consensus_al), len(unit_al)
        #     print("Len of TRF consensus %d. Len of my consensus %d." %(len_cons, len_unit))
        #     print("Edit Distance (after rotation) %d" % hd)
        #     print(consensus_al)
        #     print(unit_al)

        return classif_info

    def run(self):
        results = defaultdict(list)
        read_ids = []
        for raw_read in SeqIO.parse(self.params.input, self.inp_form):
            read_ids.append(raw_read.id)
            read = Read.FromBiopyRead(raw_read, self.params.simulated)
            read_results = self.__run_on_read(read)
            for key, val in read_results.items():
                results[key].append(val)

        out_dataframe = pd.DataFrame(results, index=read_ids)
        out_dataframe.to_csv(os.path.join(self.params.outdir,
                                          'classification_data.csv'), sep=',', na_rep='None')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input fastq", required=True)
    parser.add_argument("-o", "--outdir", help="Output directory", required=True)
    parser.add_argument("-k", help="kmer len", type=int, required=True)
    parser.add_argument("-b", "--bin-size", help="bin size", type=int, default=10)
    parser.add_argument("-M", type=int, default=10)
    parser.add_argument("-n-consensus-iter", type=int, default=1)
    parser.add_argument("-d", "--delta", help="len of flanking seq to anchor k-mer", type=int, default=50)
    parser.add_argument("--hook", help="Minimal number of positions for a hook k-mer in max bin convolution", type=int, default=1)
    parser.add_argument("--simulated", action='store_true')
    params = parser.parse_args()
    smart_makedirs(params.outdir)

    TandemReadAnalyser(params).run()


if __name__ == "__main__":
    main()
