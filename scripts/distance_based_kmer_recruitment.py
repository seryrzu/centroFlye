# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
import itertools
import sys
from collections import defaultdict
from utils.os_utils import smart_makedirs
from ncrf_parser import NCRF_Report
from read_kmer_cloud import get_reads_kmer_clouds


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncrf', help='NCRF report on reads', required=True)
    parser.add_argument('--coverage', help='Average coverage of the dataset',
                        type=int, required=True)
    parser.add_argument('--min-coverage',
                        help='minCov threshold',
                        type=int,
                        default=4)
    parser.add_argument('--outdir', help='Output directory', required=True)
    parser.add_argument('-k', type=int, default=19)
    parser.add_argument('--min-nreads', type=int, default=0)
    parser.add_argument('--max-nreads', type=int, default=sys.maxsize)
    parser.add_argument('--min-distance', type=int, default=1)
    parser.add_argument('--max-distance', type=int, default=150)
    parser.add_argument('--bottom', type=float, default=0.9)
    parser.add_argument('--top', type=float, default=3.)
    parser.add_argument('--kmer-survival-rate', type=float, default=0.34)
    parser.add_argument('--max-nonuniq', type=int, default=3)
    parser.add_argument('--verbose', action='store_true', default=True)
    params = parser.parse_args()
    return params


def get_kmer_freqs_from_ncrf_report(reads_ncrf_report,
                                    k, verbose,
                                    max_nonuniq):
    non_unique_freqs = defaultdict(int)
    all_kmers = defaultdict(int)
    for i, (r_id, record) in enumerate(reads_ncrf_report.records.items()):
        if i % 100 == 0 and verbose:
            print(i + 1, len(reads_ncrf_report.records))
        r_al = record.r_al
        r_al = r_al.replace('-', '')

        read_freq = defaultdict(int)
        for i in range(len(r_al)-k+1):
            kmer = r_al[i:i+k]
            read_freq[kmer] += 1

        for kmer, freq in read_freq.items():
            if freq > 1:
                non_unique_freqs[kmer] += 1
            if non_unique_freqs[kmer] <= max_nonuniq:
                all_kmers[kmer] += 1
            else:
                if kmer in all_kmers:
                    del all_kmers[kmer]
    return all_kmers


def get_rare_kmers(reads_ncrf_report, k,
                   bottom, top, coverage, kmer_survival_rate, max_nonuniq,
                   verbose):
    all_kmers = get_kmer_freqs_from_ncrf_report(reads_ncrf_report,
                                                k=k,
                                                verbose=verbose,
                                                max_nonuniq=max_nonuniq)

    left = bottom*coverage*kmer_survival_rate
    right = top*coverage*kmer_survival_rate

    rare_kmers = [kmer for kmer, freq in all_kmers.items()
                  if left <= freq <= right]
    rare_kmers = set(rare_kmers)
    if verbose:
        print(f'# rare kmers: {len(rare_kmers)}')
    return rare_kmers


def get_kmer_dist_map(reads_kmer_clouds, kmers,
                      min_n, max_n,
                      min_d, max_d,
                      verbose):
    def index_clouds(reads_kmer_clouds, min_n, max_n):
        indexed_reads_kmer_clouds = {}
        for r_id, kmer_clouds in \
                itertools.islice(reads_kmer_clouds.items(), min_n, max_n):
            kmers = kmer_clouds.kmers
            indexed_read_kmer_clouds = []
            for cloud in kmers:
                indexed_cloud = [kmer_index[kmer] for kmer in cloud]
                indexed_read_kmer_clouds.append(indexed_cloud)
            indexed_reads_kmer_clouds[r_id] = indexed_read_kmer_clouds
        return indexed_reads_kmer_clouds

    if verbose:
        print("Indexing")
    kmer_index = {kmer: i for i, kmer in enumerate(kmers)}
    indexed_reads_kmer_clouds = index_clouds(reads_kmer_clouds, min_n, max_n)

    if verbose:
        print("Inferring distances")
    dist_cnt = {}

    for dist in range(min_d, max_d + 1):
        dist_cnt[dist] = [defaultdict(int) for i in range(len(kmers))]
        dt = dist_cnt[dist]
        for n, (r_id, kmer_clouds) in \
                enumerate(indexed_reads_kmer_clouds.items()):
            if verbose:
                print(dist,
                      n + min_n,
                      len(kmer_clouds),
                      len(indexed_reads_kmer_clouds.items()) + min_n)
            for i, i_cloud in enumerate(kmer_clouds[:-dist]):
                j_cloud = kmer_clouds[i+dist]
                for i_index in i_cloud:
                    for j_index in j_cloud:
                        # assert i_index != j_index
                        if i_index != j_index:
                            dt[i_index][j_index] += 1
    return dist_cnt, kmer_index


def filter_dist_tuples(dist_cnt, min_coverage, rel_threshold=0.8):
    candidate_edges = {}
    for dist, dt in dist_cnt.items():
        for i_index in range(len(dt)):
            for j_index in dt[i_index]:
                freq = dt[i_index][j_index]
                if freq >= min_coverage:
                    candidate_edges[(i_index, j_index, dist)] = freq

    selected_kmers = []
    selected_edges = []
    for (i_index, j_index, cand_dist), freq in candidate_edges.items():
        all_occ = sum(dist_cnt[dist][i_index][j_index] for dist in dist_cnt)
        if freq / all_occ >= rel_threshold:
            selected_kmers.append(i_index)
            selected_kmers.append(j_index)
            selected_edges.append((cand_dist, i_index, j_index, freq))
    selected_kmers = set(selected_kmers)
    return selected_kmers, selected_edges


def output_results(kmer_index, min_coverage,
                   unique_kmers_ind, dist_edges, outdir):
    kmer_index_reversed = {}
    for kmer, index in kmer_index.items():
        kmer_index_reversed[index] = kmer

    kmers_out_fn = \
        os.path.join(outdir, f'unique_kmers_min_edge_cov_{min_coverage}.txt')
    with open(kmers_out_fn, 'w') as f:
        kmers = [kmer_index_reversed[index] for index in unique_kmers_ind]
        kmers = sorted(list(kmers))
        for kmer in kmers:
            print(kmer, file=f)
    edges_out_fn = \
        os.path.join(outdir, f'unique_edges_min_edge_cov_{min_coverage}.txt')
    with open(edges_out_fn, 'w') as f:
        for t in dist_edges:
            print(t[0], kmer_index_reversed[t[1]],
                  kmer_index_reversed[t[2]], t[3],
                  file=f)


def main():
    params = parse_args()
    smart_makedirs(params.outdir)

    reads_ncrf_report = NCRF_Report(params.ncrf)
    rare_kmers = get_rare_kmers(reads_ncrf_report,
                                k=params.k,
                                bottom=params.bottom,
                                top=params.top,
                                coverage=params.coverage,
                                kmer_survival_rate=params.kmer_survival_rate,
                                max_nonuniq=params.max_nonuniq,
                                verbose=params.verbose)

    reads_kmer_clouds = get_reads_kmer_clouds(reads_ncrf_report,
                                              n=1,
                                              k=params.k,
                                              genomic_kmers=rare_kmers)

    dist_cnt, kmer_index = get_kmer_dist_map(reads_kmer_clouds,
                                             rare_kmers,
                                             min_n=params.min_nreads,
                                             max_n=params.max_nreads,
                                             min_d=params.min_distance,
                                             max_d=params.max_distance,
                                             verbose=params.verbose)

    unique_kmers_ind, dist_edges = \
        filter_dist_tuples(dist_cnt, min_coverage=params.min_coverage)

    output_results(kmer_index=kmer_index,
                   min_coverage=params.min_coverage,
                   unique_kmers_ind=unique_kmers_ind,
                   dist_edges=dist_edges,
                   outdir=params.outdir)


if __name__ == "__main__":
    main()
