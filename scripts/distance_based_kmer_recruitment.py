#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

import argparse
import os
import math
from collections import defaultdict, Counter
from utils.os_utils import smart_makedirs
from ncrf_parser import NCRF_Report
from read_kmer_cloud import get_reads_kmer_clouds, filter_reads_kmer_clouds


def get_kmer_freqs_from_ncrf_reports(reads_ncrf_report, k=19, verbose=False):
    all_kmers = defaultdict(int)
    for i, (r_id, record) in enumerate(reads_ncrf_report.records.items()):
        if i % 100 == 0 and verbose:
            print(i + 1, len(reads_ncrf_report.records))
        r_al = record.r_al
        r_al = r_al.replace('-', '')

        read_kmers = defaultdict(int)
        for i in range(len(r_al)-k+1):
            kmer = r_al[i:i+k]
            read_kmers[kmer] += 1
        for kmer, freq in read_kmers.items():
            if freq == 1:
                all_kmers[kmer] += 1
    return all_kmers


def get_kmer_dist_map(reads_kmer_clouds, kmers, min_n, max_n, verbose=False):
    dist_tuples = Counter()
    kmer_index = {}
    index = 0
    for kmer in kmers:
        kmer_index[kmer] = index
        index += 1

    if verbose:
        print("Indexing")
    indexed_reads_kmer_clouds = {}
    for n, (r_id, kmer_clouds) in enumerate(reads_kmer_clouds.items()):
        if n <= min_n:
            continue
        if n >= max_n:
            break
        indexed_read_kmer_clouds = []
        kmers = kmer_clouds.kmers
        if verbose:
            print(n, len(kmers), len(reads_kmer_clouds.items()))
        for cloud in kmers:
            indexed_cloud = []
            for kmer in cloud:
                indexed_cloud.append(kmer_index[kmer])
            indexed_read_kmer_clouds.append(indexed_cloud)
        indexed_reads_kmer_clouds[r_id] = indexed_read_kmer_clouds

    if verbose:
        print("Inferring distances")
    for n, (r_id, kmer_clouds) in enumerate(indexed_reads_kmer_clouds.items()):
        if verbose:
            print(n, len(kmer_clouds), len(indexed_reads_kmer_clouds.items()))
        for i in range(len(kmer_clouds)):
            i_indexed_cloud = kmer_clouds[i]
            for j in range(i + 1, len(kmer_clouds)):
                j_indexed_cloud = kmer_clouds[j]
                for ikmer_index in i_indexed_cloud:
                    for jkmer_index in j_indexed_cloud:
                        #kmer_pair = min(ikmer, jkmer), max(ikmer, jkmer)
                        if ikmer_index != jkmer_index:
                            dist_tuples[(j-i, ikmer_index, jkmer_index)] += 1
    return dist_tuples, kmer_index


def filter_dist_tuples(dist_tuples, min_coverage, max_dist=300, rel_threshold=0.8):
    filtered_dist_tuples = {k: v for (k, v) in dist_tuples.items() if v >= min_coverage}
    selected_kmers = []
    selected_tuples = []
    for i, ((dist, kmer1, kmer2), val) in enumerate(filtered_dist_tuples.items()):
        all_occ = sum(dist_tuples[(sec_dist, kmer1, kmer2)] for sec_dist in range(max_dist + 1))
        if val / all_occ >= rel_threshold:
            #print(dist, kmer1, kmer2, val)
            selected_kmers.append(kmer1)
            selected_kmers.append(kmer2)
            selected_tuples.append((dist, kmer1, kmer2, val))
    selected_kmers = set(selected_kmers)
    return selected_kmers, selected_tuples


def from_filtered_kmer_dist_map_to_kmers(filtered_dist_tuples):
    kmer_pair2dist = {}
    for i, ((dist, kmer1, kmer2), val) in enumerate(filtered_dist_tuples.items()):
        all_occ = sum(dist_tuples[(sec_dist, kmer1, kmer2)] for sec_dist in range(1, 100))
        if val / all_occ > 0.95:
            #print(dist, kmer1, kmer2, val)
            kmer_pair2dist[(kmer1, kmer2)] = dist

    diff_distance = {}
    for (kmer1, kmer2), read_dist in kmer_pair2dist.items():
        if kmer1 not in genomic_kmer_positions or \
            kmer2 not in genomic_kmer_positions:
            continue
        genome_dist = abs(genomic_kmer_positions[kmer1] - genomic_kmer_positions[kmer2])
        diff_distance[(kmer1, kmer2)] = read_dist - genome_dist


    dist_kmer_score = defaultdict(lambda: [0, 0])

    for (kmer1, kmer2), dist in diff_distance.items():
        if dist == 0:
            dist_kmer_score[kmer1][0] += 1
            dist_kmer_score[kmer2][0] += 1
        else:
            dist_kmer_score[kmer1][1] += 1
            dist_kmer_score[kmer2][1] += 1

    kmers_good_dist_score = set([kmer for kmer, sc in dist_kmer_score.items() \
                     if (sc[0] >= 5 and sc[1] / (sc[0] + sc[1]) < 0.1)])

    return kmers_good_dist_score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncrf', help='NCRF report on reads', required=True)
    parser.add_argument('--coverage', help='Average coverage of the dataset', type=int, required=True)
    parser.add_argument('--min-coverages', help='Comma separated list of minCov thresholds', required=True)
    parser.add_argument('--outdir', help='Output directory', required=True)
    parser.add_argument('-k', type=int, default=19)
    parser.add_argument('--min-nreads', type=int, default=0)
    parser.add_argument('--max-nreads', type=int, default=math.inf)
    params = parser.parse_args()

    bottom, top = 0.9, 3
    mean_survival_rate = 0.34
    min_coverages = [int(x) for x in params.min_coverages.split(',')]
    smart_makedirs(params.outdir)

    reads_ncrf_report = NCRF_Report(params.ncrf)
    all_kmers = get_kmer_freqs_from_ncrf_reports(reads_ncrf_report)

    filtered_kmers = {k: v for k, v in all_kmers.items() \
                    if bottom*params.coverage*mean_survival_rate <= v <= top*params.coverage*mean_survival_rate}
    filtered_kmers_set = set(filtered_kmers.keys())
    print(f'# rare kmers: {len(filtered_kmers_set)}')

    reads_kmer_clouds = get_reads_kmer_clouds(reads_ncrf_report, n=1, k=params.k,
                                              genomic_kmers=filtered_kmers_set)
    dist_tuples, kmer_index = get_kmer_dist_map(reads_kmer_clouds,
                                                filtered_kmers_set,
                                                min_n=params.min_nreads,
                                                max_n=params.max_nreads,
                                                verbose=True)

    unique_kmers = {}
    selected_unique_tuples = {}
    for min_coverage in min_coverages:
        uk, sut = filter_dist_tuples(dist_tuples, min_coverage=min_coverage)
        unique_kmers[min_coverage] = uk
        selected_unique_tuples[min_coverage] = sut

    kmer_index_reversed = {}
    for kmer, index in kmer_index.items():
        kmer_index_reversed[index] = kmer

    for min_coverage, kmer_indexes in unique_kmers.items():
        with open(os.path.join(params.outdir, f'unique_kmers_min_edge_cov_{min_coverage}_RC.txt'), 'w') as f:
            kmers = [kmer_index_reversed[index] for index in kmer_indexes]
            kmers = sorted(list(kmers))
            for kmer in kmers:
                print(kmer, file=f)
        with open(os.path.join(params.outdir, f'unique_tuples_min_edge_cov_{min_coverage}_RC.txt'), 'w') as f:
            sut = selected_unique_tuples[min_coverage]
            for t in sut:
                print(t[0], kmer_index_reversed[t[1]], kmer_index_reversed[t[2]], t[3], file=f)



if __name__ == "__main__":
    main()
