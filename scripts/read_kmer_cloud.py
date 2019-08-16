# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import math
from collections import Counter


class ReadKMerCloud:
    def __init__(self, kmers, r_id):
        self.r_id = r_id
        self.kmers = kmers
        self.all_kmers = []
        for v in self.kmers:
            self.all_kmers += v

    @classmethod
    def fromNCRF_record(cls, ncrf_record, n, k, genomic_kmers):
        r_id = ncrf_record.r_id
        mas = ncrf_record.get_motif_alignments(n=n)

        kmer_clouds = []
        for ma in mas:
            kmers = set()
            r_al = ma.r_al.upper().replace('-', '')
            for i in range(len(r_al) - k + 1):
                r_kmer = r_al[i:i+k]
                if r_kmer in genomic_kmers:
                    kmers.add(r_kmer)
            kmer_clouds.append(kmers)
        return cls(kmers=kmer_clouds, r_id=r_id)


def get_reads_kmer_clouds(ncrf_report, n, k, genomic_kmers=None):
    kmer_clouds = {}
    for r_id, record in ncrf_report.records.items():
        kmer_clouds[r_id] = \
            ReadKMerCloud.fromNCRF_record(record, n=n, k=k,
                                          genomic_kmers=genomic_kmers)
    return kmer_clouds


def filter_reads_kmer_clouds(kmer_clouds, min_mult=2, max_mult=math.inf):
    all_kmers = Counter()
    for r_id, kmer_cloud in kmer_clouds.items():
        for pos, kmers in enumerate(kmer_cloud.kmers):
            for kmer in kmers:
                all_kmers[kmer] += 1
    for r_id, kmer_cloud in kmer_clouds.items():
        for pos, kmers in enumerate(kmer_cloud.kmers):
            kmer_clouds[r_id].kmers[pos] = \
                set(kmer for kmer in kmers
                    if max_mult >= all_kmers[kmer] >= min_mult)
    return kmer_clouds


def get_all_kmers(kmer_clouds):
    all_kmers = []
    for r_id, kmer_cloud in kmer_clouds.items():
        all_kmers += kmer_clouds.all_kmers
    all_kmers.sort()
    assert len(list(set(all_kmers))) == len(all_kmers)
    return all_kmers
