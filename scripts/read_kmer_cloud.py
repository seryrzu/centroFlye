#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

import math
from collections import defaultdict, Counter
from ncrf_parser import NCRF_Report
from utils.bio import RC


class ReadKMerCloud:
    def __init__(self, kmers, r_id):
        self.r_id = r_id
        self.kmers = kmers
        self.all_kmers = set()
        for v in self.kmers:
            self.all_kmers |= v

    @classmethod
    def fromNCRF_record(cls, ncrf_record, n, k, genomic_kmers=None):
        r_id = ncrf_record.r_id
        mas = ncrf_record.get_motif_alignments(n=n)

        all_kmers = Counter()
        for ma in mas:
            r_al = ma.r_al.upper().replace('-', '')
            for i in range(len(r_al) - k + 1):
                r_kmer = r_al[i:i+k]
                if genomic_kmers is None or (genomic_kmers is not None and r_kmer in genomic_kmers):
                    all_kmers[r_kmer] += 1
                # RC_r_kmer = RC(r_kmer)
                # if genomic_kmers is None or (genomic_kmers is not None and RC_r_kmer in genomic_kmers):
                #     all_kmers[RC(r_kmer)] += 1

        kmer_clouds = []
        for ma in mas:
            kmers, rc_kmers = set(), set()
            r_al = ma.r_al.upper().replace('-', '')
            for i in range(len(r_al) - k + 1):
                r_kmer = r_al[i:i+k]
                # RC_r_kmer = RC(r_kmer)
                if all_kmers[r_kmer] == 1: #and ncrf_record.strand == '+':
                    kmers.add(r_kmer)
                # if all_kmers[RC_r_kmer] == 1: #and ncrf_record.strand == '-':
                #     kmers.add(RC_r_kmer)
            kmer_clouds.append(kmers)

        # if ncrf_record.strand == '-':
        #     kmers = kmers[::-1]

        return cls(kmers=kmer_clouds, r_id=r_id)

def get_reads_kmer_clouds(ncrf_report, n, k, genomic_kmers=None):
    kmer_clouds = {}
    for r_id, record in ncrf_report.records.items():
        kmer_clouds[r_id] = ReadKMerCloud.fromNCRF_record(record, n=n, k=k, genomic_kmers=genomic_kmers)
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
                set(kmer for kmer in kmers if max_mult >= all_kmers[kmer] >= min_mult)
    return kmer_clouds


def get_all_kmers(kmer_clouds):
    all_kmers = []
    for r_id, kmer_cloud in kmer_clouds.items():
        all_kmers +=  kmer_clouds.all_kmers
    all_kmers.sort()
    assert len(list(set(all_kmers))) == len(all_kmers)
    return all_kmers
