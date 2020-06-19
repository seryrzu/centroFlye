# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
from collections import defaultdict

from ncrf_parser import NCRF_Report
from utils.os_utils import smart_makedirs
from cloud_contig import CloudContig, update_mapping_scores
from read_kmer_cloud import get_reads_kmer_clouds, filter_reads_kmer_clouds


class ReadPlacer:
    def __init__(self, params):
        self.params = params
        self.ncrf_report = NCRF_Report(params.ncrf)
        self.cloud_contig = CloudContig(params.min_cloud_kmer_freq)
        if params.genomic_kmers is not None:
            kmers = []
            with open(params.genomic_kmers) as f:
                for line in f:
                    kmers.append(line.strip())
            self.genomic_kmers = set(kmers)
        else:
            self.genomic_kmers = None
        smart_makedirs(params.outdir)
        self.position_outfile = \
            os.path.join(self.params.outdir, 'read_positions.csv')

    def reset_cloud_contig(self):
        self.cloud_contig = CloudContig(self.params.min_cloud_kmer_freq)

    def add_prefix_reads(self, prefix_reads, reads_kmer_clouds):
        with open(self.position_outfile, 'w') as f:
            for r_id in prefix_reads:
                read_kmer_clouds = reads_kmer_clouds[r_id]
                self.cloud_contig.add_read(read_kmer_clouds, position=0)
                print(r_id, 0, file=f)

    def add_reads(self, reads, reads_kmer_clouds,
                  min_unit, min_inters, min_prop=3):
        kmers2pos = defaultdict(list)
        for r_id in reads:
            kmer_clouds = reads_kmer_clouds[r_id]
            for i, cloud in enumerate(kmer_clouds.kmers):
                for kmer in cloud:
                    kmers2pos[kmer].append((r_id, i))

        unused_reads = set(reads)
        n_reads = len(unused_reads)
        scores = None
        freq_kmers = []
        for kmer in self.cloud_contig.freq_kmers:
            for pos in self.cloud_contig.kmer_positions[kmer]:
                freq_kmers.append((kmer, pos))
        with open(self.position_outfile, 'a') as f:
            while len(unused_reads):
                scores = update_mapping_scores(self.cloud_contig, kmers2pos,
                                               freq_kmers=freq_kmers,
                                               scores=scores)
                best_score, best_position, best_read = (-1, -1), None, None
                for r_id in unused_reads:
                    for pos in scores[r_id]:
                        score = scores[r_id][pos]
                        score = (len(score), sum(score.values()))
                        if (score > best_score and
                                score[0] >= min_unit and
                                score[0] * min_prop <= score[1] and
                                score[1] >= min_inters) or \
                            (score == best_score and pos > best_position) or \
                                (score == best_score and
                                 pos == best_position and
                                 r_id < best_read):
                            best_score = score
                            best_position = pos
                            best_read = r_id
                if best_read is None:
                    print(f"Unused reads {len(unused_reads)}, {n_reads}, "
                          f"{len(unused_reads) / n_reads}")
                    for read in unused_reads:
                        print(read, None, file=f)
                    return
                print(best_score, best_position, best_read)
                print("")
                print(best_read, best_position,
                      best_score[0], best_score[1],
                      file=f)
                best_read_cloud = reads_kmer_clouds[best_read]

                freq_kmers = self.cloud_contig.add_read(best_read_cloud,
                                                        position=best_position)
                unused_reads.remove(best_read)

    def run(self):
        left_PT_reads, FT_reads, right_PT_reads = \
            self.ncrf_report.classify(
                large_threshold=self.params.prefix_threshold)

        print(f'Left: {len(left_PT_reads)}')
        print(f'FT: {len(FT_reads)}')
        print(f'Right: {len(right_PT_reads)}')

        print("Reading kmer clouds from reads")
        reads_kmer_clouds = \
            get_reads_kmer_clouds(self.ncrf_report,
                                  n=self.params.n_motif,
                                  k=self.params.k_cloud,
                                  genomic_kmers=self.genomic_kmers)
        print("Filtering kmer clouds from reads")
        reads_kmer_clouds = \
            filter_reads_kmer_clouds(reads_kmer_clouds,
                                     min_mult=self.params.min_kmer_mult)
        print("Adding prefix reads")
        self.add_prefix_reads(left_PT_reads, reads_kmer_clouds)
        print(self.cloud_contig.max_pos)

        print("Adding inner reads")
        self.add_reads(FT_reads, reads_kmer_clouds,
                       min_unit=self.params.min_unit,
                       min_inters=self.params.min_inters)
        print(self.cloud_contig.max_pos)

        print("\nNow adding suffix reads")
        self.add_reads(right_PT_reads, reads_kmer_clouds,
                       min_unit=self.params.min_unit,
                       min_inters=self.params.min_inters)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncrf',
                        help='NCRF report on reads',
                        required=True)
    parser.add_argument('--genomic-kmers',
                        help='Unique genomic kmers if known',
                        required=True)
    parser.add_argument('--n-motif',
                        help='Number of motifs stuck together',
                        default=1,
                        type=int)
    parser.add_argument('--k-cloud',
                        help='Size of k-mer for k-mer cloud',
                        default=19,
                        type=int)
    parser.add_argument('--min-cloud-kmer-freq',
                        help='Minimal frequency of a kmer in the cloud',
                        default=2,
                        type=int)
    parser.add_argument('--min-kmer-mult',
                        help='Minimal frequency of a kmer in input',
                        default=2,
                        type=int)
    parser.add_argument('--min-unit',
                        help='Score[0]',
                        default=2,
                        type=int)
    parser.add_argument('--min-inters',
                        help='Score[1]',
                        default=10,
                        type=int)
    parser.add_argument('--prefix-threshold',
                        help='Min pre/suffix length for read classification',
                        default=50000,
                        type=int)
    parser.add_argument('--outdir',
                        help='Output directory',
                        required=True)
    params = parser.parse_args()

    ReadPlacer(params).run()


if __name__ == "__main__":
    main()
