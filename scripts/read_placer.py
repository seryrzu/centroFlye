#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

import argparse
import os
import copy
import sys
import heapq
from collections import Counter, defaultdict
from ncrf_parser import NCRF_Report
from utils.bio import RC
from utils.os_utils import smart_makedirs

from cloud_contig import CloudContig
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
        self.position_outfile = os.path.join(self.params.outdir, 'read_positions.csv')

    def add_left_PT_reads(self, PT_reads, reads_kmer_clouds):
        with open(self.position_outfile, 'w') as f:
            for r_id in PT_reads:
                read_kmer_clouds = reads_kmer_clouds[r_id]
                self.cloud_contig.add_read(read_kmer_clouds, position=0)
                print(r_id, 0, file=f)

    def add_FT_reads(self, FT_reads, reads_kmer_clouds):
        #unused_FT_reads = set(r_id for r_id in FT_reads if len(reads_kmer_clouds[r_id].kmers) >= 10)
        unused_FT_reads = set(FT_reads)
        n_reads = len(unused_FT_reads)
        with open(self.position_outfile, 'a') as f:
            while len(unused_FT_reads):
                rough_scores = {}
                for FT_read in unused_FT_reads:
                    kmer_cloud = reads_kmer_clouds[FT_read]
                    rough_score = self.cloud_contig.calc_rough_inters_score(kmer_cloud)
                    rough_scores[FT_read] = rough_score
                # selected_FT_reads = heapq.nlargest(50, rough_scores,
                #                                    key=rough_scores.get)
                # selected_FT_reads = [k for k in selected_FT_reads if rough_scores[k] >= 30]
                # selected_FT_reads = [k for k in unused_FT_reads if rough_scores[k] >= 30]


                best_score, best_position, best_read = (-1, -1), None, None
                # for FT_read in selected_FT_reads:
                for FT_read in unused_FT_reads:
                    kmer_clouds = reads_kmer_clouds[FT_read]
                    score, position = self.cloud_contig.calc_inters_score(read_kmer_cloud=kmer_clouds)
                    if score == (0, 0):
                        continue
                    # print(score, position, FT_read)
                    if score > best_score or (score == best_score and position > best_position):
                        best_score = score
                        best_position = position
                        best_read = FT_read
                if best_read is None:
                    print(f"Unused reads {len(unused_FT_reads)}, {n_reads}, {len(unused_FT_reads) / n_reads}")
                    for read in unused_FT_reads:
                        print(read, None, file=f)
                    return
                print(best_score, best_position, best_read)
                print("")
                print(best_read, best_position, best_score[0], best_score[1], file=f)
                best_read_cloud = reads_kmer_clouds[best_read]
                self.cloud_contig.add_read(best_read_cloud, position=best_position)
                unused_FT_reads.remove(best_read)

    def add_right_PT_reads(self, PT_reads, reads_kmer_clouds):
        reads_kmer_clouds = copy.deepcopy(reads_kmer_clouds)
        # reverse kmers clouds to add to new cloud_contig
        for r_id in PT_reads:
            read_kmer_cloud = reads_kmer_clouds[r_id]
            kmers = read_kmer_cloud.kmers
            reads_kmer_clouds[r_id].kmers = kmers[::-1]
        # create new cloud contig
        PT_cloud_contig = CloudContig(self.params.min_cloud_kmer_freq)
        for r_id in PT_reads:
            read_kmer_clouds = reads_kmer_clouds[r_id]
            PT_cloud_contig.add_read(read_kmer_clouds, position=0)
        # Get kmer clouds from new cloud contig and reverse them back
        PT_clouds = []
        for i in range(PT_cloud_contig.max_pos):
            PT_clouds.append(PT_cloud_contig.clouds[i])
        PT_clouds = ReadKMerCloud(kmers=PT_clouds[::-1], r_id='cloud_contig')
        # Find the best position for new kmer clouds
        min_position = self.cloud_contig.max_pos - PT_cloud_contig.max_pos
        best_score, best_pos = self.cloud_contig.calc_inters_score(PT_clouds, min_position)
        print(best_score, best_pos)

        # Add reads with corresponding positions
        with open(self.position_outfile, 'a') as f:
            for r_id in PT_reads:
                read_kmer_clouds = reads_kmer_clouds[r_id]
                position = best_pos + PT_cloud_contig.max_pos - len(read_kmer_clouds.kmers)
                self.cloud_contig.add_read(read_kmer_clouds, position=position)
                print(r_id, position)
                print(r_id, position, file=f)


    def run(self):
        print("Reading kmer clouds from reads")
        reads_kmer_clouds = get_reads_kmer_clouds(self.ncrf_report,
                                                  n=self.params.n_motif,
                                                  k=self.params.k_cloud,
                                                  genomic_kmers=self.genomic_kmers)
        print("Filtering kmer clouds from reads")
        reads_kmer_clouds = filter_reads_kmer_clouds(reads_kmer_clouds,
                                                     min_mult=self.params.min_kmer_mult)

        left_PT_reads, FT_reads, right_PT_reads = self.ncrf_report.classify()

        print(f'Left: {len(left_PT_reads)}')
        print(f'FT: {len(FT_reads)}')
        print(f'Right: {len(right_PT_reads)}')

        print("Adding left reads")
        self.add_left_PT_reads(left_PT_reads, reads_kmer_clouds)
        print(self.cloud_contig.max_pos)
        for pos, cloud in self.cloud_contig.clouds.items():
            freq = {k: v for k, v in cloud.items() if v >= self.params.min_cloud_kmer_freq}
            print(pos, len(freq), sum(freq.values()))
        print("")

        self.add_FT_reads(FT_reads, reads_kmer_clouds)
        print(self.cloud_contig.max_pos)
        for pos, cloud in self.cloud_contig.clouds.items():
            freq = {k: v for k, v in cloud.items() if v >= self.params.min_cloud_kmer_freq}
            print(pos, len(freq), sum(freq.values()))
        print("")

        print("\nNow right PT")
        # self.add_right_PT_reads(right_PT_reads, reads_kmer_clouds)
        self.add_FT_reads(right_PT_reads, reads_kmer_clouds)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncrf', help='NCRF report on reads', required=True)
    parser.add_argument('--genomic-kmers', help='Unique genomic kmers if known', required=True)
    parser.add_argument('--n-motif', help='Number of motifs stuck together', default=1, type=int)
    parser.add_argument('--k-cloud', help='Size of k-mer for k-mer cloud', default=19, type=int)
    parser.add_argument('--min-cloud-kmer-freq', help='Minimal frequency of a kmer in the cloud', default=2, type=int)
    parser.add_argument('--min-kmer-mult', help='Minimal frequency of a kmer in input', default=2, type=int)
    parser.add_argument('--outdir', help='Output directory', required=True)
    params = parser.parse_args()

    ReadPlacer(params).run()


if __name__ == "__main__":
    main()
