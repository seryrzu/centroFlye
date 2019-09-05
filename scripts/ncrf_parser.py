# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import regex as re
from utils.bio import RC
from collections import defaultdict, namedtuple
import numpy as np


class NCRF_Report:
    class NCRF_Record:
        def __init__(self, r_id, r_len, r_al_len, r_st, r_en, r_al,
                     motif, strand, m_al_len, al_score, m_al):
            self.r_id = r_id
            self.r_len = int(r_len)
            self.r_al_len = int(r_al_len)
            self.r_st = int(r_st)
            self.r_en = int(r_en)
            self.r_al = r_al
            self.motif = motif
            self.strand = strand
            self.m_al_len = int(m_al_len)
            self.al_score = int(al_score)
            self.m_al = m_al

        def get_motif_alignments(self, n=1, overlapped=False):
            motif_alignments = []
            motif = self.motif
            # we already RC the alignment on reading step
            # if self.strand == '-':
            #     motif = RC(motif)
            motif = ''.join(f'{base}([-]*)' for base in motif)
            motif *= n
            m_al = self.m_al.upper()
            all_matches = re.finditer(motif, m_al, overlapped=overlapped)
            all_matches = list(all_matches)
            if len(all_matches) == 0:
                return []

            # [start, end)
            MotifAlignment = \
                namedtuple('MotifAlignment',
                           ['r_id', 'start', 'end', 'r_al', 'm_al'])

            coords = [match.start() for match in all_matches]
            coords.append(all_matches[-1].end())
            if coords[0] > len(self.motif) * 0.2:
                coords.insert(0, 0)
            if coords[-1] < len(self.r_al) - len(self.motif) * 0.2:
                coords.append(len(self.r_al))
            for st, en in zip(coords[:-1], coords[1:]):
                ma = MotifAlignment(r_id=self.r_id,
                                    start=st, end=en,
                                    r_al=self.r_al[st:en],
                                    m_al=self.m_al[st:en])
                motif_alignments.append(ma)
            return motif_alignments

    def __init__(self, report_fn, min_record_len=5000):
        self.records = {}
        self.positions_all_alignments = defaultdict(list)
        # self.longest_alignment_index = {}
        self.read_lens = {}
        with open(report_fn, 'r') as f:
            lines = [x.strip() for x in f.readlines()]
            lines = list(filter(lambda x: len(x) and x[0] != '#', lines))
            read_records = [lines[i:i+2] for i in range(0, len(lines), 2)]
            seen_r_ids = []
            for read_record in read_records:
                fst, snd = read_record

                fst_search = re.search('^([^ ]+)\s+(\d+)\s+(\d+)bp\s+(\d+)-(\d+)\s+(.+)$', fst)
                snd_search = re.search('^([^+-]+)([+-])\s+(\d+)bp\s+score=(\d+)\s+(.+)$', snd)
                r_id, r_len, r_al_len, r_st, r_en, r_al = fst_search.groups()
                motif, strand, m_al_len, al_score, m_al = snd_search.groups()
                r_al_len = int(r_al_len)
                r_st, r_en = int(r_st), int(r_en)
                r_len = int(r_len)
                m_al_len = int(m_al_len)
                al_score = int(al_score)
                seen_r_ids.append(r_id)

                # if r_al_len < min_alignment_len:
                #     continue

                self.positions_all_alignments[r_id].append((r_st, r_en, strand))
                self.read_lens[r_id] = r_len

                if r_id not in self.records or self.records[r_id].r_al_len < r_al_len:
                    if r_al_len < min_record_len:
                        continue
                    # if strand is negative, we reverse the alignment but maintain information
                    # that the stand of alignment is originally negative
                    if strand == '-':
                        # strand = '+'
                        r_st, r_en = r_len - r_en, r_len - r_st
                        r_al = RC(r_al)
                        m_al = RC(m_al)
                    self.records[r_id] = self.NCRF_Record(r_id=r_id,
                                                          r_len=r_len,
                                                          r_al_len=r_al_len,
                                                          r_st=r_st,
                                                          r_en=r_en,
                                                          r_al=r_al,
                                                          motif=motif,
                                                          strand=strand,
                                                          m_al_len=m_al_len,
                                                          al_score=al_score,
                                                          m_al=m_al)
        for r_id in self.positions_all_alignments:
            self.positions_all_alignments[r_id].sort()
        seen_r_ids = list(set(seen_r_ids))
        self.discarded_reads = []
        for r_id in seen_r_ids:
            if r_id not in self.records:
                self.discarded_reads.append(r_id)

    def classify(self, large_threshold, small_threshold=1000):
        prefix_reads, suffix_reads, internal_reads = [], [], []
        for r_id, record in self.records.items():
            r_len = self.read_lens[r_id]
            if record.strand == '+':
                first_alignment = self.positions_all_alignments[r_id][0]
                last_alignment = self.positions_all_alignments[r_id][-1]
                left_pos = first_alignment[0]
                right_pos = last_alignment[1]
            else:
                first_alignment = self.positions_all_alignments[r_id][-1]
                last_alignment = self.positions_all_alignments[r_id][0]
                left_pos = r_len - first_alignment[1]
                right_pos = r_len - last_alignment[0]
            if left_pos > large_threshold \
                    and right_pos > r_len - small_threshold \
                    and right_pos == self.records[r_id].r_en:

                prefix_reads.append(r_id)
            elif right_pos < r_len - large_threshold \
                    and left_pos < small_threshold \
                    and left_pos == self.records[r_id].r_st:
                suffix_reads.append(r_id)
            else:
                internal_reads.append(r_id)
        return prefix_reads, internal_reads, suffix_reads

    def get_efficiency(self):
        efficiency = {}
        total_length = 0
        total_used_length = 0
        for r_id, alignments in self.positions_all_alignments.items():
            all_alignments_len = \
                sum(alignment[1] - alignment[0] + 1 \
                    for alignment in alignments)
            total_length += all_alignments_len
            if r_id not in self.records:
                efficiency[r_id] = 0
            else:
                record = self.records[r_id]
                record_len = record.r_en - record.r_st + 1
                total_used_length += record_len
                efficiency[r_id] = record_len / all_alignments_len
        global_efficiency = total_used_length / total_length
        return efficiency, global_efficiency


    '''
    :param n = number of motifs stuck together (seek for motif*n)
    '''
    def get_motif_alignments(self, n=1):
        motifs_alignments = {}
        for read_id, record in self.records.items():
            motifs_alignments[read_id] = record.get_motif_alignments(n=n)
        return motifs_alignments


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncrf", help="NCRF report")
    parser.add_argument('--prefix-threshold',
                        help='Min pre/suffix length for read classification',
                        default=50000,
                        type=int)
    params = parser.parse_args()

    ncrf_report = NCRF_Report(params.ncrf)
    print(f"# records: {len(ncrf_report.records)}")
    prefix_reads, internal_reads, suffix_reads = \
        ncrf_report.classify(large_threshold=params.prefix_threshold)
    print(f"# prefix reads {len(prefix_reads)}")
    for read in prefix_reads:
        print(read[:8])
    print(f"# internal reads {len(internal_reads)}")
    print(f"# suffix reads {len(suffix_reads)}")
    for read in suffix_reads:
        print(read[:8])

    efficiency, global_efficiency = ncrf_report.get_efficiency()
    print(f'Mean efficiency = {np.mean(list(efficiency.values()))}')
    print(f'Median efficiency = {np.median(list(efficiency.values()))}')
    print(f'Global efficiency = {global_efficiency}')


if __name__ == "__main__":
    main()
