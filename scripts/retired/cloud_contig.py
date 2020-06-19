# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import defaultdict, Counter


class CloudContig:
    def __init__(self, min_cloud_kmer_freq):
        self.max_pos = 0
        self.min_cloud_kmer_freq = max(1, min_cloud_kmer_freq)

        self.clouds = defaultdict(Counter)
        self.freq_clouds = defaultdict(set)
        self.freq_kmers = set()
        self.kmer_positions = defaultdict(set)
        self.read_positions = {}
        self.coverage = defaultdict(int)

    def update_max_pos(self):
        if len(self.clouds):
            self.max_pos = max(self.clouds.keys())
        else:
            self.max_pos = 0

    def add_read(self, read_kmer_clouds, position):
        self.read_positions[read_kmer_clouds.r_id] = position
        new_freq_kmers = []
        for i, cloud in enumerate(read_kmer_clouds.kmers):
            self.coverage[i + position] += 1
            self.clouds[i + position]
            for kmer in cloud:
                self.kmer_positions[kmer].add(i+position)
                self.clouds[i+position][kmer] += 1
                if self.clouds[i+position][kmer] == self.min_cloud_kmer_freq:
                    self.freq_clouds[i+position].add(kmer)
                    self.freq_kmers.add(kmer)
                    new_freq_kmers.append((kmer, i+position))
        self.update_max_pos()
        assert len(set(new_freq_kmers)) == len(new_freq_kmers)
        return new_freq_kmers

    def calc_rough_inters_score(self, read_kmer_cloud):
        return len(read_kmer_cloud.all_kmers & self.freq_kmers)

    def calc_inters_score(self, read_kmer_cloud,
                          min_position=0, max_position=None,
                          min_unit=2, min_inters=10,
                          verbose=False):
        if max_position is None:
            max_position = self.max_pos
        best_score, best_pos = (0, 0), None
        kmers = read_kmer_cloud.kmers
        positions = [pos for pos in range(min_position, max_position + 1)]
        for pos in positions:
            score = [0, 0]
            max_i = min(self.max_pos-pos+1, len(kmers))
            for i in range(max_i):
                assert pos + i <= self.max_pos
                freq_cloud = self.freq_clouds[pos + i]
                inters = freq_cloud & set(kmers[i])
                if verbose:
                    print(pos, i, inters, len(inters),
                          len(freq_cloud), len(kmers[i]))
                score[0] += len(inters) >= 1
                score[1] += len(inters)
            if verbose:
                print(f'pos: {pos}, i: {i}, score: {score}')
            score = tuple(score)
            # we want to take the rightmost best, so >= instead of >
            if score[0] >= min_unit and \
                    score[1] >= min_inters and \
                    score >= best_score:
                best_score = score
                best_pos = pos
        return best_score, best_pos

    def get_spread_kmers(self, max_npos=5):
        kmers = []
        for kmer in self.freq_kmers:
            positions = self.kmer_positions[kmer]
            if len(positions) > max_npos:
                kmers.append(kmer)
        return set(kmers)


def update_mapping_scores(cloud_contig, kmers2pos, freq_kmers, scores=None):
    if scores is None:
        scores = defaultdict(lambda: defaultdict(Counter))
    for kmer, cc_kmer_pos in freq_kmers:
        if kmer in kmers2pos:
            for r_id, pos in kmers2pos[kmer]:
                if cc_kmer_pos >= pos:
                    scores[r_id][cc_kmer_pos-pos][pos] += 1
    return scores


def map_reads(cloud_contig, reads_kmer_clouds,
              threshold=(5, 10), verbose=False):
    scores = {}
    pos = {}
    for i, (r_id, kmer_clouds) in enumerate(reads_kmer_clouds.items()):
        max_pos = cloud_contig.max_pos - len(kmer_clouds.kmers) + 1
        best_score, best_pos = \
            cloud_contig.calc_inters_score(kmer_clouds,
                                           max_position=max_pos)
        if (best_pos == 0) or best_score > threshold:
            scores[r_id] = best_score
            pos[r_id] = best_pos
            if verbose:
                print(f"{i+1} / {len(reads_kmer_clouds)},"
                      f"pos = {pos[r_id]}, score = {scores[r_id]},"
                      f"id = {r_id}")
    return pos, scores


def map_reads_fast(cloud_contig, reads_kmer_clouds,
                   threshold=(5, 10), verbose=False, debug=False):
    kmers2pos = defaultdict(list)
    for r_id, kmer_clouds in reads_kmer_clouds.items():
        for i, cloud in enumerate(kmer_clouds.kmers):
            for kmer in cloud:
                kmers2pos[kmer].append((r_id, i))

    freq_kmers = []
    for kmer in cloud_contig.freq_kmers:
        for pos in cloud_contig.kmer_positions[kmer]:
            freq_kmers.append((kmer, pos))

    scores = update_mapping_scores(cloud_contig, kmers2pos, freq_kmers)
    positions = {}
    for i, (r_id, kmer_clouds) in enumerate(reads_kmer_clouds.items()):
        best_score, best_pos = (0, 0), None
        for pos, score in scores[r_id].items():
            if pos + len(kmer_clouds.kmers) > len(cloud_contig.clouds):
                continue
            score = (len(score), sum(score.values()))
            if score[0] < threshold[0] or score[1] < threshold[1]:
                continue
            if (score > best_score) or \
                    (score == best_score and pos > best_pos):
                best_pos = pos
                best_score = score
        if best_pos is not None:
            positions[r_id] = best_pos
            if debug:
                max_pos = cloud_contig.max_pos - len(kmer_clouds.kmers) + 1
                best_score_slow, best_pos_slow = \
                    cloud_contig.calc_inters_score(kmer_clouds,
                                                   max_position=max_pos,
                                                   min_unit=threshold[0],
                                                   min_inters=threshold[1])
                if best_score_slow != best_score or best_pos_slow != best_pos:
                    print(r_id[:8], best_score,
                          best_score_slow, best_pos, best_pos_slow)
    return positions, scores
