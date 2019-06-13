#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

from collections import defaultdict, Counter


class CloudContig:
    def __init__(self, min_cloud_kmer_freq):
        self.clouds = defaultdict(lambda: Counter()) # map from position to dictionary of kmer counts
        self.freq_clouds = defaultdict(lambda: set())
        self.freq_kmers = set()
        self.prohibited_kmers = set()
        self.kmer_positions = {}
        self.kmer2reads = defaultdict(list)
        self.read_positions = {}
        self.max_pos = 0
        self.min_cloud_kmer_freq = min_cloud_kmer_freq
        self.coverage = defaultdict(int)
        self.pos2kmers = defaultdict(list)

    def update_max_pos(self):
        if len(self.clouds):
            self.max_pos = max(self.clouds.keys())
        else:
            self.max_pos = 0

    def update_freq_clouds(self):
        self.freq_kmers = set()
        for pos, cloud in self.clouds.items():
            self.freq_clouds[pos] = \
                set(k for k, v in cloud.items() if v >= self.min_cloud_kmer_freq)
            self.freq_kmers |= self.freq_clouds[pos]


    def add_read(self, read_kmer_clouds, position):
        self.read_positions[read_kmer_clouds.r_id] = position
        for i, cloud in enumerate(read_kmer_clouds.kmers):
            self.coverage[i + position] += 1
            self.pos2kmers[i + position].append((read_kmer_clouds.r_id, i))
            self.clouds[i+position]
            for kmer in cloud:
                # if kmer in self.prohibited_kmers:
                #    continue
                if kmer not in self.kmer_positions:
                    self.kmer_positions[kmer] = i+position
                    self.clouds[i+position][kmer] = 1
                    self.kmer2reads[kmer].append(read_kmer_clouds.r_id)
                else:
                    kmer_pos = self.kmer_positions[kmer]
                    if kmer_pos == i+position:
                        self.clouds[kmer_pos][kmer] += 1
                        self.kmer2reads[kmer].append(read_kmer_clouds.r_id)
                    else:
                        # del self.clouds[kmer_pos][kmer]
                        # del self.kmer2reads[kmer]
                        # self.kmer_positions.pop(kmer, None)
                        # self.prohibited_kmers.add(kmer)
                        # new code after this line
                        if kmer in self.freq_kmers:
                            continue
                        del self.clouds[kmer_pos][kmer]
                        del self.kmer2reads[kmer]
                        self.kmer_positions.pop(kmer, None)

                        self.kmer_positions[kmer] = i+position
                        self.clouds[i+position][kmer] = 1
                        self.kmer2reads[kmer].append(read_kmer_clouds.r_id)
        self.update_max_pos()
        self.update_freq_clouds()

    def calc_rough_inters_score(self, read_kmer_cloud):
        return len(read_kmer_cloud.all_kmers & self.freq_kmers)

    def calc_inters_score(self, read_kmer_cloud,
                          min_position=None, max_position=None,
                          min_unit=2, min_inters=10,
                          verbose=False):
        best_score, best_pos = (0, 0), None
        kmers = read_kmer_cloud.kmers
        positions = list(self.clouds.keys())
        positions.sort()
        if min_position is not None:
            positions = [pos for pos in positions if pos >= min_position]
        if max_position is not None:
            positions = [pos for pos in positions if pos <= max_position]
        for pos in positions:
            score = [0, 0]
            max_i = min(self.max_pos-pos, len(kmers))
            for i in range(max_i):
                assert pos + i <= self.max_pos
                cloud = self.clouds[pos + i]
                freq_cloud = self.freq_clouds[pos + i]
                # inters = {k: cloud[k] for k in kmers[i] if cloud[k] >= self.min_cloud_kmer_freq}
                inters = freq_cloud & kmers[i]
                # inters_freqs = {k: cloud[k] for k in inters}
                if verbose:
                    print(pos, i, inters, len(inters))
                score[0] += len(inters) >= 1
                # score[1] += sum(inters.values())
                score[1] += len(inters)
                # score += len(inters)
                # score += sum(inters.values())
            if verbose:
                print(f'pos: {pos}, i: {i}, score: {score}')
            score = tuple(score)
            # we want to take the rightmost best, so >= instead of >
            if score[0] >= min_unit and score[1] >= min_inters and score >= best_score:
                best_score = score
                best_pos = pos
        return best_score, best_pos


def map_reads(cloud_contig, reads_kmer_clouds, threshold=(5, 10), verbose=False):
    scores = {}
    pos = {}
    for i, (r_id, kmer_clouds) in enumerate(reads_kmer_clouds.items()):
        max_pos = cloud_contig.max_pos - len(kmer_clouds.kmers) + 1
        best_score, best_pos = cloud_contig.calc_inters_score(kmer_clouds,
                                                              max_position=max_pos)
        if (best_pos == 0) or best_score > threshold:
            scores[r_id] = best_score
            pos[r_id] = best_pos
            if verbose:
                print(f"{i+1} / {len(reads_kmer_clouds)}, pos = {pos[r_id]}, score = {scores[r_id]}, id = {r_id}")
    return pos, scores
