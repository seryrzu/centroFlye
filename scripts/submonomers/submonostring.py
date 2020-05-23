# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from collections import Counter

from copy import deepcopy
import edlib

from monomers.monostring import MonoInstance, MonoString
from utils.kmers import get_kmer_index_seq
from utils.various import list2str

logger = logging.getLogger("centroFlye.submonomers.submonostring")


class SubmonoInstance(MonoInstance):
    # submonomers can be empty if it is unidentified
    def __init__(self, submonomers, dist2submonomers,
                 monomer, strand, seq_id, nucl_segment,
                 st, en, seq_len,
                 reliability):
        super().__init__(monomer, strand, seq_id, nucl_segment,
                         st, en, seq_len, reliability)
        self.submonomers = submonomers
        self.dist2submonomers = dist2submonomers

    @classmethod
    def from_monoinstance(cls, monoinstance, submonomers, dist2submonomers):
        return cls(submonomers=submonomers,
                   dist2submonomers=dist2submonomers,
                   monomer=monoinstance.monomer,
                   strand=monoinstance.strand,
                   seq_id=monoinstance.seq_id,
                   nucl_segment=monoinstance.nucl_segment,
                   st=monoinstance.st,
                   en=monoinstance.en,
                   seq_len=monoinstance.seq_len,
                   reliability=monoinstance.reliability)

    def is_unequivocal(self):
        return len(self.submonomers) == 1

    def is_unidentified(self):
        return len(self.submonomers) == 0

    def is_identified(self):
        identified = not self.is_unidentified()
        if identified:
            assert self.is_reliable
        return identified

    def is_equivocal(self):
        # return True if len(self.submonomers) > 1
        return self.is_identified() and not self.is_unequivocal()

    def make_unequivocal(self, submonoindex):
        submonomers_dict = {submonomer.submonoindex: submonomer
                            for submonomer in self.submonomers}
        assert submonoindex in submonomers_dict.keys()
        self.submonomers = [submonomers_dict[submonoindex]]


class SubmonoString(MonoString):
    submonogap_symb = '¿'    # if monomer is identifed, but no submonomer
    submonoequiv_symb = '*'  # if several submonomers are identified

    def __init__(self, seq_id, monoinstances, raw_monostring, nucl_sequence,
                 monomer_db, is_reversed,
                 submonoinstances, raw_submonostring, submonomer_db):
        super().__init__(seq_id=seq_id,
                         monoinstances=monoinstances,
                         raw_monostring=raw_monostring,
                         nucl_sequence=nucl_sequence,
                         monomer_db=monomer_db,
                         is_reversed=is_reversed)
        assert len(submonoinstances) == len(raw_submonostring)
        assert len(raw_monostring) == len(raw_submonostring)
        self.submonoinstances = submonoinstances
        self.raw_submonostring = raw_submonostring
        self.submonomer_db = submonomer_db

    @classmethod
    def from_monostring(cls, monostring, submonomer_db,
                        max_k=5, ext=10, max_diff=3):
        raw_submonostring = []
        submonoinstances = []

        monogap_symb = MonoString.gap_symb
        submonogap_symb = SubmonoString.submonogap_symb
        submonoequiv_symb = SubmonoString.submonoequiv_symb

        for i, mi in enumerate(monostring.monoinstances):
            assert len(submonoinstances) == len(raw_submonostring) == i
            if not mi.is_reliable():
                raw_submonostring.append(monogap_symb)
                sminst = SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                           submonomers={},
                                                           dist2submonomers={})
                submonoinstances.append(sminst)
                continue

            st, en = mi.st - ext, mi.en + ext
            nucl_segment = monostring.get_nucl_segment(st, en)
            monoindex = mi.get_monoindex()
            submonomers = \
                submonomer_db.get_submonomers_by_mono_index(monoindex)

            dist_smis = []
            for sm in submonomers:
                alignment = edlib.align(sm.seq, nucl_segment,
                                        k=max_k, mode='HW')
                dist = alignment['editDistance']
                if dist == -1:
                    continue
                dist_smis.append((dist, sm.submono_index))

            dist_smis.sort()
            if len(dist_smis) == 0:
                raw_submonostring.append(submonogap_symb)
                sminst = SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                           submonomers={},
                                                           dist2submonomers={})
                submonoinstances.append(sminst)
                continue

            distances, smindexes = list(zip(*dist_smis))
            best_dist, best_smindex = distances[0], smindexes[0]
            if len(distances) == 1:
                raw_submonostring.append(best_smindex)
                sm = submonomer_db.submonomers[best_smindex]
                submonomers = {best_smindex: sm}
                d2s = {best_smindex: best_dist}
                sminst = \
                    SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                      submonomers=submonomers,
                                                      dist2submonomers=d2s)
                submonoinstances.append(sminst)
                continue

            sec_best_dist, sec_best_smindex = distances[1], smindexes[1]
            # len(distances) > 1
            if best_dist == sec_best_dist == 0:
                best_sm = submonomer_db.submonomers[best_smindex]
                sec_best_sm = submonomer_db.submonomers[sec_best_smindex]
                logger.error(f'distance 0 to at least 2 segments:\n'
                             f'seq id = {monostring.seq_id}\n'
                             f'index in the string = {i}\n'
                             f'monoindex = {monoindex}\n'
                             f'nucl segment = {nucl_segment}\n'
                             f'best submonomer = {best_smindex}, '
                             f'{best_sm.seq}'
                             f'sec best submonomer = {sec_best_smindex}',
                             f'{sec_best_sm.seq}')
            if best_dist == 0 and sec_best_dist > 0:
                raw_submonostring.append(best_smindex)
                sm = submonomer_db.submonomers[best_smindex]
                submonomers = {best_smindex: sm}
                # distance to submonomer
                d2s = {best_smindex: best_dist}  # best_dist == 0
                sminst = \
                    SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                      submonomers=submonomers,
                                                      dist2submonomers=d2s)
                submonoinstances.append(sminst)
            else:
                submonomers = {}
                d2s = {}  # distance to submonomer
                for distance, smindex in zip(distances, smindexes):
                    if distance - best_dist > max_diff:
                        # distances is sorted
                        break
                    d2s[smindex] = distance
                    sm = submonomer_db.submonomers[smindex]
                    submonomers[smindex] = sm
                raw_submonostring.append(submonoequiv_symb)
                sminst = \
                    SubmonoInstance.from_monoinstance(monoinstance=mi,
                                                      submonomers=submonomers,
                                                      dist2submonomers=d2s)
                submonoinstances.append(sminst)
        submonostring = cls(seq_id=monostring.seq_id,
                            monoinstances=monostring.monoinstances,
                            raw_monostring=monostring.raw_monostring,
                            nucl_sequence=monostring.nucl_sequence,
                            monomer_db=monostring.monomer_db,
                            is_reversed=monostring.is_reversed,
                            submonoinstances=submonoinstances,
                            raw_submonostring=raw_submonostring,
                            submonomer_db=submonomer_db)
        return submonostring

    def get_gap_symbols(self):
        return {self.gap_symb, self.submonogap_symb, self.submonoequiv_symb}

    def get_kmer_index(self, mink, maxk):
        return get_kmer_index_seq(seq=self.raw_submonostring,
                                  mink=mink, maxk=maxk,
                                  ignored_chars=self.get_gap_symbols())

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.raw_submonostring[sub.start:sub.stop:sub.step]
            return sublist
        return self.raw_submonostring[sub]

    def to_tsv(self, filename, mode='w', sep='\t', print_header=True):
        seq_id = self.seq_id
        is_reversed = self.is_reversed
        with open(filename, mode=mode) as f:
            if print_header:
                header = ['seq_id', 'is_seq_reversed', 'submonochar',
                          'st', 'en', 'monoindex', 'monomer_id',
                          'mono_strand', 'mono_reliability',
                          'submono_is_identified', 'submono_is_unequivocal',
                          'nucl_segment']
                print(sep.join(header), file=f)

            for i, c in enumerate(self.raw_submonostring):
                sminst = self.submonoinstances[i]

                line = []
                line.append(seq_id)
                line.append(is_reversed)

                line.append(c)

                st, en = sminst.st, sminst.en
                line.append(st)
                line.append(en)

                monoindex = sminst.get_monoindex()
                line.append(monoindex)
                monomer_id = sminst.get_monoid()
                line.append(monomer_id)

                strand = sminst.strand
                line.append(strand.value)

                reliability = sminst.reliability
                line.append(reliability.value)

                sm_is_identified = sminst.is_identified()
                line.append(sm_is_identified)

                sm_is_unequivocal = sminst.is_unequivocal()
                line.append(sm_is_unequivocal)

                nucl_segment = sminst.nucl_segment
                line.append(nucl_segment)

                line = list2str(line, sep=sep)
                print(line, file=f)


class CorrectedSubmonoString(SubmonoString):
    def __init__(self, seq_id, monoinstances, raw_monostring, nucl_sequence,
                 monomer_db, is_reversed,
                 submonoinstances, raw_submonostring, submonomer_db,
                 cor_pos2submonoindex):
        super().__init__(seq_id=seq_id,
                         monoinstances=monoinstances,
                         raw_monostring=raw_monostring,
                         nucl_sequence=nucl_sequence,
                         monomer_db=monomer_db,
                         is_reversed=is_reversed,
                         submonoinstances=submonoinstances,
                         raw_submonostring=raw_submonostring,
                         submonomer_db=submonomer_db)
        self.cor_pos2submonoindex = cor_pos2submonoindex

        for pos, smindex in cor_pos2submonoindex.items():
            assert smindex == raw_submonostring[pos]

    @classmethod
    def from_submonostring(cls, submonostring, kmer_index, min_score=3):

        def correct_with_k(cor_submonostring, kmer_index_k, k):

            def gap_pos(string, gap_symb):
                return [i for i, c in enumerate(string)
                        if c == submonostring.submonoequiv_symb]

            def generate_all_closest_kmers(kmer, pos, smindexes):
                for smindex in smindexes:
                    cor_kmer = kmer
                    cor_kmer[pos] = smindex
                    yield tuple(cor_kmer)

            qpos = gap_pos(cor_submonostring, equiv_symb)
            pos2cor = {}
            while True:
                qpos = gap_pos(cor_submonostring, equiv_symb)
                if len(qpos) == 0:
                    break
                nonclumpedkmer = {}
                for i in qpos:
                    for j in range(max(0, i-k+1),
                                   min(i, len(cor_submonostring)-k+1)):
                        kmer = cor_submonostring[j:j+k]
                        assert len(kmer) == k
                        assert kmer[i-j] == equiv_symb
                        nogap = gap_symb not in kmer
                        oneequiv = Counter(kmer)[equiv_symb] == 1
                        nosubgap = subgap_symb not in kmer
                        if nogap and oneequiv and nosubgap:
                            nonclumpedkmer[i] = j
                            break
                cnt = 0
                for i, j in nonclumpedkmer.items():
                    kmer = cor_submonostring[j:j+k]
                    smindexes = \
                        submonostring.submonoinstances[i].dist2submonomers.keys()
                    score_kmers = []
                    for cor_kmer in \
                            generate_all_closest_kmers(kmer=kmer,
                                                       pos=i-j,
                                                       smindexes=smindexes):
                        if kmer_index_k[cor_kmer] > 0:
                            score = kmer_index_k[cor_kmer]
                            score_kmers.append((score, cor_kmer))
                    if len(score_kmers) == 1:
                        score, cor_kmer = score_kmers[0]
                        if score >= min_score:
                            cor_submonostring[i] = cor_kmer[i-j]
                            pos2cor[i] = cor_kmer[i-j]
                            cnt += 1
                if cnt == 0:
                    break
            return cor_submonostring, pos2cor

        gap_symb = submonostring.gap_symb
        subgap_symb = submonostring.submonogap_symb
        equiv_symb = submonostring.submonoequiv_symb

        cor_submonostring = deepcopy(submonostring.raw_submonostring)
        cor_pos2submonoindex = {}

        mink = min(kmer_index.keys())
        maxk = max(kmer_index.keys())

        for k in range(maxk, mink-1, -1):
            cor_submonostring, pos2cor = correct_with_k(cor_submonostring,
                                                        kmer_index[k],
                                                        k=k)
            cor_pos2submonoindex.update(pos2cor)

        submonoinstances = submonostring.submonoinstances
        cor_submonostring = cls(seq_id=submonostring.seq_id,
                                monoinstances=submonostring.monoinstances,
                                raw_monostring=submonostring.raw_monostring,
                                nucl_sequence=submonostring.nucl_sequence,
                                monomer_db=submonostring.monomer_db,
                                is_reversed=submonostring.is_reversed,
                                submonoinstances=submonoinstances,
                                raw_submonostring=cor_submonostring,
                                submonomer_db=submonostring.submonomer_db,
                                cor_pos2submonoindex=cor_pos2submonoindex)
        return cor_submonostring
