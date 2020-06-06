# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter
import logging

import numpy as np

from monomers.monostring import MonoString
from sequence_graph.db_graph import DeBruijnGraph
from utils.kmers import get_kmer_index, def_get_min_mult, \
    def_get_frequent_kmers, correct_kmers
from utils.bio import compress_homopolymer
from utils.various import fst_iterable

logger = logging.getLogger("centroFlye.monomers.monostring_set")


class MonoStringSet:
    def __init__(self, monostrings, monostrings_filt_out, monomer_db, mode):
        self.monostrings = monostrings
        self.monostrings_filt_out = monostrings_filt_out
        self.monomer_db = monomer_db
        self.mode = mode
        self.hybrids_corrected = False
        self.get_stats()

    @classmethod
    def from_sd_report(cls, report, sequences, monomer_db, mode,
                       correct_hybrids=False):
        def get_raw_monostrings(report=report,
                                sequences=sequences,
                                monomer_db=monomer_db):
            monostrings = {}
            for seq_id, seq_record in report.groupby('s_id'):
                seq_record = seq_record.sort_values(by=['s_st'])
                nucl_sequence = sequences[seq_id]
                monostring = \
                    MonoString.from_sd_record(seq_id=seq_id,
                                              monomer_db=monomer_db,
                                              sd_record=seq_record,
                                              nucl_sequence=nucl_sequence)
                monostrings[seq_id] = monostring
            return monostrings

        def filter_monostrings(monostrings,
                               max_lowercase=0.1,
                               max_unreliable=0.4):
            good_monostrings, bad_monostrings = {}, {}
            for seq_id, ms in monostrings.items():
                if ms.get_perc_lowercase() < max_lowercase and \
                   ms.get_perc_unreliable() < max_unreliable:
                    good_monostrings[seq_id] = ms
                else:
                    bad_monostrings[seq_id] = ms
            return good_monostrings, bad_monostrings

        assert mode in ['ont', 'hifi', 'assembly']
        unfiltered_monostrings = get_raw_monostrings()
        monostrings, monostrings_filt_out = \
            filter_monostrings(unfiltered_monostrings)
        monostring_set = cls(monostrings=monostrings,
                             monostrings_filt_out=monostrings_filt_out,
                             monomer_db=monomer_db,
                             mode=mode)
        if correct_hybrids:
            monostring_set.correct_likely_hybrids()
        return monostring_set

    def get_nucl_seq(self, seq_id):
        if seq_id in self.monostrings:
            return self.monostrings[seq_id].nucl_sequence
        else:
            assert seq_id in self.monostrings_filt_out
            return self.monostrings_filt_out[seq_id].nucl_sequence

    def __getitem__(self, sub):
        return self.monostrings[sub]

    def classify_monomerinstances(self, only_reliable=True):
        monoindexes = self.monomer_db.get_monoindexes()
        all_monomerinstances_dict = \
            {monoindex: [] for monoindex in monoindexes}
        for ms in self.monostrings.values():
            mi_dict = \
                ms.classify_monomerinstances(only_reliable=only_reliable)
            for monoindex in monoindexes:
                all_monomerinstances_dict[monoindex] += mi_dict[monoindex]
        return all_monomerinstances_dict

    def get_monomerinstances_by_monoindex(self, mono_index,
                                          only_reliable=True):
        all_monomerinstances_dict = \
            self.classify_monomerinstances(only_reliable=only_reliable)
        return all_monomerinstances_dict[mono_index]

    def get_stats(self, return_stats=False):
        def get_ngap_symbols(monostrings, compr_hmp=False):
            cnt = 0
            for monostring in monostrings.values():
                string = monostring.raw_monostring
                if compr_hmp:
                    string = compress_homopolymer(string, return_list=True)
                char_counter = Counter(string)
                cnt += char_counter[monostring.gap_symb]
            return cnt

        stats = {}
        monostrings = self.monostrings
        monostrings_filt_out = self.monostrings_filt_out
        monostrings_lens = [len(monostr) for monostr in monostrings.values()]
        stats['nmonostrings'] = len(monostrings_lens)
        stats['nfiltered_out'] = len(monostrings_filt_out)
        stats['min_len'] = np.min(monostrings_lens)
        stats['max_len'] = np.max(monostrings_lens)
        stats['mean_len'] = np.mean(monostrings_lens)
        stats['tot_len'] = np.sum(monostrings_lens)
        stats['ngaps'] = get_ngap_symbols(monostrings)
        stats['pgaps'] = stats['ngaps'] / stats['tot_len']
        stats['ngap_runs'] = get_ngap_symbols(monostrings, compr_hmp=True)

        logger.info(f'# monostrings: {stats["nmonostrings"]}')
        logger.info(f'# filtered monostrings: {stats["nfiltered_out"]}')
        logger.info(f'Min length = {stats["min_len"]}')
        logger.info(f'Mean length = {stats["mean_len"]}')
        logger.info(f'Max length = {stats["max_len"]}')
        logger.info(f'Total length = {stats["tot_len"]}')

        logger.info(f'#(%) Gap symbols = {stats["ngaps"]} ({stats["pgaps"]})')
        logger.info(f'#Gap runs = {stats["ngap_runs"]}')
        if return_stats:
            return stats

    def get_kmer_index(self, mink, maxk, positions=False):
        assert len(self.monostrings) > 0
        gap_symb = fst_iterable(self.monostrings.values()).gap_symb
        return get_kmer_index(seqs=self.monostrings,
                              mink=mink, maxk=maxk,
                              ignored_chars=set([gap_symb]),
                              positions=positions)

    def keys(self):
        return self.monostrings.keys()

    def values(self):
        return self.monostrings.values()

    def items(self):
        return self.monostrings.items()

    def correct_likely_hybrids(self, k=100,
                               get_frequent_kmers=def_get_frequent_kmers,
                               get_min_mult=def_get_min_mult):
        if self.hybrids_corrected:
            return

        while True:
            kmer_index_w_pos = self.get_kmer_index(mink=k,
                                                   maxk=k,
                                                   positions=True)
            kmer_index_w_pos = kmer_index_w_pos[k]
            kmer_index_wo_pos = {kmer: len(pos)
                                for kmer, pos in kmer_index_w_pos.items()}

            min_mult = get_min_mult(k=k, mode=self.mode)
            frequent_kmers = get_frequent_kmers(kmer_index=kmer_index_wo_pos,
                                                string_set=self,
                                                min_mult=min_mult)

            db = DeBruijnGraph.from_kmers(kmers=frequent_kmers.keys(),
                                        kmer_coverages=frequent_kmers,
                                        min_tip_cov=min_mult)

            bubbles = db.detect_bubbles()

            if len(bubbles) == 0:
                break

            readpos_to_change = set()
            for str1, str2 in bubbles:
                diff = [i for i, (a, b) in enumerate(zip(str1, str2)) if a != b]
                for d in diff:
                    bottom = max(0, d-k)
                    top = 1 + min(d+k, len(str1)-k)
                    for p in range(bottom, top):
                        i = d-p
                        assert 0 <= i < k
                        kmer = tuple(str1[p:p+k])
                        assert kmer in kmer_index_w_pos
                        for (s_id, s) in kmer_index_w_pos[kmer]:
                            readpos_to_change.add((s_id, s+i, str2[d]))
                            assert self.monostrings[s_id][s+i] == str1[d]
            for s_id, p, mono_index in readpos_to_change:
                self.monostrings[s_id][p] = mono_index

        self.hybrids_corrected = True

    def correct_likely_hybrids_depr(self, k=101):
        # This function is not used
        kmer_index_w_pos = self.get_kmer_index(mink=k, maxk=k, positions=True)
        kmer_index_w_pos = kmer_index_w_pos[k]

        k2 = k//2
        while True:
            # kmer_index_wo_pos = {kmer: len(pos)
            #                      for kmer, pos in kmer_index_w_pos.items()}

            # # TODO add mode
            # min_mult = def_get_min_mult(k=k, mode='ont')
            # frequent_kmers = \
            #     def_get_frequent_kmers(kmer_index_wo_pos,
            #                            string_set=self,
            #                            min_mult=min_mult)
            # frequent_kmers_w_pos = {kmer: kmer_index_w_pos[kmer]
            #                         for kmer in frequent_kmers}
            top_subst = correct_kmers(kmer_index_w_pos=kmer_index_w_pos,
                                      string_set=self)
            if top_subst is None:
                break
            kmer, mono_index = top_subst
            print(kmer, mono_index)
            for s_id, pos in kmer_index_w_pos[kmer].copy():
                print(s_id, pos)
                monostring = self.monostrings[s_id]

                bottom = max(0, pos-k2)
                top = 1 + min(pos+k2, len(monostring)-k)
                for p in range(bottom, top):
                    i = pos-p+k2
                    p_kmer = tuple(monostring[p:p+k])
                    mut_kmer = list(p_kmer)
                    mut_kmer[i] = mono_index
                    mut_kmer = tuple(mut_kmer)
                    if (s_id, p) in kmer_index_w_pos[p_kmer]:
                        kmer_index_w_pos[mut_kmer].add((s_id, p))
                        kmer_index_w_pos[p_kmer].remove((s_id, p))
                to_remove = []
                for kmer in kmer_index_w_pos:
                    if len(kmer_index_w_pos[kmer]) == 0:
                        to_remove.append(kmer)
                for kmer in to_remove:
                    del kmer_index_w_pos[kmer]

                change_pos = pos + k2
                monostring.raw_monostring[change_pos] = mono_index

            # assert_kmer_index_w_pos = self.get_kmer_index(mink=k,
            #                                     maxk=k,
            #                                     positions=True)
            # assert_kmer_index_w_pos = assert_kmer_index_w_pos[k]
            # for kmer in kmer_index_w_pos:
            #     if kmer_index_w_pos[kmer] != assert_kmer_index_w_pos[kmer]:
            #         print(kmer_index_w_pos[kmer] - assert_kmer_index_w_pos[kmer])
            #         print(assert_kmer_index_w_pos[kmer] - kmer_index_w_pos[kmer])
            # assert kmer_index_w_pos == assert_kmer_index_w_pos
            print("")

        self.hybrids_corrected = True
