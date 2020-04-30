# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from enum import Enum
from itertools import count
import logging

import numpy as np

from utils.bio import RC

logger = logging.getLogger("centroFlye.monomers.monostring")


class Strand(Enum):
    FORWARD = '+'
    REVERSE = '-'

    @staticmethod
    def switch(strand):
        if strand == Strand.FORWARD:
            return Strand.REVERSE
        else:
            assert strand == Strand.REVERSE
            return Strand.FORWARD


class Reliability(Enum):
    RELIABLE = '+'
    UNRELIABLE = '?'


class MonoInstance:
    def __init__(self, monomer_db, index, strand, seq_id, nucl_segment,
                 st, en, seq_len,
                 reliability):
        assert 0 <= index < monomer_db.get_size()
        assert en - st == len(nucl_segment)
        self.monomer_db = monomer_db
        self.index = index
        self.strand = strand
        self.seq_id = seq_id
        self.nucl_segment = nucl_segment
        self.st = st
        self.en = en
        self.seq_len = seq_len
        self.reliability = reliability

    def get_ref_seq(self):
        return self.monomer_db.seq[self.index]

    def get_name(self):
        return self.monomer_db.index2names[self.index]

    def is_lowercase(self):
        return self.strand == Strand.REVERSE

    def reverse(self):
        self.nucl_segment = RC(self.nucl_segment)
        self.strand = Strand.switch(self.strand)
        # [st; en)
        self.st, self.en = self.seq_len - self.en, self.seq_len - self.st


class MonoString:
    # string is stored as a list because |monomer_db| often exceeds |ascii|
    gap_symb = '?'
    reverse_symb = "'"

    def __init__(self, seq_id, monoinstances, string, nucl_sequence,
                 monomer_db, is_reversed,
                 pref_cut, suf_cut):
        self.seq_id = seq_id
        self.monoinstances = monoinstances
        self.string = string
        self.nucl_sequence = nucl_sequence
        self.monomer_db = monomer_db
        self.is_reversed = is_reversed
        self.pref_cut = pref_cut
        self.suf_cut = suf_cut
        assert_monostring_validity(self)

    @classmethod
    def from_sd_record(cls, seq_id, monomer_db, sd_record, nucl_sequence):
        def trim_unreliable_monomers(sd_record):
            rels = [Reliability(raw_rel)
                    for raw_rel in sd_record.reliability]
            st = 0
            while st < len(rels) and rels[st] == Reliability.UNRELIABLE:
                st += 1
            en = len(rels) - 1
            while en >= 0 and rels[en] == Reliability.UNRELIABLE:
                en -= 1

            # Ignore the first and last monomer, as they may be incomplete
            len_rels = len(rels)
            st = max(st, 1)
            en = min(en, len(rels) - 2)

            # trim unreliable monomers
            sd_record = sd_record[st:en+1]
            rels = rels[st:en+1]

            pref_cut = st
            suf_cut = len_rels - 1 - en
            logger.debug(f'{pref_cut} prefix monomers are UNRELIABLE')
            logger.debug(f'{suf_cut} suffix monomers are UNRELIABLE')
            return sd_record, rels, pref_cut, suf_cut

        def get_monoinstances(sd_record, reliabilities,
                              pref_cut, suf_cut):
            def id2index_strand(monomer_id, monomer_db=monomer_db):
                if monomer_id[-1] == cls.reverse_symb:
                    monomer_id = monomer_id[:-1]
                    strand = Strand.REVERSE
                else:
                    strand = Strand.FORWARD
                index = monomer_db.id2index[monomer_id]
                return index, strand

            starts = sd_record.s_st.to_list()
            ends = [en + 1 for en in sd_record.s_en]

            ids = sd_record.monomer.to_list()
            indexes_strands = map(id2index_strand, ids)
            indexes, strands = zip(*indexes_strands)

            monoinstances = []
            for i, st, en, mono_id, rel, strand, index in \
                    zip(count(), starts, ends, ids,
                        reliabilities, strands, indexes):
                nucl_segment = nucl_sequence[st:en]
                monoinstance = MonoInstance(monomer_db=monomer_db,
                                            index=index,
                                            strand=strand,
                                            seq_id=seq_id,
                                            nucl_segment=nucl_segment,
                                            st=st,
                                            en=en,
                                            seq_len=len(nucl_sequence),
                                            reliability=rel)
                monoinstances.append(monoinstance)
            return monoinstances

        def reverse_if_needed(monoinstances, nucl_sequence,
                              max_lowercase=0.5):
            is_lowercase = [monoinstance.is_lowercase()
                            for monoinstance in monoinstances]
            perc_lower_case = np.mean(is_lowercase)
            is_reversed = perc_lower_case > max_lowercase
            if is_reversed:
                # reverse monoinstance
                monoinstances.reverse()
                for monoinstance in monoinstances:
                    monoinstance.reverse()
                # reverse nucl_sequence
                nucl_sequence = RC(nucl_sequence)
            return monoinstances, nucl_sequence, is_reversed

        def get_string(monoinstance):
            string = []
            for monoinstance in monoinstances:
                mono_index = monoinstance.index
                if monoinstance.reliability == Reliability.RELIABLE:
                    if monoinstance.strand == Strand.FORWARD:
                        string.append(mono_index)
                    else:
                        assert monoinstance.strand == Strand.REVERSE
                        string.append(mono_index + monomer_db.get_size())
                else:
                    assert monoinstance.reliability == Reliability.UNRELIABLE
                    string.append(cls.gap_symb)
            return string

        logger.debug(f'Constructing monostring for sequence {seq_id}')

        sd_record, reliabilities, pref_cut, suf_cut = \
            trim_unreliable_monomers(sd_record)

        monoinstances = get_monoinstances(sd_record=sd_record,
                                          reliabilities=reliabilities,
                                          pref_cut=pref_cut,
                                          suf_cut=suf_cut)

        monoinstances, nucl_sequence, is_reversed = \
            reverse_if_needed(monoinstances, nucl_sequence)
        string = get_string(monoinstances)

        logger.debug(f'Finished constructing monostring for sequence {seq_id}')
        logger.debug(f'    length of string = {len(string)}')
        logger.debug(f'    string: {string}')

        monostring = cls(seq_id=seq_id,
                         monoinstances=monoinstances,
                         string=string,
                         nucl_sequence=nucl_sequence,
                         monomer_db=monomer_db,
                         is_reversed=is_reversed,
                         pref_cut=pref_cut,
                         suf_cut=suf_cut)
        return monostring

    def __len__(self):
        return self.string.__len__()

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.string[sub.start:sub.stop:sub.step]
            return sublist
        return self.string[sub]


def assert_monostring_validity(monostring):
    string = monostring.string
    monomer_db = monostring.monomer_db
    monomer_db_size = monomer_db.get_size()
    monoinsts = monostring.monoinstances
    for i, monoinstance in enumerate(monoinsts):
        mono_index = monoinstance.index
        if monoinstance.strand == Strand.REVERSE:
            mono_index += monomer_db_size
        if monoinstance.reliability == Reliability.RELIABLE:
            assert mono_index == string[i]

    nucl_sequence = monostring.nucl_sequence
    for mi in monoinsts:
        if nucl_sequence[mi.st:mi.en] != mi.nucl_segment:
            print(mi.st, mi.en)
        assert nucl_sequence[mi.st:mi.en] == mi.nucl_segment
