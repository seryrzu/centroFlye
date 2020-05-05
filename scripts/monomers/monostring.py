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
    def __init__(self, monomer, strand, seq_id, nucl_segment,
                 st, en, seq_len,
                 reliability):
        assert en - st == len(nucl_segment)
        self.monomer = monomer
        self.strand = strand
        self.seq_id = seq_id
        self.nucl_segment = nucl_segment
        self.st = st
        self.en = en
        self.seq_len = seq_len
        self.reliability = reliability

    def get_monoindex(self):
        return self.monomer.mono_index

    def get_ref_seq(self):
        return self.monomer.seq

    def get_monoid(self):
        return self.monomer.mono_id

    def is_lowercase(self):
        return self.strand == Strand.REVERSE

    def is_reliable(self):
        return self.reliability is Reliability.RELIABLE

    def reverse(self):
        self.nucl_segment = RC(self.nucl_segment)
        self.strand = Strand.switch(self.strand)
        # [st; en)
        self.st, self.en = self.seq_len - self.en, self.seq_len - self.st


class MonoString:
    # monostring is stored as a list because |monomer_db| often exceeds |ascii|
    gap_symb = '?'
    reverse_symb = "'"

    def __init__(self, seq_id, monoinstances, raw_monostring, nucl_sequence,
                 monomer_db, is_reversed):
        self.seq_id = seq_id
        self.monoinstances = monoinstances
        self.raw_monostring = raw_monostring
        self.nucl_sequence = nucl_sequence
        self.monomer_db = monomer_db
        self.is_reversed = is_reversed
        assert_monostring_validity(self)

    @classmethod
    def from_sd_record(cls, seq_id, monomer_db, sd_record, nucl_sequence):
        def get_monoinstances(sd_record):
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
            reliabilities = [Reliability(raw_rel)
                             for raw_rel in sd_record.reliability]

            ids = sd_record.monomer.to_list()
            indexes_strands = map(id2index_strand, ids)
            indexes, strands = zip(*indexes_strands)

            monoinstances = []
            for i, st, en, rel, strand, mono_index in \
                    zip(count(), starts, ends,
                        reliabilities, strands, indexes):
                monomer = monomer_db.monomers[mono_index]
                nucl_segment = nucl_sequence[st:en]
                monoinstance = MonoInstance(monomer=monomer,
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
                            for monoinstance in monoinstances
                            if monoinstance.is_reliable()]
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
                mono_index = monoinstance.get_monoindex()
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

        logger.debug(f'Constructing raw_monostring for sequence {seq_id}')

        # Trim first and last monomer because they are often unreliable
        sd_record = sd_record[1:-1]

        monoinstances = get_monoinstances(sd_record=sd_record)

        monoinstances, nucl_sequence, is_reversed = \
            reverse_if_needed(monoinstances, nucl_sequence)
        string = get_string(monoinstances)

        logger.debug(f'Finished constructing raw_monostring for sequence {seq_id}')
        logger.debug(f'    length of string = {len(string)}')
        logger.debug(f'    string: {string}')

        monostring = cls(seq_id=seq_id,
                         monoinstances=monoinstances,
                         raw_monostring=string,
                         nucl_sequence=nucl_sequence,
                         monomer_db=monomer_db,
                         is_reversed=is_reversed)
        return monostring

    def __len__(self):
        return self.raw_monostring.__len__()

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.raw_monostring[sub.start:sub.stop:sub.step]
            return sublist
        return self.raw_monostring[sub]

    def get_perc_reliable(self):
        is_reliable = [monoinstance.is_reliable()
                       for monoinstance in self.monoinstances]
        perc_reliable = np.mean(is_reliable)
        return perc_reliable

    def get_perc_unreliable(self):
        return 1 - self.get_perc_reliable()

    def get_perc_lowercase(self):
        is_lowercase = [monoinstance.is_lowercase()
                        for monoinstance in self.monoinstances
                        if monoinstance.is_reliable()]
        perc_lowercase = np.mean(is_lowercase)
        return perc_lowercase

    def get_perc_uppercase(self):
        return 1 - self.get_perc_lowercase()

    def classify_monomerinstances(self, only_reliable=True):
        monoindexes = self.monomer_db.get_monoindexes()
        monomerinstances_dict = {monoindex: [] for monoindex in monoindexes}
        for mi in self.monoinstances:
            if (not only_reliable) or (only_reliable and mi.is_reliable()):
                monoindex = mi.get_monoindex()
                monomerinstances_dict[monoindex].append(mi)
        return monomerinstances_dict

    def get_monomerinstances_by_monoindex(self, mono_index,
                                          only_reliable=True):
        monomerinstances_dict = \
            self.classify_monomerinstances(only_reliable=only_reliable)
        return monomerinstances_dict[mono_index]

    def get_nucl_segment(self, st, en):
        assert 0 <= st < en < len(self.nucl_sequence)
        return self.nucl_sequence[st:en]


def assert_monostring_validity(monostring):
    string = monostring.raw_monostring
    monomer_db = monostring.monomer_db
    monomer_db_size = monomer_db.get_size()
    monoinsts = monostring.monoinstances
    for i, monoinstance in enumerate(monoinsts):
        mono_index = monoinstance.get_monoindex()
        if monoinstance.strand == Strand.REVERSE:
            mono_index += monomer_db_size
        if monoinstance.reliability == Reliability.RELIABLE:
            assert mono_index == string[i]

    nucl_sequence = monostring.nucl_sequence
    for mi in monoinsts:
        if nucl_sequence[mi.st:mi.en] != mi.nucl_segment:
            print(mi.st, mi.en)
        assert nucl_sequence[mi.st:mi.en] == mi.nucl_segment
