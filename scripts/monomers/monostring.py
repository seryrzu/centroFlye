# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from enum import Enum
import logging

from utils.bio import RC

logger = logging.getLogger("centroFlye.monomers.monostring")


class Strand(Enum):
    FORWARD = '+'
    REVERSE = '-'


class Reliability(Enum):
    RELIABLE = '+'
    UNRELIABLE = '?'


class MonoInstance:
    def __init__(self, monomer_db, index, strand, seq_id, segment, st, en,
                 reliability):
        assert 0 <= index < monomer_db.get_total_monomers()
        self.monomer_db = monomer_db
        self.index = index
        self.strand = strand
        self.seq_id = seq_id
        self.segment = segment
        self.st = st
        self.en = en
        self.reliability = reliability

    def get_ref_seq(self):
        return self.monomer_db.seq[self.index]

    def get_name(self):
        return self.monomer_db.index2names[self.index]


class MonoString:
    # string is stored as a list because |monomer_db| often exceeds |ascii|
    gap_symb = '?'
    reverse_symb = "'"

    def __init__(self, seq_id, index2monoinstance, string, nucl_sequence):
        self.seq_id = seq_id
        self.index2monoinstance = index2monoinstance
        self.string = string
        self.nucl_sequence = nucl_sequence

    @classmethod
    def from_sd_record(cls, seq_id, monomer_db, sd_record, sequence,
                       max_gap=100):
        def id2index_strand(monomer_id, monomer_db=monomer_db):
            if monomer_id[-1] == cls.reverse_symb:
                monomer_id = monomer_id[:-1]
                strand = Strand.REVERSE
            else:
                strand = Strand.FORWARD
            index = monomer_db.id2index[monomer_id]
            return index, strand

        logger.debug(f'Constructing monostring for sequence {seq_id}')

        # Ignore the first and last monomer, as they may be incomplete
        # We remove last element now and later discard the first
        sd_record = sd_record[:-1]

        starts = sd_record.s_st.to_list()
        ends = [en + 1 for en in sd_record.s_en]

        # Now we discard the first monomer
        sd_record = sd_record[1:]
        ids = sd_record.monomer.to_list()
        reliabilities = [Reliability(raw_rel)
                         for raw_rel in sd_record.reliability]

        indexes_strands = map(id2index_strand, ids)
        indexes, strands = zip(*indexes_strands)

        coord1 = zip(starts[:-1], ends[:-1])
        coord2 = zip(starts[1:], ends[1:])
        mean_monomer_len = monomer_db.get_mean_monomer_len()
        n_monomers = monomer_db.get_total_monomers()

        index2monoinstance = {}
        string = []
        for (st1, en1), (st2, en2), id2, rel2, strand2, index2 in \
                zip(coord1, coord2, ids, reliabilities, strands, indexes):
            gap_len = st2 - en1
            if gap_len > max_gap:
                gaprun_len = int(round(gap_len / mean_monomer_len))
                string += cls.gap_symb * gaprun_len

            segment = sequence[st2:en2]
            if strand2 == Strand.REVERSE:
                segment = RC(segment)
            else:
                assert strand2 == Strand.FORWARD
            monoinstance = MonoInstance(monomer_db=monomer_db,
                                        index=index2,
                                        strand=strand2,
                                        seq_id=seq_id,
                                        segment=segment,
                                        st=st2,
                                        en=en2,
                                        reliability=rel2)
            cur_string_index = len(string)
            index2monoinstance[cur_string_index] = monoinstance
            if rel2 == Reliability.RELIABLE:
                string_monoindex = index2
                if strand2 == Strand.REVERSE:
                    # if strand is negative, index := index + |monomers|
                    string_monoindex += n_monomers
                string.append(string_monoindex)
            elif rel2 == Reliability.UNRELIABLE:
                string.append(cls.gap_symb)
            else:
                assert False

        monostring = cls(seq_id=seq_id,
                         index2monoinstance=index2monoinstance,
                         string=string,
                         nucl_sequence=sequence)
        return monostring
