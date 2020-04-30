# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from monomers.monostring import MonoString
from monomers.monostring_correction import correct_monostrings

logger = logging.getLogger("centroFlye.monomers.monostring_set")


class MonoStringSet:
    def __init__(self, monostrings, monomer_db):
        self.monostrings = monostrings
        self.monomer_db = monomer_db

    @classmethod
    def from_sd_report(cls, report, sequences, monomer_db):
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
        raw_monostrings = get_raw_monostrings()
        monostrings = correct_monostrings(raw_monostrings)
        monostring_set = cls(monostrings=monostrings, monomer_db=monomer_db)
        return monostring_set
