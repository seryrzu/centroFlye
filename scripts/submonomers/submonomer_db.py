# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from submonomers.submonomers_extraction import submonomer_db_extraction

logger = logging.getLogger("centroFlye.submonomers.submonomer_db")


class Submonomer:
    def __init__(self, submono_index, monomer, seq):
        self.submono_index = submono_index
        self.monomer = monomer
        self.seq = seq

    def get_monomerid(self):
        return self.monomer.monomer_id

    def get_monoindex(self):
        return self.monomer.mono_index


class SubmonomerDB:
    def __init__(self, submonomers, monomer_db):
        self.submonomers = submonomers
        self.monomer_db = monomer_db

        for i, sm in enumerate(self.submonomers):
            assert self.submonomers[i].submono_index == i

    @classmethod
    def from_monostring_set(cls, monostring_set, coverage):
        return submonomer_db_extraction(monostring_set=monostring_set,
                                        coverage=coverage)

    def get_submonomers_by_mono_index(self, mono_index):
        selected_submonomers = []
        for submonomer in self.submonomers:
            if submonomer.get_monoindex() == mono_index:
                selected_submonomers.append(submonomer)
        return selected_submonomers

    def get_submonomers_by_monomer_id(self, monomer_id):
        mono_index = self.monomer_db.id2index[monomer_id]
        return self.get_submonomers_by_mono_index(mono_index)
