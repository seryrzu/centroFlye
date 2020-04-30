# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

import numpy as np

from utils.bio import read_bio_seqs
from utils.os_utils import expandpath

logger = logging.getLogger("centroFlye.monomers.monomer_db")


class Monomer:
    def __init__(self, monomer_id, mono_index, seq):
        self.monomer_id = monomer_id
        self.mono_index = mono_index
        self.seq = seq


class MonomerDB:
    def __init__(self, id2index, index2id, monomers):
        self.id2index = id2index
        self.index2id = index2id
        self.monomers = monomers

    @classmethod
    def from_fasta_file(cls, fn):
        fn = expandpath(fn)
        logger.info(f'Creating Monomer DataBase from {fn}')
        raw_monomers = read_bio_seqs(fn)
        id2index = {}
        index2id = {}
        monomers = []
        for i, (monomer_id, monomer_seq) in enumerate(raw_monomers.items()):
            monomer = Monomer(monomer_id=monomer_id,
                              mono_index=i,
                              seq=monomer_seq)
            monomers.append(monomer)
            id2index[monomer_id] = i
            index2id[i] = monomer_id
            logger.debug(f'Monomer: index = {i}, id = {monomer_id}')
            logger.debug(f'         monomer sequence = {monomer_seq}')

        monomer_db = cls(id2index=id2index,
                         index2id=index2id,
                         monomers=monomers)
        logger.info(f'Finished Creating Monomer DataBase')
        return monomer_db

    def get_monomer_by_id(self, mono_id):
        return self.monomers[self.id2index[mono_id]]

    def get_seq_by_index(self, mono_index):
        assert 0 <= mono_index < len(self.monomers)
        return self.monomers[mono_index].seq

    def get_seq_by_id(self, monomer_id):
        mono_index = self.id2index[monomer_id]
        return self.monomers[mono_index].seq

    def get_ids(self):
        ids_generator = self.id2index.keys()
        return list(ids_generator)

    def get_size(self):
        return len(self.monomers)
