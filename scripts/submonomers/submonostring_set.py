# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from collections import Counter

from joblib import Parallel, delayed
from submonomers.submonostring import SubmonoString, CorrectedSubmonoString

logger = logging.getLogger('centroFlye.submonomers.submonostring_set')


class SubmonoStringSet:
    def __init__(self, submonostrings, submonomer_db):
        self.submonostrings = submonostrings
        self.submonomer_db = submonomer_db

    @classmethod
    def from_monostringset(cls, monostring_set, submonomer_db, n_jobs=30):
        submonostrings_list = Parallel(n_jobs=n_jobs, verbose=60,
                                       backend='multiprocessing')\
            (delayed(SubmonoString.from_monostring)
            (monostring, submonomer_db)
            for monostring in monostring_set.monostrings.values())
        submonostrings = {}
        for submonostring in submonostrings_list:
            seq_id = submonostring.seq_id
            submonostrings[seq_id] = submonostring
        submonostring_set = cls(submonostrings=submonostrings,
                                submonomer_db=submonomer_db)
        return submonostring_set

    def get_submonokmer_index(self, mink, maxk):
        assert 0 < mink <= maxk
        kmer_index = {}
        for k in range(mink, maxk+1):
            kmer_index[k] = Counter()
        for submonostring in self.submonostrings.values():
            sms_index = submonostring.get_submonokmer_index(mink=mink,
                                                            maxk=maxk)
            for k in range(mink, maxk+1):
                for kmer, cnt in sms_index[k].items():
                    kmer_index[k][kmer] += cnt
        return kmer_index

    def __getitem__(self, sub):
        return self.submonostrings[sub]


class CorrectedSubmonoStringSet:
    def __init__(self, cor_submonostrings, submonomer_db):
        self.cor_submonostrings = cor_submonostrings
        self.submonomer_db = submonomer_db

    @classmethod
    def from_submonostring_set(cls, submonostring_set, mink=3, maxk=50):
        logger.info(f'Constructing kmer index for mink={mink}, maxk={maxk}')
        kmer_index = submonostring_set.get_submonokmer_index(mink=mink,
                                                             maxk=maxk)
        logger.info('Finished constructing kmer index')
        logger.info('Correcting submonostrings')
        cor_submonostrings = {}
        for seq_id, submonostring in submonostring_set.submonostrings.items():
            cor_submonostring = \
                CorrectedSubmonoString.from_submonostring(submonostring,
                                                          kmer_index)
            cor_submonostrings[seq_id] = cor_submonostring
        logger.info('Finished correcting submonostrings')
        cor_submonostring_set = \
            cls(cor_submonostrings=cor_submonostrings,
                submonomer_db=submonostring_set.submonomer_db)
        return cor_submonostring_set

    def __getitem__(self, sub):
        return self.cor_submonostrings[sub]
