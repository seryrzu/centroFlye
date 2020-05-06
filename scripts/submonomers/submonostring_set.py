# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from collections import Counter

from joblib import Parallel, delayed
from submonomers.submonostring import SubmonoString

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
