# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from joblib import Parallel, delayed
from submonomers.submonostring import SubmonoString, CorrectedSubmonoString
from utils.kmers import get_kmer_index
from utils.various import fst_iterable

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

    def __getitem__(self, sub):
        return self.submonostrings[sub]

    def to_tsv(self, filename, sep='\t', print_header=True):
        with open(filename, mode='w') as f:
            if print_header:
                header = ['seq_id', 'is_seq_reversed', 'submonochar',
                          'st', 'en', 'monoindex', 'monomer_id',
                          'mono_strand', 'mono_reliability',
                          'submono_is_identified', 'submono_is_unequivocal',
                          'nucl_segment']
                print(sep.join(header), file=f)

        for submonostring in self.submonostrings.values():
            submonostring.to_tsv(filename=filename,
                                 mode='a',
                                 sep=sep,
                                 print_header=False)

    def get_kmer_index(self, mink, maxk):
        assert len(self.submonostrings) > 0
        gap_symbs = \
            fst_iterable(self.submonostrings.values()).get_gap_symbols()
        return get_kmer_index(seqs=self.monostrings.values(),
                              mink=mink, maxk=maxk,
                              ignored_chars=gap_symbs)


class CorrectedSubmonoStringSet:
    def __init__(self, cor_submonostrings, submonomer_db):
        self.cor_submonostrings = cor_submonostrings
        self.submonomer_db = submonomer_db

    @classmethod
    def from_submonostring_set(cls, submonostring_set, mink=3, maxk=50):
        logger.info(f'Constructing kmer index for mink={mink}, maxk={maxk}')
        kmer_index = submonostring_set.get_kmer_index(mink=mink,
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

    def to_tsv(self, filename, sep='\t', print_header=True):
        # TODO remove this repetitive code. Compare to SubmonoStringSet
        with open(filename, mode='w') as f:
            if print_header:
                header = ['seq_id', 'is_seq_reversed', 'submonochar',
                          'st', 'en', 'monoindex', 'monomer_id',
                          'mono_strand', 'mono_reliability',
                          'submono_is_identified', 'submono_is_unequivocal',
                          'nucl_segment']
                print(sep.join(header), file=f)

        for cor_submonostring in self.cor_submonostrings.values():
            cor_submonostring.to_tsv(filename=filename,
                                     mode='a',
                                     sep=sep,
                                     print_header=False)
