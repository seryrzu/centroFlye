# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter
import logging

logger = logging.getLogger("centroFlye.utils.kmers")


def get_kmer_index_seq(seq, mink, maxk, ignored_chars=None):
    if ignored_chars is None:
        ignored_chars = set()
    assert 0 < mink <= maxk
    kmerindex = {}
    for k in range(mink, maxk+1):
        kmerindex[k] = Counter()
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            if len(set(kmer) & ignored_chars) == 0:
                kmer = tuple(kmer)
                kmerindex[k][kmer] += 1
    return kmerindex


def get_kmer_index(seqs, mink, maxk, ignored_chars=None):
    if ignored_chars is None:
        ignored_chars = set()
    seqs = list(seqs)  # in case seqs is a iterable
    logger.info(f'Extracting kmer index mink={mink}, maxk={maxk}')
    assert 0 < mink <= maxk
    kmer_index = {k: Counter() for k in range(mink, maxk+1)}
    nseqs = len(seqs)
    for i, seq in enumerate(seqs):
        logger.debug(f'{i+1} / {nseqs}')
        ms_index = get_kmer_index_seq(seq, mink=mink, maxk=maxk,
                                      ignored_chars=ignored_chars)
        for k in range(mink, maxk+1):
            for kmer, cnt in ms_index[k].items():
                kmer_index[k][kmer] += cnt
    logger.info(f'Finished extracting kmer index')
    return kmer_index


def get_kmer_freqs_in_index(kmer_index, kmer_list):
    frequences = []
    for kmer in kmer_list:
        frequences.append(kmer_index[kmer])
    return frequences
