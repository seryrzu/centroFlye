# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter, defaultdict
import logging

logger = logging.getLogger("centroFlye.utils.kmers")


def get_kmer_index_seq(seq, mink, maxk, ignored_chars=None, positions=False):
    if ignored_chars is None:
        ignored_chars = set()
    assert 0 < mink <= maxk
    kmerindex = {}
    for k in range(mink, maxk+1):
        if positions:
            kmerindex[k] = defaultdict(list)
        else:
            kmerindex[k] = Counter()
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            if len(set(kmer) & ignored_chars) == 0:
                kmer = tuple(kmer)
                if positions:
                    kmerindex[k][kmer].append(i)
                else:
                    kmerindex[k][kmer] += 1
    return kmerindex


def get_kmer_index(seqs, mink, maxk, ignored_chars=None, positions=False):
    if ignored_chars is None:
        ignored_chars = set()
    logger.info(f'Extracting kmer index mink={mink}, maxk={maxk}. '
                f'Positions={positions}')
    assert 0 < mink <= maxk
    if positions:
        kmer_index = {k: defaultdict(list) for k in range(mink, maxk+1)}
    else:
        kmer_index = {k: Counter() for k in range(mink, maxk+1)}
    nseqs = len(seqs)
    for i, (s_id, seq) in enumerate(seqs.items()):
        logger.debug(f'{i+1} / {nseqs}')
        ms_index = get_kmer_index_seq(seq, mink=mink, maxk=maxk,
                                      ignored_chars=ignored_chars,
                                      positions=positions)
        for k in range(mink, maxk+1):
            if positions:
                for kmer, pos in ms_index[k].items():
                    for p in pos:
                        kmer_index[k][kmer].append((s_id, p))
            else:
                for kmer, cnt in ms_index[k].items():
                    kmer_index[k][kmer] += cnt
    logger.info(f'Finished extracting kmer index')
    return kmer_index


def get_kmer_freqs_in_index(kmer_index, kmer_list):
    frequences = []
    for kmer in kmer_list:
        frequences.append(kmer_index[kmer])
    return frequences
