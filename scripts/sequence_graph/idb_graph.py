# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter
import logging
import os

import networkx as nx

from sequence_graph.db_graph import DeBruijnGraph
from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.sequence_graph.idb_graph")


def def_get_min_mult(k, mode):
    if mode == 'assembly':
        return 1
    if mode == 'hifi':
        return 3
    assert mode == 'ont'
    if k < 100:
        return 25
    elif 100 <= k < 150:
        return 20
    elif 150 <= k < 200:
        return 15
    elif 200 <= k < 250:
        return 10
    elif 250 <= k < 300:
        return 5
    elif 300 <= k:
        return 3


def def_get_frequent_kmers(kmer_index, min_mult):
    return {kmer: cnt for kmer, cnt in kmer_index.items()
            if cnt >= min_mult}


def get_idb(string_set,
            mink, maxk,
            outdir,
            mode='ont',
            get_min_mult=None,
            get_frequent_kmers=None,
            all_kmer_index=None,
            step=1):

    assert mode in ['ont', 'hifi', 'assembly']
    if get_min_mult is None:
        get_min_mult = def_get_min_mult
    if get_frequent_kmers is None:
        get_frequent_kmers = def_get_frequent_kmers

    if all_kmer_index is None:
        all_kmer_index = string_set.get_kmer_index(mink=mink, maxk=maxk)
    else:
        assert all(k in all_kmer_index.keys()
                   for k in range(mink, maxk+1, step))

    smart_makedirs(outdir)
    dbs = {}
    all_frequent_kmers = {}

    contig_kmers = {}
    complex_kp1mers = {}
    for k in range(mink, maxk+1, step):
        min_mult = get_min_mult(k=k, mode=mode)
        kmer_index = all_kmer_index[k]
        frequent_kmers = get_frequent_kmers(kmer_index=kmer_index,
                                            min_mult=min_mult)
        # extending frequent kmers with contig kmers
        for kmer, cnt in contig_kmers.items():
            if kmer not in frequent_kmers:
                frequent_kmers[kmer] = cnt

        # extending frequent kmers with k+1-mers that pass through complex
        # nodes
        for kmer, cnt in complex_kp1mers.items():
            if kmer in frequent_kmers:
                assert cnt == frequent_kmers[kmer]
        frequent_kmers.update(complex_kp1mers)

        all_frequent_kmers[k] = frequent_kmers

        logger.info(f'k={k}')
        logger.info(f'#frequent kmers = {len(frequent_kmers)}')
        logger.info(f'min_mult = {min_mult}')

        db = DeBruijnGraph.from_kmers(kmers=frequent_kmers.keys(),
                                      kmer_coverages=frequent_kmers)

        ncc = nx.number_weakly_connected_components(db.nx_graph)
        logger.info(f'#cc = {ncc}')
        for cc in nx.weakly_connected_components(db.nx_graph):
            logger.info(len(cc))

        dot_file = os.path.join(outdir, f'db_k{k}.dot')
        db.write_dot(outfile=dot_file, export_pdf=False)

        dot_compact_file = os.path.join(outdir, f'db_k{k}_compact.dot')
        db.write_dot(outfile=dot_compact_file,
                     export_pdf=k==maxk,
                     compact=True)

        dbs[k] = db

        if k < maxk:
            contigs, _ = db.get_contigs()
            contig_kmers = Counter()
            for contig in contigs:
                for i in range(len(contig)-(k+1)+1):
                    kmer = contig[i:i+k+1]
                    contig_kmers[kmer] += 1

            complex_kp1mers = \
                db.get_paths_thru_complex_nodes(all_kmer_index[k+1])
    return dbs, all_frequent_kmers
