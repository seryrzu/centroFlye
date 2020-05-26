# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


import argparse
import logging
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

import networkx as nx
import pandas as pd
import seaborn as sns

from standard_logger import get_logger
from utils.bio import read_bio_seqs, calc_identity
from utils.git import get_git_revision_short_hash
from utils.os_utils import smart_makedirs
from utils.various import list2str


logger = logging.getLogger("centroFlye.monomers.monomer_db")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", required=True, help="Sequences")
    parser.add_argument("--outdir", required=True, )
    parser.add_argument("--max-ident", type=float, default=0.9)
    params = parser.parse_args()
    return params


def cluster_sequences(sequences, max_ident, logger=None):
    identities = {}
    graph = nx.Graph()
    rev_map = {seq: s_id for s_id, seq in sequences.items()}
    for s1 in rev_map:
        s_id1 = rev_map[s1]
        graph.add_node(s_id1)
        for s2 in rev_map:
            if s1 <= s2:
                continue
            s_id2 = rev_map[s2]
            ident = calc_identity(s1, s2)
            identities[(s_id1, s_id2)] = ident
            identities[(s_id2, s_id1)] = ident
            if ident > max_ident:
                graph.add_edge(s_id1, s_id2)
                logger.info(f'Ident {ident:0.2} b/w {s_id1} and {s_id2}')
    clusters = list(nx.connected_components(graph))
    logger.info(f'Extracted {len(clusters)} clusters')
    n_nontrivial_clusters = sum(len(cluster) > 1 for cluster in clusters)
    logger.info(f'Extracted {n_nontrivial_clusters} clusters of size > 1')
    return clusters, identities


def export_identities_heatmap(identities, outdir):
    sns.set(rc={'figure.figsize':(20,20)})
    identities = {k: round(ident*100, 2)
                  for k, ident in identities.items()}
    identities = pd.Series(list(identities.values()),
                           index=pd.MultiIndex.from_tuples(identities.keys()))
    identities = identities.unstack()
    heatmap = sns.heatmap(identities,
                          vmin=60, vmax=100,
                          annot=True,
                          cmap="YlGnBu").get_figure()
    heatmap_fn = os.path.join(outdir, 'identity_heatmap.pdf')
    heatmap.savefig(heatmap_fn, format='pdf')


def export_clusters(clusters, sequences, outdir, sep='\t', logger=None):
    outfile = os.path.join(outdir, 'clusters.tsv')
    with open(outfile, 'w') as f:
        header = ['cluster', 's_id', 'sequence']
        print(list2str(header), file=f)
        for i, cluster in enumerate(clusters):
            for s_id in cluster:
                outlist = [i, s_id, sequences[s_id]]
                print(list2str(outlist), file=f)


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfn = os.path.join(params.outdir, 'cluster_sequences.log')
    logger = get_logger(logfn,
                        logger_name='centroFlye: cluster sequences')

    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    sequences = read_bio_seqs(params.sequences)
    clusters, identities = \
        cluster_sequences(sequences, params.max_ident, logger=logger)
    export_identities_heatmap(identities, params.outdir)

    export_clusters(clusters, sequences, params.outdir, logger=logger)


if __name__ == "__main__":
    main()
