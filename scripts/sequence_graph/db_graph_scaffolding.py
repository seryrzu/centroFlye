# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import defaultdict, Counter
import logging
import os
import statistics
from subprocess import check_call

import edlib
from joblib import Parallel, delayed
import networkx as nx
from tqdm import tqdm

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from config.config import config
from utils.bio import read_bio_seq, write_bio_seqs
from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.sequence_graph.db_graph_scaffolding")


def monoscaffolds2scaffolds(db, monoreads, outdir,
                            n_threads=config['common']['threads']):
    mappings = db.map_strings(monoreads)
    scaffolds, edge_scaffolds = \
        scaffolding(db, mappings, outdir=os.path.join(outdir, 'scaffolding'))
    locations = map_monoreads2scaffolds(monoreads, scaffolds)
    covered_scaffolds = \
        cover_scaffolds_w_reads(locations, monoreads, scaffolds,
                                outdir=os.path.join(outdir, 'coverage'))
    read_pseudounits = extract_read_pseudounits(covered_scaffolds,
                                                scaffolds,
                                                monoreads,
                                                min_coverage=0)
    polish_scaffolds(scaffolds, read_pseudounits,
                     outdir=os.path.join(outdir, 'polishing'),
                     n_iter=config['polishing']['n_iter'],
                     n_threads=n_threads,
                     min_cov=config['polishing']['min_cov'],
                     flye_bin='flye')


def scaffolding(db, mappings,
                min_connections=config['scaffolding']['min_connections'],
                outdir=None):
    def find_connections():
        connections = defaultdict(lambda: defaultdict(int))
        for r_id, mapping in mappings.items():
            if mapping is None or not mapping.valid:
                continue
            path = mapping.epath
            mapped_edges = set(path)
            inters = mapped_edges & long_edges
            if len(inters) > 1:
                indexes = []
                for long_edge in inters:
                    index = path.index(long_edge)
                    indexes.append(index)
                indexes.sort()
                for i, j in zip(indexes[:-1], indexes[1:]):
                    long_edge_pair = (path[i], path[j])
                    connection = tuple(path[i:j+1])
                    connections[long_edge_pair][connection] += 1
        logger.info('Connections in scaffolding')
        for ee, ee_connections in connections.items():
            logger.info(f'Pair of long edges = {ee}')
            for connection, val in ee_connections.items():
                logger.info(f'Connection = {connection},  #reads = {val}')
            logger.info('')
        return connections

    def build_scaffold_graph(connections):
        scaffold_graph = nx.DiGraph()
        for long_edge in long_edges:
            scaffold_graph.add_node(long_edge)

        for (e1, e2), connection_counts in connections.items():
            n_support_connections = sum(connection_counts.values())
            if n_support_connections >= min_connections:
                scaffold_graph.add_edge(e1, e2,
                                        connections=connection_counts)
        if outdir is not None:
            smart_makedirs(outdir)
            dotfile = os.path.join(outdir, 'scaffolding_graph.dot')
            nx.drawing.nx_pydot.write_dot(scaffold_graph, dotfile)
            pdffile = os.path.join(outdir, 'scaffolding_graph.pdf')
            check_call(['dot', '-Tpdf', dotfile, '-o', pdffile])
        return scaffold_graph

    def select_lists(scaffold_graph):
        longedge_scaffolds = []
        for cc in nx.weakly_connected_components(scaffold_graph):
            cc_sg = scaffold_graph.subgraph(cc)
            if nx.is_directed_acyclic_graph(cc_sg):
                longest_path = nx.dag_longest_path(cc_sg)
                longedge_scaffolds.append(longest_path)

        logger.info('Longedge scaffolds')
        for longedge_scaffold in longedge_scaffolds:
            logger.info(f'\t{longedge_scaffold}')
        return longedge_scaffolds

    def get_longest_extensions(longedge_scaffold):
        left_edge = longedge_scaffold[0]
        right_edge = longedge_scaffold[-1]
        lst_left_ext = []
        lst_right_ext = []
        for r_id, mapping in mappings.items():
            if mapping is None or not mapping.valid:
                continue
            path = mapping.epath
            try:
                left_index = path.index(left_edge)
                left_ext = path[:left_index]
                if len(left_ext) > len(lst_left_ext):
                    lst_left_ext = left_ext
            except ValueError:
                pass
            try:
                right_index = path.index(right_edge)
                right_ext = path[right_index+1:]
                if len(right_ext) > len(lst_right_ext):
                    lst_right_ext = right_ext
            except ValueError:
                pass
        return lst_left_ext, lst_right_ext

    def get_edge_scaffolds(longedge_scaffolds, connections):
        edge_scaffolds = []
        for longedge_scaffold in longedge_scaffolds:
            edge_scaffold = [longedge_scaffold[0]]
            for e1, e2 in zip(longedge_scaffold[:-1],
                              longedge_scaffold[1:]):
                edge_connect = connections[(e1, e2)]
                path = max(edge_connect, key=edge_connect.get)
                edge_scaffold += path[1:]
            left_ext, right_ext = get_longest_extensions(longedge_scaffold)
            edge_scaffold = left_ext + edge_scaffold + right_ext
            edge_scaffolds.append(edge_scaffold)

        logger.info(f'# edge scaffolds {len(edge_scaffolds)}:')
        for i, edge_scaffold in enumerate(edge_scaffolds):
            logger.info(f'{i}: {edge_scaffold}')

        if outdir is not None:
            scaffolds_fn = os.path.join(outdir, 'edge_scaffolds.txt')
            with open(scaffolds_fn, 'w') as f:
                for i, edge_scaffold in enumerate(edge_scaffolds):
                    print(f'{i}:\t{edge_scaffold}', file=f)

        return edge_scaffolds

    def get_sequence_scaffolds(edge_scaffolds):
        scaffolds = []
        for edge_scaffold in edge_scaffolds:
            scaffold = db.get_path(edge_scaffold)
            scaffolds.append(scaffold)
        return scaffolds

    long_edges = db.get_unique_edges(mappings=mappings)
    connections = find_connections()
    scaffold_graph = build_scaffold_graph(connections)
    longedge_scaffolds = select_lists(scaffold_graph)
    edge_scaffolds = get_edge_scaffolds(longedge_scaffolds,
                                        connections)

    scaffolds = get_sequence_scaffolds(edge_scaffolds)
    return scaffolds, edge_scaffolds


def map_monoreads2scaffolds(monoreads, scaffolds,
                            max_nloc=config['polishing']['max_nloc'],
                            gap_symb_matching=True):
    def map_monoread2scaffold(monoread, scaffold, add_matches):
        if gap_symb_matching:
            size = monoread.monomer_db.get_size()
            add_matches = [(monoread.gap_symb, i) for i in range(size)]
        else:
            add_matches = None
        align = edlib.align(monoread,
                            scaffold,
                            mode='HW',
                            task='path',
                            k=0,
                            additionalEqualities=add_matches)
        locs = align['locations']
        return locs

    all_locations = {}
    for i, scaffold in enumerate(scaffolds):
        all_locations[i] = {}
        for s_id, monoread in monoreads.items():
            locs = map_monoread2scaffold(monoread, scaffold, gap_symb_matching)
            if 0 < len(locs) <= max_nloc:
                all_locations[i][s_id] = locs
    return all_locations


def cover_scaffolds_w_reads(locations, monoreads, scaffolds, outdir=None):
    if outdir is not None:
        smart_makedirs(outdir)
    all_coverages = []
    n_scaffolds = len(scaffolds) - 1
    for s_i, scaffold in enumerate(scaffolds):
        coverage = [{} for i in range(len(scaffold) + 1)]
        for s_id, read_locs in locations[s_i].items():
            for st, en in read_locs:
                monoread = monoreads[s_id]
                if en - st + 1 != len(monoread):
                    continue
                for i, minst in enumerate(monoread.monoinstances):
                    p = st + i
                    coverage[p][minst.seq_id] = minst
        if outdir is not None:
            num_coverage = []
            for p in range(len(scaffold)):
                num_coverage.append(len(coverage[p]))
            plt.plot(num_coverage)
            plt.title(f'Coverage, scaffold {s_i} / {n_scaffolds}', fontsize=20)
            plt.xlabel('scaffold (monomers)', fontsize=18)
            plt.ylabel('coverage', fontsize=18)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig(os.path.join(outdir, f'scaffold_{s_i}.pdf'),
                        format='pdf')
            plt.close()

        all_coverages.append(coverage)
    return all_coverages


def partition_pseudounits(monostring):
    pseudounits = []
    i = 0
    while i < len(monostring):
        j = 0
        monomer_cnt = Counter()
        while i+j < len(monostring):
            monomer = monostring[i+j]
            monomer_cnt[monomer] += 1
            if monomer_cnt[monomer] > 1:
                break
            j += 1
        pseudounit = (i, i + j - 1)
        pseudounits.append(pseudounit)
        i += j
    return pseudounits


def extract_read_pseudounits(covered_scaffolds, scaffolds,
                             monostrings, min_coverage=0):
    all_read_pseudounits = []
    for s_i, scaffold in enumerate(scaffolds):
        scaf_pseudounits = partition_pseudounits(scaffold)
        covered_scaffold = covered_scaffolds[s_i]
        read_pseudounits = []
        for s, e in scaf_pseudounits:
            s_monomers = covered_scaffold[s]
            e_monomers = covered_scaffold[e]
            r_ids = set(s_monomers) & set(e_monomers)
            segms = {}
            for r_id in r_ids:
                st = s_monomers[r_id].st
                en = e_monomers[r_id].en
                segm = monostrings[r_id].nucl_sequence[st:en]
                segms[f'{r_id}|{st}|{en}'] = segm
            read_pseudounits.append((s, e, segms))
        all_read_pseudounits.append(read_pseudounits)
    return all_read_pseudounits


def polish_scaffolds(scaffolds, read_pseudounits, outdir,
                     n_iter=config['polishing']['n_iter'],
                     n_threads=config['common']['threads'],
                     min_cov=config['polishing']['min_cov'],
                     flye_bin='flye'):

    smart_makedirs(outdir)
    for i, scaffold in enumerate(scaffolds):

        scaf_outdir = os.path.join(outdir, f'scaffold_{i}')
        smart_makedirs(scaf_outdir)

        cmds = []
        for j, (s, e, pseudounits) in enumerate(read_pseudounits[i]):
            if len(pseudounits) == 0:
                continue

            pseudounit_outdir = os.path.join(scaf_outdir, f'pseudounit_{j}')
            smart_makedirs(pseudounit_outdir)

            reads_fn = os.path.join(pseudounit_outdir, 'reads.fasta')
            write_bio_seqs(reads_fn, pseudounits)

            template_fn = os.path.join(pseudounit_outdir, 'template.fasta')
            template_id, template = "", None
            r_units_lens = [len(read) for read in pseudounits.values()]
            med_len = statistics.median_high(r_units_lens)
            for r_id in sorted(pseudounits):
                read = pseudounits[r_id]
                if len(read) == med_len:
                    template_id = r_id
                    template = read
                    break
            assert len(pseudounits[template_id]) == med_len
            assert len(template) == med_len

            write_bio_seqs(template_fn,
                           {template_id: template})
            if len(pseudounits) <= min_cov:
                logger.info(f'Unit {j} will not be polished, '
                            f'using {template_id}')
                polished_pseudounit_fn = \
                    os.path.join(pseudounit_outdir, f'polished_{n_iter}.fasta')
                write_bio_seqs(polished_pseudounit_fn,
                               {template_id: template})
            else:
                cmd = [flye_bin,
                       '--nano-raw', reads_fn,
                       '--polish-target', template_fn,
                       '-i', n_iter,
                       '-t', 1,
                       '-o', pseudounit_outdir]
                cmd = [str(x) for x in cmd]
                cmds.append(cmd)
        Parallel(n_jobs=n_threads)(delayed(check_call)(cmd)
                                   for cmd in tqdm(cmds))
        polished_scaffold = []
        for j, (s, e, pseudounits) in enumerate(read_pseudounits[i]):
            pseudounit_outdir = os.path.join(scaf_outdir, f'pseudounit_{j}')
            polished_pseudounit_fn = \
                os.path.join(pseudounit_outdir, f'polished_{n_iter}.fasta')
            if os.path.isfile(polished_pseudounit_fn):
                polished_pseudounit = read_bio_seq(polished_pseudounit_fn)
                polished_scaffold.append(polished_pseudounit)

        polished_scaffold = ''.join(polished_scaffold)
        polished_scaffold_fn = os.path.join(scaf_outdir, f'scaffold_{i}.fasta')
        write_bio_seqs(polished_scaffold_fn,
                       {f'scaffold_{i}_niter_{n_iter}': polished_scaffold})
