# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from collections import defaultdict, Counter
import statistics
import subprocess

import edlib
from joblib import Parallel, delayed
import networkx as nx
import numpy as np
import os
from tqdm import tqdm

from sequence_graph.sequence_graph import SequenceGraph
from utils.bio import read_bio_seq, write_bio_seqs, parse_cigar
from utils.nx_all_simple_paths_multigraph import all_simple_edge_paths
from utils.os_utils import smart_makedirs
from utils.various import fst_iterable


logger = logging.getLogger("centroFlye.sequence_graph.db_graph")


class DeBruijnGraph(SequenceGraph):
    coverage = 'coverage'

    def __init__(self, nx_graph, nodeindex2label, nodelabel2index, k,
                 collapse=True):
        super().__init__(nx_graph=nx_graph,
                         nodeindex2label=nodeindex2label,
                         nodelabel2index=nodelabel2index,
                         collapse=collapse)
        self.k = k  # length of an edge in the uncompressed graph

    @classmethod
    def _generate_label(cls, par_dict):
        cov = par_dict[cls.coverage]
        length = par_dict[cls.length]

        mean_cov = np.mean(cov)
        label = f'len={length}\ncov={mean_cov:0.2f}'
        return label

    @classmethod
    def from_kmers(cls, kmers, kmer_coverages=None,
                   min_tip_cov=1, collapse=True):
        def add_kmer(kmer, coverage=1, color='black'):
            prefix, suffix = kmer[:-1], kmer[1:]

            if prefix in nodelabel2index:
                prefix_node_ind = nodelabel2index[prefix]
            else:
                prefix_node_ind = len(nodelabel2index)
                nodelabel2index[prefix] = prefix_node_ind
                nodeindex2label[prefix_node_ind] = prefix

            if suffix in nodelabel2index:
                suffix_node_ind = nodelabel2index[suffix]
            else:
                suffix_node_ind = len(nodelabel2index)
                nodelabel2index[suffix] = suffix_node_ind
                nodeindex2label[suffix_node_ind] = suffix

            length = 1
            coverage = [coverage]
            label = cls._generate_label({cls.length: length,
                                         cls.coverage: coverage})
            nx_graph.add_edge(prefix_node_ind, suffix_node_ind,
                              string=kmer,
                              length=length,
                              coverage=coverage,
                              label=label,
                              color=color)

        def remove_lowcov_tips():
            while True:
                edges_to_remove = []
                for s, e, key, data in nx_graph.edges(keys=True, data=True):
                    edge = (s, e, key)
                    cov = data[cls.coverage]
                    # We save coverage as list due to subsequent collapsing
                    # At this stage we did not collapse yet
                    assert len(cov) == 1
                    cov = cov[0]
                    indegree = nx_graph.in_degree(s)
                    outdegree = nx_graph.out_degree(e)
                    is_tip = (indegree == 0) or (outdegree == 0)
                    if is_tip and cov < min_tip_cov:
                        edges_to_remove.append((edge, indegree, outdegree))

                if len(edges_to_remove) == 0:
                    break

                for (s, e, key), indegree, outdegree in edges_to_remove:
                    nx_graph.remove_edge(s, e, key)

            isolates = list(nx.isolates(nx_graph))
            nx_graph.remove_nodes_from(isolates)
            for isolate in isolates:
                label = nodeindex2label[isolate]
                del nodeindex2label[isolate]
                del nodelabel2index[label]

        nx_graph = nx.MultiDiGraph()
        nodeindex2label = {}
        nodelabel2index = {}
        kmers = [tuple(kmer) for kmer in kmers]
        for kmer in kmers:
            if kmer_coverages is None:
                add_kmer(kmer)
            else:
                add_kmer(kmer, coverage=kmer_coverages[kmer])

        assert len(kmers)
        k = len(kmers[0])
        assert all(len(kmer) == k for kmer in kmers)

        remove_lowcov_tips()

        db_graph = cls(nx_graph=nx_graph,
                       nodeindex2label=nodeindex2label,
                       nodelabel2index=nodelabel2index,
                       k=k,
                       collapse=collapse)
        return db_graph

    def _add_edge(self, node, color, string,
                  in_node, out_node,
                  in_data, out_data,
                  edge_len):
        in_cov = in_data[self.coverage]
        out_cov = out_data[self.coverage]
        cov = sorted(in_cov + out_cov)
        assert len(cov) == edge_len
        label = self._generate_label({self.length: edge_len,
                                      self.coverage: np.mean(cov)})
        self.nx_graph.add_edge(in_node, out_node,
                               string=string,
                               coverage=cov,
                               label=label,
                               length=edge_len,
                               color=color)

    def get_complex_nodes(self):
        complex_nodes = []
        for node in self.nx_graph.nodes():
            indegree = self.nx_graph.in_degree(node)
            outdegree = self.nx_graph.out_degree(node)
            if indegree > 1 and outdegree > 1:
                complex_nodes.append(node)
        return complex_nodes

    def get_paths_thru_complex_nodes(self, kmer_index, min_mult=2):
        complex_nodes = self.get_complex_nodes()
        k = self.k

        selected_kp1mers = {}
        for node in complex_nodes:
            for in_edge in self.nx_graph.in_edges(node,
                                                  keys=True, data=True):
                for out_edge in self.nx_graph.out_edges(node,
                                                        keys=True, data=True):
                    in_kmer = in_edge[3][self.string][-k:]
                    out_kmer = out_edge[3][self.string][:k]

                    assert in_kmer[1:] == \
                        out_kmer[:-1] == \
                        self.nodeindex2label[node]

                    kp1 = in_kmer + (out_kmer[-1],)
                    if kmer_index[kp1] >= min_mult:
                        selected_kp1mers[kp1] = kmer_index[kp1]
        return selected_kp1mers

    def get_all_kmers(self):
        kmers = {}
        for s, e, key, data in self.nx_graph.edges(keys=True, data=True):
            coverage = data[self.coverage]
            string = data[self.string]
            assert len(coverage) == len(string)-self.k+1
            for i in range(len(string)-self.k+1):
                kmer = string[i:i+self.k]
                kmer = tuple(kmer)
                kmers[kmer] = coverage[i]
        return kmers

    def detect_bubbles(self, max_diff=1, cutoff=1):
        paths_set = set()
        for n1 in self.nx_graph.nodes:
            for n2 in self.nx_graph.nodes:
                if n1 == n2:
                    continue
                paths = all_simple_edge_paths(self.nx_graph, n1, n2,
                                              cutoff=cutoff)
                paths = list(paths)
                if n1 == 212 and n2 == 1346:
                    print(paths)
                if len(paths) > 1:
                    for path in paths:
                        path = tuple(path)
                        paths_set.add(path)

        all_lens2paths = defaultdict(lambda: defaultdict(list))
        for path in paths_set:
            fst_edge, lst_edge = path[0], path[-1]
            s = fst_edge[0]
            e = lst_edge[1]
            string = self.get_path(path)
            all_lens2paths[(s, e)][len(string)].append(path)

        # all_selected_paths = set()
        for s, e in list(all_lens2paths):
            filt_lens2paths = \
                {k: v for k, v in all_lens2paths[(s, e)].items()
                 if len(v) == 2}
            if len(filt_lens2paths):
                all_lens2paths[(s, e)] = filt_lens2paths
                # for path in filt_lens2paths.values():
                #     path = tuple(path)
                #     all_selected_paths.add(path)
            else:
                del all_lens2paths[(s, e)]

        path2strCov = defaultdict(list)
        for (s, e), len2paths in all_lens2paths.items():
            print(s, e)
            for length, paths in len2paths.items():
                for path in paths:
                    string = self.get_path(path)
                    path_coverages = []
                    for edge in path:
                        edge_data = self.nx_graph.get_edge_data(*edge)
                        coverage = edge_data[self.coverage]
                        coverage = np.mean(coverage)
                        path_coverages.append(coverage)
                    path_coverage = np.min(path_coverages)
                    print(path, path_coverage)
                    path2strCov[(s, e)].append((string, path_coverage))

        bubbles = []
        for (s, e), ((str1, cov1), (str2, cov2)) in path2strCov.items():
            diff = [i for i, (a, b) in enumerate(zip(str1, str2)) if a != b]
            print(s, e, len(diff), str1[diff[0]], str2[diff[0]])
            if len(diff) > max_diff:
                continue
            if cov1 >= cov2:
                bubbles.append((str2, str1))
                # self.nx_graph.remove_edge(s, e, 1)
            else:
                bubbles.append((str1, str2))
                # self.nx_graph.remove_edge(s, e, 0)
        # self.collapse_nonbranching_paths()
        return bubbles


def scaffolding(db, mappings, min_connections, long_edges):
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
        for ee, ee_connections in connections.items():
            print(ee)
            for connection, val in ee_connections.items():
                print(connection, val)
            print("\n")
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
        return scaffold_graph

    def select_lists(scaffold_graph):
        longedge_scaffolds = []
        for cc in nx.weakly_connected_components(scaffold_graph):
            cc_sg = scaffold_graph.subgraph(cc)
            if nx.is_directed_acyclic_graph(cc_sg):
                longest_path = nx.dag_longest_path(cc_sg)
                longedge_scaffolds.append(longest_path)
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
                # print(e1, e2)
                edge_connect = connections[(e1, e2)]
                path = max(edge_connect, key=edge_connect.get)
                edge_scaffold += path[1:]
            left_ext, right_ext = get_longest_extensions(longedge_scaffold)
            edge_scaffold = left_ext + edge_scaffold + right_ext
            edge_scaffolds.append(edge_scaffold)
        return edge_scaffolds

    def get_sequence_scaffolds(edge_scaffolds):
        scaffolds = []
        for edge_scaffold in edge_scaffolds:
            scaffold = db.get_path(edge_scaffold)
            scaffolds.append(scaffold)
        return scaffolds

    long_edges = set(long_edges)
    connections = find_connections()
    scaffold_graph = build_scaffold_graph(connections)
    longedge_scaffolds = select_lists(scaffold_graph)

    edge_scaffolds = get_edge_scaffolds(longedge_scaffolds,
                                        connections)

    scaffolds = get_sequence_scaffolds(edge_scaffolds)
    return scaffolds, edge_scaffolds


def map_monoreads2scaffolds(monoreads, scaffolds, max_nloc=1):
    def map_monoread2scaffold(monoread, scaffold):
        size = monoread.monomer_db.get_size()
        add_matches = [(monoread.gap_symb, i) for i in range(size)]
        align = edlib.align(monoread,
                            scaffold,
                            mode='HW',
                            task='path',
                            k=0,
                            additionalEqualities=add_matches)
        locs = align['locations']
        # if len(locs) == 2 and locs[0][1] >= 11924 and locs[1][0] <= 14393:
            # f_s, f_e = locs[0]
            # s_s, s_e = locs[1]
            # f_s = monoassembly.monoinstances[f_s].st
            # f_e = monoassembly.monoinstances[f_e].en
            # s_s = monoassembly.monoinstances[s_s].st
            # s_e = monoassembly.monoinstances[s_e].en
            # r_s = monoread.monoinstances[0].st
            # r_e = monoread.monoinstances[-1].en
            # strand = '+' if not monoread.is_reversed else '-'
            # print(s_id, locs)
            # cigar, cnt, a1, a2 = parse_cigar(align['cigar'], monoread, scaffold[locs[0][0]:locs[0][1]+1])
            # for c1, c2 in zip(a1, a2):
            #     print(c1, c2)


        #     align = edlib.align(monoread,
        #                         scaffold,
        #                         mode='HW',
        #                         task='locations',
        #                         k=10,
        #                         additionalEqualities=add_matches)
        #     locs = align['locations']
        #     print(locs)
        #     for st, en in locs:
        #         diff = [(a, b) for a, b in zip(monoread, scaffold[st:en+1])
        #                 if a != b and a != '?']
        #         print(diff)
        return locs

    all_locations = {}
    for i, scaffold in enumerate(scaffolds):
        all_locations[i] = {}
        for s_id, monoread in monoreads.items():
            locs = map_monoread2scaffold(monoread, scaffold)
            if 0 < len(locs) <= max_nloc:
                all_locations[i][s_id] = locs
    return all_locations


def cover_scaffolds_w_reads(locations, monoreads, scaffolds):
    all_coverages = []
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


def polish(scaffolds, read_pseudounits, outdir, monomer_db,
           n_iter=2, n_threads=30, flye_bin='flye'):

    def get_template(raw_monostring, monomer_db):
        print(raw_monostring)
        template = [fst_iterable(monomer_db.get_seqs_by_index(mono_index))
                    for mono_index in raw_monostring]
        return ''.join(template)

    smart_makedirs(outdir)
    for i, scaffold in enumerate(scaffolds):

        scaf_outdir = os.path.join(outdir, f'scaffold_{i}')
        smart_makedirs(scaf_outdir)

        cmds = []
        for j, (s, e, pseudounits) in enumerate(read_pseudounits[i]):
            print(j, len(read_pseudounits[i]))
            pseudounit_outdir = os.path.join(scaf_outdir, f'pseudounit_{j}')
            smart_makedirs(pseudounit_outdir)

            reads_fn = os.path.join(pseudounit_outdir, 'reads.fasta')
            write_bio_seqs(reads_fn, pseudounits)

            template_fn = os.path.join(pseudounit_outdir, 'template.fasta')
            # template_id, template_read = "", None
            # r_units_lens = [len(read) for read in pseudounits.values()]
            # med_len = statistics.median_high(r_units_lens)
            # for r_id in sorted(pseudounits):
            #     read = pseudounits[r_id]
            #     if len(read) == med_len:
            #         template_id = r_id
            #         template_read = read
            #         break
            # assert len(pseudounits[template_id]) == med_len
            # assert len(template_read) == med_len
            template = get_template(scaffold[s:e+1], monomer_db)
            template_id = f'{s}|{e}|'
            template_id += '_'.join([str(m) for m in scaffold[s:e+1]])
            write_bio_seqs(template_fn,
                           {template_id: template})

            cmd = [flye_bin,
                   '--nano-raw', reads_fn,
                   '--polish-target', template_fn,
                   '-i', n_iter,
                   '-t', 1,
                   '-o', pseudounit_outdir]
            cmd = [str(x) for x in cmd]
            print(' '.join(cmd))
            # subprocess.check_call(cmd)
            cmds.append(cmd)
        Parallel(n_jobs=n_threads)(delayed(subprocess.check_call)(cmd)
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
