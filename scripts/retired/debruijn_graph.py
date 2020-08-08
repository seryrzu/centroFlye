# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import defaultdict, Counter
from copy import deepcopy
from itertools import groupby
import os
import subprocess
import statistics

import networkx as nx
import numpy as np

from utils.bio import read_bio_seq, read_bio_seqs, write_bio_seqs, RC
from utils.os_utils import smart_makedirs


class DeBruijnGraph:
    def __init__(self, k, max_uniq_cov=60, min_uniq_len=700):
        self.graph = nx.MultiDiGraph()
        self.k = k
        self.node_mapping = {}
        self.rev_node_mapping = {}
        self.max_uniq_cov = max_uniq_cov
        self.min_uniq_len = min_uniq_len

    def add_kmer(self, kmer, coverage=1):
        prefix, suffix = kmer[:-1], kmer[1:]

        if prefix in self.node_mapping:
            prefix_node_ind = self.node_mapping[prefix]
        else:
            prefix_node_ind = len(self.node_mapping)
            self.node_mapping[prefix] = prefix_node_ind
            self.rev_node_mapping[prefix_node_ind] = prefix

        if suffix in self.node_mapping:
            suffix_node_ind = self.node_mapping[suffix]
        else:
            suffix_node_ind = len(self.node_mapping)
            self.node_mapping[suffix] = suffix_node_ind
            self.rev_node_mapping[suffix_node_ind] = suffix

        self.graph.add_edge(
            prefix_node_ind,
            suffix_node_ind,
            edge_kmer=kmer,
            length=1,
            coverages=[coverage],
            label=f'len=1\ncov={coverage}',
            color='black')

    def add_kmers(self, kmers, coverage=None):
        for kmer in kmers:
            if coverage is None:
                self.add_kmer(kmer)
            else:
                self.add_kmer(kmer, coverage=coverage[kmer])

    def index_edges(self, min_k=2):
        try:
            return self.db_index
        except AttributeError:
            all_index = {}
            for k in range(min_k, self.k+1):
                index = defaultdict(list)
                for e_ind, edge in enumerate(self.graph.edges(keys=True)):
                    edge_seq = self.graph.get_edge_data(*edge)['edge_kmer']
                    for i in range(len(edge_seq)-k+1):
                        kmer = edge_seq[i:i+k]
                        index[kmer].append((e_ind, i))
                index_unique = {kmer: pos[0]
                                for kmer, pos in index.items()
                                if len(pos) == 1}
                all_index[k] = index_unique
            self.db_index = all_index
            return all_index

    def collapse_nonbranching_paths(self):
        def node_on_nonbranching_path(graph, node):
            in_degree = graph.in_degree(node)
            out_degree = graph.out_degree(node)
            if graph.in_degree(node) != 1 or \
                    graph.out_degree(node) != 1:
                return False
            in_edges = list(graph.in_edges(node))
            in_edge = in_edges[0]
            out_edges = list(graph.out_edges(node))
            out_edge = out_edges[0]
            return in_edge != out_edge

        for node in list(self.graph.nodes()):
            if node_on_nonbranching_path(self.graph, node):
                in_edge = list(self.graph.in_edges(node, keys=True))[0]
                out_edge = list(self.graph.out_edges(node, keys=True))[0]
                # in_edge_color = self.graph.edges[in_edge]['color']
                # out_edge_color = self.graph.edges[out_edge]['color']
                in_edge_kmer = self.graph.edges[in_edge]['edge_kmer']
                out_edge_kmer = self.graph.edges[out_edge]['edge_kmer']
                in_edge_cov = self.graph.edges[in_edge]['coverages']
                out_edge_cov = self.graph.edges[out_edge]['coverages']

                in_node = in_edge[0]
                out_node = out_edge[1]

                new_kmer = in_edge_kmer + \
                    out_edge_kmer[-(len(out_edge_kmer)-self.k+1):]
                new_coverages = sorted(in_edge_cov + out_edge_cov)
                new_edge_len = len(new_coverages)
                new_edge_med_cov = np.median(new_coverages)
                if new_edge_len + self.k - 1 >= self.min_uniq_len and \
                        new_edge_med_cov <= self.max_uniq_cov:
                    new_edge_color = 'blue'
                else:
                    new_edge_color = 'black'
                self.graph.add_edge(
                    in_node, out_node, edge_kmer=new_kmer,
                    coverages=new_coverages, length=new_edge_len,
                    label=f'len={new_edge_len}\ncov={new_edge_med_cov}',
                    color=new_edge_color)
                self.graph.remove_node(node)

    def get_edges(self):
        self.collapse_nonbranching_paths()
        contigs, coverages = [], []
        for u, v, keys in self.graph.edges(data=True):
            contigs.append(keys['edge_kmer'])
            coverages.append(np.median(keys['coverages']))
        return contigs, coverages

    def get_path(self, list_edges):
        path = []
        path += self.graph.get_edge_data(*list_edges[0])['edge_kmer']
        for edge in list_edges[1:]:
            new_edge_seq = self.graph.get_edge_data(*edge)['edge_kmer']
            assert tuple(path[-(self.k-1):]) == new_edge_seq[:self.k-1]
            path += new_edge_seq[self.k-1:]

        # process cycled path
        if list_edges[0][0] == list_edges[-1][1]:
            path = path[:-(self.k-1)]
        return tuple(path)

    def get_edgepath2coords(self, list_edges):
        edgepath2coords = {}
        str_coord = 0
        path = self.get_path(list_edges)
        for i, edge_id in enumerate(list_edges):
            edge_str = self.graph.get_edge_data(*edge_id)['edge_kmer']
            for j in range(len(edge_str)):
                assert path[str_coord] == edge_str[j]
                edgepath2coords[(i, j)] = str_coord
                str_coord += 1
            str_coord -= (self.k - 1)
            edgepath2coords[i] = str_coord
        return edgepath2coords

    def get_contigs(self):
        def get_longest_valid_outpaths(graph):
            def get_valid_outpath_edge(edge, taken_edges=None):
                if taken_edges is None:
                    taken_edges = set()
                path = [edge]
                out_node = edge[1]
                out_degree_out_node = graph.out_degree(out_node)
                if out_degree_out_node == 1:
                    out_edge = list(
                        graph.out_edges(
                            out_node, keys=True))[0]
                    if out_edge not in taken_edges:
                        taken_edges.add(edge)
                        out_edge_path = \
                            get_valid_outpath_edge(out_edge,
                                                   taken_edges=taken_edges)
                        path += out_edge_path
                return path

            outpaths = {}
            for edge in graph.edges(keys=True):
                if edge not in outpaths:
                    outpaths[edge] = get_valid_outpath_edge(edge)
                    for i, e in enumerate(outpaths[edge][1:]):
                        outpaths[e] = outpaths[edge][i+1:]
            return outpaths

        self.collapse_nonbranching_paths()
        outpaths = get_longest_valid_outpaths(self.graph)

        rev_graph = self.graph.reverse()
        rev_inpaths = get_longest_valid_outpaths(rev_graph)
        inpaths = {}
        for rev_edge in rev_inpaths:
            edge = (rev_edge[1], rev_edge[0], rev_edge[2])
            inpaths[edge] = rev_inpaths[rev_edge][::-1]
            inpaths[edge] = [(e[1], e[0], e[2]) for e in inpaths[edge]]

        valid_paths = []
        for edge in outpaths:
            valid_path = inpaths[edge]
            seen_edges = set(inpaths[edge])
            for e in outpaths[edge][1:]:
                if e not in seen_edges:
                    valid_path.append(e)
                    seen_edges.add(e)
                else:
                    break
            valid_path = tuple(valid_path)
            valid_paths.append(valid_path)

        valid_paths = list(set(valid_paths))

        # Select only paths that are not a subpath of other paths
        selected_paths = []
        for path1 in valid_paths:
            has_duplicate = False
            for path2 in valid_paths:
                if path1 == path2:
                    continue
                for i in range(len(path2)-len(path1)+1):
                    if path1 == path2[i:i+len(path1)]:
                        has_duplicate = True
                        break
                if has_duplicate:
                    break
            if not has_duplicate:
                selected_paths.append(path1)
        valid_paths = selected_paths

        contigs = []
        for path in valid_paths:
            contigs.append(self.get_path(path))
        contigs = list(set(contigs))
        return contigs, valid_paths

    def map_reads(self, monoreads, gap_symbs='?*¿', verbose=True):
        def split_list(lst, es):
            es = set(gap_symbs)
            spl_lst = []
            spl = []
            for x in lst:
                if x not in es:
                    spl.append(x)
                else:
                    spl_lst.append(spl)
                    spl = []
            if len(spl):
                spl_lst.append(spl)
            return spl_lst

        self.index_edges()
        mapping = {}
        db_edges = list(self.graph.edges(keys=True))
        # monoreads = {r_id: monoread.string
        #              for r_id, monoread in monoreads.items()}
        for i, (r_id, string) in enumerate(monoreads.items()):
            if verbose:
                print(i+1, len(monoreads))

            split_strings = split_list(string, gap_symbs)
            read_coords = []
            cumm_len = 0
            for split_ind, split_string in enumerate(split_strings):
                for i in range(len(split_string)-self.k+1):
                    kmer = tuple(split_string[i:i+self.k])
                    if kmer in self.db_index[len(kmer)]:
                        read_coords.append((self.db_index[len(kmer)][kmer],
                                            cumm_len + i))
                        string_kmer = string[cumm_len+i:cumm_len+i+min(self.k, len(split_string))]
                        string_kmer = tuple(string_kmer)
                        assert string_kmer == kmer
                cumm_len += len(split_string) + 1

            if len(read_coords) == 0:
                mapping[r_id] = None
                continue

            (edge_ind, edge_coord), _ = read_coords[0]
            edge_path = [edge_ind]
            for (cur_edge_ind, cur_edge_coord), _ in read_coords[1:]:
                if edge_ind != cur_edge_ind or cur_edge_coord <= edge_coord:
                    edge_path.append(cur_edge_ind)
                edge_ind = cur_edge_ind
                edge_coord = cur_edge_coord
            path = [db_edges[edge_ind] for edge_ind in edge_path]

            valid_path = True
            for e1, e2 in zip(path[:-1], path[1:]):
                if e1[1] != e2[0]:
                    valid_path = False
                    break
            mapping[r_id] = (read_coords[0], read_coords[-1],
                                valid_path, path)
        return mapping

    def get_long_edges(self):
        edges = {}
        for edge in self.graph.edges(data=True, keys=True):
            edge_color = edge[-1]['color']
            if edge_color == 'blue':
                edges[edge[:-1]] = edge[-1]['edge_kmer']
        return edges


def get_all_kmers(strings, k, gap_symbs='?*¿'):
    all_kmers = Counter()
    read_kmer_locations = defaultdict(list)
    for r_id, string in strings.items():
        for i in range(len(string)-k+1):
            kmer = tuple(string[i:i+k])
            if all(gap_symb not in kmer for gap_symb in gap_symbs):
                all_kmers[kmer] += 1
                read_kmer_locations[kmer].append((r_id, i))
    return all_kmers, read_kmer_locations


def get_complex_nodes(graph):
    complex_nodes = []
    for node in graph.nodes():
        indegree, outdegree = graph.in_degree(node), graph.out_degree(node)
        if indegree > 1 and outdegree > 1:
            complex_nodes.append(node)
    return complex_nodes


def get_paths_thru_complex_nodes(db, strings, min_mult=2):
    complex_nodes = get_complex_nodes(db.graph)
    k = db.k
    all_kp1mers, _ = get_all_kmers(strings, k=k+1)
    selected_kp1mers = {}
    for node in complex_nodes:
        for in_edge in db.graph.in_edges(node, keys=True, data=True):
            for out_edge in db.graph.out_edges(node, keys=True, data=True):
                in_kmer = in_edge[3]['edge_kmer'][-k:]
                out_kmer = out_edge[3]['edge_kmer'][:k]
                assert in_kmer[1:] == out_kmer[:-1]
                kp1 = in_kmer + tuple([out_kmer[-1]])
                if all_kp1mers[kp1] >= min_mult:
                    selected_kp1mers[kp1] = all_kp1mers[kp1]
    return selected_kp1mers


def get_frequent_kmers(strings, k, min_mult=5):
    all_kmers, read_kmer_locations = get_all_kmers(strings, k)
    frequent_kmers = {kmer: cnt for kmer, cnt in all_kmers.items()
                      if cnt >= min_mult}
    frequent_kmers_read_pos = {
        kmer: read_kmer_locations[kmer] for kmer in frequent_kmers}
    return frequent_kmers, frequent_kmers_read_pos


def iterative_graph(monostrings, min_k, max_k, outdir,
                    def_min_mult=None, step=1, starting_graph=None,
                    verbose=True):
    def get_min_mult(k):
        if def_min_mult is not None:
            return def_min_mult
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
    smart_makedirs(outdir)
    dbs, uncompr_dbs, all_contigs = {}, {}, {}
    all_frequent_kmers, all_frequent_kmers_read_pos = {}, {}
    strings = {r_id: monostring for r_id, monostring in monostrings.items()}
    input_strings = strings.copy()
    complex_kp1mers = {}

    min_mult = get_min_mult(min_k)
    if starting_graph is not None:
        contigs, contig_paths = starting_graph.get_contigs()
        for i in range(len(contigs)):
            for j in range(min_mult):
                input_strings[f'contig_k{min_k}_i{i}_j{j}'] = contigs[i]

        complex_kp1mers = get_paths_thru_complex_nodes(starting_graph, strings)

    for k in range(min_k, max_k+1, step):
        min_mult = get_min_mult(k)
        frequent_kmers, frequent_kmers_read_pos = \
            get_frequent_kmers(input_strings, k=k, min_mult=min_mult)
        frequent_kmers.update(complex_kp1mers)
        if verbose:
            print(f'\nk={k}')
            print(f'#frequent kmers = {len(frequent_kmers)}')
            print(f'min_mult = {min_mult}')
        all_frequent_kmers[k] = frequent_kmers
        all_frequent_kmers_read_pos[k] = frequent_kmers_read_pos

        db = DeBruijnGraph(k=k)
        db.add_kmers(frequent_kmers, coverage=frequent_kmers)
        uncompr_dbs[k] = deepcopy(db)

        db.collapse_nonbranching_paths()
        if verbose and nx.number_weakly_connected_components(db.graph) > 1:
            print(f'#cc = {nx.number_weakly_connected_components(db.graph)}')
            for cc in nx.weakly_connected_components(db.graph):
                print(len(cc))
            # break
        dbs[k] = db

        dot_file = os.path.join(outdir, f'db_k{k}.dot')
        nx.drawing.nx_pydot.write_dot(db.graph, dot_file)

        compact_graph = deepcopy(db.graph)
        for u, v, data in compact_graph.edges(data=True):
            data['coverages'] = None
            data['edge_kmer'] = None

        dot_compact_file = os.path.join(outdir, f'db_k{k}_compact.dot')
        nx.drawing.nx_pydot.write_dot(compact_graph, dot_compact_file)

        contigs, contig_paths = db.get_contigs()
        all_contigs[k] = contigs

        input_strings = strings.copy()
        for i in range(len(contigs)):
            for j in range(min_mult):
                input_strings[f'contig_k{k}_i{i}_j{j}'] = contigs[i]

        complex_kp1mers = get_paths_thru_complex_nodes(db, strings)

    return all_contigs, dbs, uncompr_dbs, all_frequent_kmers, all_frequent_kmers_read_pos


def scaffolding(db, mappings, min_connections=1, additional_edges=list()):
    def find_connections(additional_edges):
        long_edges = db.get_long_edges()
        long_edges_ids = list(long_edges.keys())
        long_edges_ids += additional_edges
        long_edges_ids = set(long_edges_ids)

        connections = defaultdict(lambda: defaultdict(int))
        for r_id, mapping in mappings.items():
            if mapping is None:
                continue
            _, _, valid_path, path = mapping
            if valid_path:
                mapped_edges = set(path)
                inters = mapped_edges & long_edges_ids
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
        for long_edge in db.get_long_edges():
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
            print(cc)
            cc_sg = scaffold_graph.subgraph(cc)
            if nx.is_directed_acyclic_graph(cc_sg):
                top_sort = list(nx.topological_sort(cc_sg))
                longest_path = nx.dag_longest_path(cc_sg)
                longedge_scaffolds.append(longest_path)
        return longedge_scaffolds

    def get_longest_extensions(longedge_scaffold):
        left_edge = longedge_scaffold[0]
        right_edge = longedge_scaffold[-1]
        lst_left_ext = []
        lst_right_ext = []
        for r_id, mapping in mappings.items():
            if mapping is None:
                continue
            _, _, valid_path, path = mapping
            if not valid_path:
                continue
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

    connections = find_connections(additional_edges=additional_edges)
    # for (e1, e2) in connections:
    #     print(e1, e2)
    #     print(connections[(e1, e2)])
    scaffold_graph = build_scaffold_graph(connections)
    nx.drawing.nx_pydot.write_dot(scaffold_graph, 'scaffold_graph.dot')
    longedge_scaffolds = select_lists(scaffold_graph)
    edge_scaffolds = get_edge_scaffolds(longedge_scaffolds,
                                        connections)
    scaffolds = get_sequence_scaffolds(edge_scaffolds)
    return scaffolds, edge_scaffolds


def read2scaffolds(db, scaffold_paths, mappings, monoreads):
    edgescaffolds2coords = [db.get_edgepath2coords(scaffold_path)
                            for scaffold_path in scaffold_paths]
    r2s = defaultdict(list)

    for r_id, mapping in mappings.items():
        if mapping is None:
            continue
        (e_st, r_st), (e_en, r_en), valid_path, read_path = mapping
        if not valid_path:
            continue
        for sc_index, scaffold_path in enumerate(scaffold_paths):
            edgescaffold2coords = edgescaffolds2coords[sc_index]
            for i in range(len(scaffold_path) - len(read_path) + 1):
                if scaffold_path[i:i+len(read_path)] == read_path:
                    # if r_id in r2s:
                    #     print(r_id, sc_index, i)
                    #     ambiguous_mappings.add(r_id)
                    r2s[r_id].append((
                        sc_index,
                        edgescaffold2coords[i, e_st[1]],
                        edgescaffold2coords[i+len(read_path)-1, e_en[1]+db.k-1]
                    ))
    r2s = {k: v[0] for k, v in r2s.items() if len(v) == 1}
    return r2s


def cover_scaffolds_w_reads(r2s, mappings, scaffold_seqs, monoreads, k):
    coverage = [[{} for i in range(len(scaffold_seq))] for scaffold_seq in scaffold_seqs]
    for r_id, (scaf_id, s_st, s_en) in r2s.items():
        (_, r_st), (_, r_en), valid_path, path = mappings[r_id]
        if not valid_path:
            continue
        # TODO: fixit
        if s_en - s_st != r_en - r_st + k - 1:
            continue
        cov_scaf = coverage[scaf_id]
        r_mono2nucl = monoreads[r_id].mono2nucl
        for i in range(s_en - s_st + 1):
            if r_st + i in r_mono2nucl:
                cov_scaf[s_st+i][r_id] = r_mono2nucl[r_st+i]
            else:
                # a corrected gap
                pass
    return coverage


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


def extract_read_pseudounits(scaf_read_coverage, scaffold_seqs,
                             monostrings, min_coverage=0):
    read_pseudounits, pseudounits = [], []
    for i, scaffold_seq in enumerate(scaffold_seqs):
        read_pseudounits.append(list())
        scaf_read_pseudounits = read_pseudounits[-1]
        sr_coverage = scaf_read_coverage[i]
        scaf_pseudounits = partition_pseudounits(scaffold_seq)
        pseudounits.append(scaf_pseudounits)
        for u_st, u_en in scaf_pseudounits:
            s_cov = sr_coverage[u_st]
            e_cov = sr_coverage[u_en]
            s_rids = set(s_cov.keys())
            e_rids = set(e_cov.keys())
            r_ids = s_rids & e_rids
            if len(r_ids) < min_coverage:
                continue
            scaf_read_pseudounits.append(dict())
            for r_id in r_ids:
                coords = s_cov[r_id][1:] + e_cov[r_id][1:]
                st, en = min(coords), max(coords)
                strand = monostrings[r_id].strand
                scaf_read_pseudounits[-1][r_id] = (st, en, strand)
    return pseudounits, read_pseudounits


def polish(scaffolds, pseudounits, read_pseudounits, reads, monomers, outdir,
           n_iter, n_threads, flye_bin='flye'):
    def get_template(scaffold, st, en):
        return ''.join(monomers[m_id] for m_id in scaffold[st:en+1])
    monomers = {m_id[0]: monomer
                for m_id, monomer in monomers.items()
                if m_id[-1] != "'"}
    smart_makedirs(outdir)
    for i, (scaffold, scaf_pseudounits) in enumerate(zip(scaffolds,
                                                         pseudounits)):
        scaf_outdir = os.path.join(outdir, f'scaffold_{i}')
        smart_makedirs(scaf_outdir)

        polished_scaffold = []
        for j, (s_st, s_en) in enumerate(scaf_pseudounits):
            pseudounit_outdir = os.path.join(scaf_outdir, f'pseudounit_{j}')
            smart_makedirs(pseudounit_outdir)

            # template = get_template(scaffold, s_st, s_en)
            # template_id = f'scaffold_{i}_template_{j}_{scaffold[s_st:s_en+1]}'
            # write_bio_seqs(template_fn, {template_id: template})

            pseudounit_reads = {}
            for r_id, (r_st, r_en, strand) in read_pseudounits[i][j].items():
                read_segm_id = f's_{i}_t_{j}_{r_id[0]}_{r_st}_{r_en+1}'
                pseudounit_read = reads[r_id[0]][r_st:r_en+1]
                if strand == '-':
                    pseudounit_read = RC(pseudounit_read)
                pseudounit_reads[read_segm_id] = pseudounit_read
            reads_fn = os.path.join(pseudounit_outdir, 'reads.fasta')
            write_bio_seqs(reads_fn, pseudounit_reads)

            template_fn = os.path.join(pseudounit_outdir, 'template.fasta')
            template_id, template_read = "", None
            r_units_lens = [len(read) for read in pseudounit_reads.values()]
            med_len = statistics.median_high(r_units_lens)
            for r_id in sorted(pseudounit_reads.keys()):
                read = pseudounit_reads[r_id]
                if len(read) == med_len:
                    template_id = r_id
                    template_read = read
                    break
            assert len(pseudounit_reads[template_id]) == med_len
            assert len(template_read) == med_len
            write_bio_seqs(template_fn,
                           {template_id: template_read})

            cmd = [flye_bin,
                   '--nano-raw', reads_fn,
                   '--polish-target', template_fn,
                   '-i', n_iter,
                   '-t', n_threads,
                   '-o', pseudounit_outdir]
            cmd = [str(x) for x in cmd]
            print(' '.join(cmd))
            subprocess.check_call(cmd)

            try:
                polished_pseudounit_fn = \
                    os.path.join(pseudounit_outdir,
                                 f'polished_{n_iter}.fasta')
                polished_pseudounit = read_bio_seq(polished_pseudounit_fn)
                polished_scaffold.append(polished_pseudounit)
            except FileNotFoundError:
                polished_scaffold.append(template)

        polished_scaffold = ''.join(polished_scaffold)
        polished_scaffold_fn = os.path.join(scaf_outdir, f'scaffold_{i}.fasta')
        write_bio_seqs(polished_scaffold_fn,
                       {f'scaffold_{i}_niter_{n_iter}': polished_scaffold})


# 3-colored graph for assembly debug
def graph3col_uncompr(gr_assembly, gr_reads):
    gr = nx.MultiDiGraph()
    k = gr_assembly.k

    nodes_assembly = set(gr_assembly.rev_node_mapping[i] for i in gr_assembly.graph.nodes)
    nodes_reads = set(gr_reads.rev_node_mapping[i] for i in gr_reads.graph.nodes)
    nodes_both = nodes_assembly & nodes_reads
    nodes_assembly = nodes_assembly - nodes_both
    nodes_reads = nodes_reads - nodes_both
    nodes_all = nodes_assembly | nodes_reads | nodes_both
    print(len(nodes_both), len(nodes_assembly), len(nodes_reads))
    nodes_index = {node: i for i, node in enumerate(nodes_all)}

    edges_assembly = set(edge[2]['edge_kmer'] for edge in gr_assembly.graph.edges(data=True))
    edges_assembly_cov = {edge[2]['edge_kmer']: edge[2]['coverages']
                          for edge in gr_assembly.graph.edges(data=True)}
    edges_reads = set(edge[2]['edge_kmer'] for edge in gr_reads.graph.edges(data=True))
    edges_reads_cov = {edge[2]['edge_kmer']: edge[2]['coverages']
                       for edge in gr_reads.graph.edges(data=True)}
    edges_both = edges_assembly & edges_reads
    edges_assembly = edges_assembly - edges_both
    edges_reads = edges_reads - edges_both
    edges_all = edges_assembly | edges_reads | edges_both
    print(len(edges_both), len(edges_assembly), len(edges_reads))
    edges_index = {edge: i for i, edge in enumerate(edges_all)}

    def add_nodes(nodes, color):
        for node in nodes:
            gr.add_node(nodes_index[node], kmer=node, color=color)

    add_nodes(nodes_reads, color='red')
    add_nodes(nodes_assembly, color='blue')
    add_nodes(nodes_both, color='green')

    def add_edges(edges, color):
        for edge in edges:
            stkmer = edge[:k-1]
            enkmer = edge[-k+1:]
            strlen = len(edges_assembly_cov[edge]) \
                     if edge in edges_assembly_cov \
                     else len(edges_reads_cov[edge])
            label = ['len' + str(strlen)]

            assembly_cov, reads_cov = [], []
            if edge in edges_assembly_cov:
                assembly_cov = edges_assembly_cov[edge]
                label.append('assembly_cov' + str(np.median(assembly_cov)))
            if edge in edges_reads_cov:
                reads_cov = edges_reads_cov[edge]
                label.append('reads_cov' + str(np.median(reads_cov)))

            label = '\n'.join(label)
            gr.add_edge(nodes_index[edge[:k-1]],
                        nodes_index[edge[-k+1:]],
                        string=edge,
                        assembly_cov=assembly_cov,
                        reads_cov=reads_cov,
                        color=color,
                        label=label)
    add_edges(edges_reads, color='red')
    add_edges(edges_assembly, color='blue')
    add_edges(edges_both, color='green')
    return gr


def collapse_nonbranching_paths_3color(graph, k):
    def node_on_nonbranching_path(graph, node):
        colors = nx.get_edge_attributes(graph, 'color')
        if graph.in_degree(node) != 1 or graph.out_degree(node) != 1:
            return False
        in_edge = list(graph.in_edges(node, keys=True))[0]
        out_edge = list(graph.out_edges(node, keys=True))[0]
        return in_edge != out_edge \
            and graph.edges[in_edge]['color'] == graph.edges[out_edge]['color']

    for node in list(graph.nodes()):
        if node_on_nonbranching_path(graph, node):
            in_edge = list(graph.in_edges(node, keys=True))[0]
            out_edge = list(graph.out_edges(node, keys=True))[0]
            in_edge_color = graph.edges[in_edge]['color']
            out_edge_color = graph.edges[out_edge]['color']
            assert in_edge_color == out_edge_color
            in_edge_kmer = graph.edges[in_edge]['string']
            out_edge_kmer = graph.edges[out_edge]['string']

            in_edge_cov_reads = graph.edges[in_edge]['reads_cov']
            out_edge_cov_reads = graph.edges[out_edge]['reads_cov']
            in_edge_cov_assembly = graph.edges[in_edge]['assembly_cov']
            out_edge_cov_assembly = graph.edges[out_edge]['assembly_cov']

            in_node = in_edge[0]
            out_node = out_edge[1]

            new_kmer = in_edge_kmer + \
                out_edge_kmer[-(len(out_edge_kmer)-k+1):]
            new_cov_reads = sorted(in_edge_cov_reads + out_edge_cov_reads)
            new_cov_assembly = sorted(in_edge_cov_assembly + out_edge_cov_assembly)
            if len(new_cov_reads) and len(new_cov_assembly):
                assert len(new_cov_reads) == len(new_cov_assembly)
            new_edge_len = max(len(new_cov_reads), len(new_cov_assembly))

            new_edge_color = in_edge_color # == out_edge_color

            new_edge_label = ['len' + str(new_edge_len)]
            if len(new_cov_assembly):
                new_edge_label.append('assembly_cov' + str(np.median(new_cov_assembly)))
            if len(new_cov_reads):
                new_edge_label.append('reads_cov' + str(np.median(new_cov_reads)))
            new_edge_label = '\n'.join(new_edge_label)

            graph.add_edge(
                in_node, out_node,
                string=new_kmer,
                assembly_cov=new_cov_assembly,
                reads_cov=new_cov_reads,
                length=new_edge_len,
                label=new_edge_label,
                color=new_edge_color)
            graph.remove_node(node)
    return graph


def graph3col(gr_assembly, gr_reads):
    gr = graph3col_uncompr(gr_assembly, gr_reads)
    k = gr_assembly.k
    compr_gr = collapse_nonbranching_paths_3color(gr, k=k)
    return compr_gr