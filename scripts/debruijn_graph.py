from collections import defaultdict, Counter
from itertools import groupby
import os

import networkx as nx
import numpy as np

from utils.os_utils import smart_makedirs


class DeBruijnGraph:
    def __init__(self, k, max_uniq_cov=30, min_uniq_len=1000):
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
            return nx.number_of_nodes(graph) > 1 \
                and graph.in_degree(node) == 1 \
                and graph.out_degree(node) == 1

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
                new_coverages = in_edge_cov + out_edge_cov
                new_coverages.sort()
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
        path = ''
        path += self.graph.get_edge_data(*list_edges[0])['edge_kmer']
        for edge in list_edges[1:]:
            new_edge_seq = self.graph.get_edge_data(*edge)['edge_kmer']
            assert path[-(self.k-1):] == new_edge_seq[:self.k-1]
            path += new_edge_seq[self.k-1:]

        # process cycled path
        if list_edges[0][0] == list_edges[-1][1]:
            path = path[:-(self.k-1)]
        return path

    def get_contigs(self):
        def get_longest_valid_outpaths(graph):
            def get_valid_outpath_edge(edge, taken_edges=set()):
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

    def map_reads(self, monomer_strings, verbose=True):
        print("Indexing started")
        self.index_edges()
        print("Indexing completed")
        mapping = {}
        db_edges = list(self.graph.edges(keys=True))
        for i, (r_id, string) in enumerate(monomer_strings.items()):
            if verbose:
                print(i+1, len(monomer_strings))

            split_strings = list(filter(lambda string: len(string), string.split('=')))
            split_lens = [0] + [len(split_string) for split_string in split_strings]
            cum_split_lens = np.cumsum(split_lens)
            read_coords = []
            for split_ind, split_string in enumerate(split_strings):
                for i in range(len(split_string)-self.k+1):
                    kmer = split_string[i:i+self.k]
                    if kmer in self.db_index[len(kmer)]:
                        read_coords.append(self.db_index[len(kmer)][kmer])

            path = [x[0] for x in read_coords]
            path = [x[0] for x in groupby(path)]
            path = [db_edges[edge_ind] for edge_ind in path]

            valid_path = True
            for e1, e2 in zip(path[:-1], path[1:]):
                if e1[1] != e2[0]:
                    valid_path = False
                    break
            if len(read_coords):
                mapping[r_id] = (read_coords[0], read_coords[-1], valid_path, path)
            else:
                mapping[r_id] = None
        return mapping


def get_all_kmers(strings, k, gap_symb='?'):
    all_kmers = Counter()
    read_kmer_locations = defaultdict(list)
    for r_id, string in strings.items():
        for i in range(len(string)-k+1):
            kmer = string[i:i+k]
            if gap_symb not in kmer:
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
                kp1 = in_kmer + out_kmer[-1]
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


def iterative_graph(strings, min_k, max_k, outdir,
                    min_mult=5, step=1, starting_graph=None, verbose=True):
    smart_makedirs(outdir)
    dbs, all_contigs = {}, {}
    all_frequent_kmers, all_frequent_kmers_read_pos = {}, {}
    input_strings = strings.copy()
    complex_kp1mers = {}

    if starting_graph is not None:
        contigs, contig_paths = starting_graph.get_contigs()
        for i in range(len(contigs)):
            for j in range(min_mult):
                input_strings[f'contig_k{min_k}_i{i}_j{j}'] = contigs[i]

        complex_kp1mers = get_paths_thru_complex_nodes(starting_graph, strings)

    for k in range(min_k, max_k+1, step):
        frequent_kmers, frequent_kmers_read_pos = \
            get_frequent_kmers(input_strings, k=k, min_mult=min_mult)
        frequent_kmers.update(complex_kp1mers)
        if verbose:
            print(f'\nk={k}')
            print(f'#frequent kmers = {len(frequent_kmers)}')
        all_frequent_kmers[k] = frequent_kmers
        all_frequent_kmers_read_pos[k] = frequent_kmers_read_pos

        db = DeBruijnGraph(k=k)
        db.add_kmers(frequent_kmers, coverage=frequent_kmers)

        db.collapse_nonbranching_paths()
        if nx.number_weakly_connected_components(db.graph) > 1:
            print(f'#cc = {nx.number_weakly_connected_components(db.graph)}')
            for cc in nx.weakly_connected_components(db.graph):
                print(len(cc))
            # break
        dbs[k] = db

        dot_file = os.path.join(outdir, f'db_k{k}.dot')
        # pdf_file = os.path.join(outdir, f'db_k{k}.pdf')
        nx.drawing.nx_pydot.write_dot(db.graph, dot_file)
        # os.system(f"dot -Tpdf {dot_file} -o {pdf_file}")

        contigs, contig_paths = db.get_contigs()
        all_contigs[k] = contigs

        input_strings = strings.copy()
        for i in range(len(contigs)):
            for j in range(min_mult):
                input_strings[f'contig_k{k}_i{i}_j{j}'] = contigs[i]

        complex_kp1mers = get_paths_thru_complex_nodes(db, strings)

    return all_contigs, dbs, all_frequent_kmers, all_frequent_kmers_read_pos
