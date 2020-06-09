# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from collections import defaultdict
import networkx as nx
import numpy as np

from sequence_graph.sequence_graph import SequenceGraph
from utils.nx_all_simple_paths_multigraph import all_simple_edge_paths


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
