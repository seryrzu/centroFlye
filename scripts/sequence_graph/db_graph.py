# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

from collections import defaultdict
import networkx as nx
import numpy as np

from sequence_graph.sequence_graph import SequenceGraph

logger = logging.getLogger("centroFlye.sequence_graph.db_graph")


class DeBruijnGraph(SequenceGraph):
    coverage = 'coverage'

    def __init__(self, nx_graph, nodeindex2label, nodelabel2index, k,
                 collapse=True, resolve_bubbles=True):
        super().__init__(nx_graph=nx_graph,
                         nodeindex2label=nodeindex2label,
                         nodelabel2index=nodelabel2index,
                         collapse=collapse)
        self.k = k  # length of an edge in the uncompressed graph
        if resolve_bubbles:
            self.resolve_bubbles()

    @classmethod
    def _generate_label(cls, par_dict):
        cov = par_dict[cls.coverage]
        length = par_dict[cls.length]

        mean_cov = np.mean(cov)
        label = f'len={length}\ncov={mean_cov:0.2f}'
        return label

    @classmethod
    def from_kmers(cls, kmers, kmer_coverages=None,
                   min_tip_cov=1, collapse=True, resolve_bubbles=False):
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
                       collapse=collapse,
                       resolve_bubbles=resolve_bubbles)
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

    def resolve_bubbles(self, max_diff=1):
        edge2strCov = defaultdict(list)
        for s, e, key, data in self.nx_graph.edges(keys=True, data=True):
            assert len(edge2strCov[(s, e)]) == key
            string = data[self.string]
            coverage = np.mean(data[self.coverage])
            edge2strCov[(s, e)].append((string, coverage))
        edge2strCov = {edge: tuple(v) for edge, v in edge2strCov.items()
                       if len(v) == 2 and len(v[0][0]) == len(v[1][0])}

        resolved_bubbles = []
        for (s, e), ((str1, cov1), (str2, cov2)) in edge2strCov.items():
            diff = [i for i, (a, b) in enumerate(zip(str1, str2)) if a != b]
            if len(diff) > max_diff:
                continue
            if cov1 >= cov2:
                self.nx_graph.remove_edge(s, e, 1)
                resolved_bubbles.append((str2, str1))
            else:
                self.nx_graph.remove_edge(s, e, 0)
                resolved_bubbles.append((str1, str2))
        self.collapse_nonbranching_paths()
        return resolved_bubbles
