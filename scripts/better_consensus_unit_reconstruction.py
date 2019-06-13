#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

import argparse

import numpy as np
from collections import defaultdict
import heapq

from utils.bio import read_bio_seq, write_bio_seqs
from ncrf_parser import NCRF_Report

import edlib
import networkx as nx


class DeBruijnGraph:
    def __init__(self, k):
        self.graph = nx.MultiDiGraph()
        self.k = k

    def add_kmer(self, kmer, color='black', coverage=1):
        self.graph.add_edge(kmer[:-1], kmer[1:], edge_kmer=kmer, color=color, coverages=[coverage])

    def add_kmers(self, kmers, color='black', coverage=None):
        for kmer in kmers:
            if coverage is None:
                self.add_kmer(kmer, color=color)
            else:
                self.add_kmer(kmer, color=color, coverage=coverage[kmer])

    def remove_tips(self):
        while True:
            nodes2delete = []
            for node in self.graph.nodes:
                if self.graph.out_degree(node) == 0 and self.graph.in_degree(node) == 0:
                    continue
                if self.graph.out_degree(node) == 0:
                    in_edge = list(self.graph.in_edges(node, keys=True))[0]
                    #in_edge_color = self.graph.edges[in_edge]['color']
                    #if in_edge_color == 'blue':
                    #    continue
                    nodes2delete.append(node)
                elif self.graph.in_degree(node) == 0:
                    out_edge = list(self.graph.out_edges(node, keys=True))[0]
                    #out_edge_color = self.graph.edges[out_edge]['color']
                    #if out_edge_color == 'blue':
                    #    continue
                    nodes2delete.append(node)
            if len(nodes2delete):
                self.graph.remove_nodes_from(nodes2delete)
            else:
                return

    def collapse_nonbranching_paths(self, respect_color=True):
        def node_on_nonbranching_path(graph, node):
            return nx.number_of_nodes(graph) > 1 and graph.in_degree(node) == 1 and graph.out_degree(node) == 1

        for node in list(self.graph.nodes()):
            if node_on_nonbranching_path(self.graph, node):
                in_edge = list(self.graph.in_edges(node, keys=True))[0]
                out_edge = list(self.graph.out_edges(node, keys=True))[0]
                in_edge_color = self.graph.edges[in_edge]['color']
                out_edge_color = self.graph.edges[out_edge]['color']
                in_edge_kmer = self.graph.edges[in_edge]['edge_kmer']
                out_edge_kmer = self.graph.edges[out_edge]['edge_kmer']
                in_edge_cov = self.graph.edges[in_edge]['coverages']
                out_edge_cov = self.graph.edges[out_edge]['coverages']

                in_node = in_edge[0]
                out_node = out_edge[1]
                #print(in_edge_kmer, out_edge_kmer)
                if (not respect_color) or (in_edge_color == out_edge_color):
                    new_kmer = in_edge_kmer + out_edge_kmer[-(len(out_edge_kmer)-self.k+1):]
                    new_coverages = in_edge_cov + out_edge_cov
                    new_coverages.sort()
                    self.graph.add_edge(in_node, out_node,
                                        edge_kmer=new_kmer,
                                        color=in_edge_color,
                                        coverages=new_coverages)
                    self.graph.remove_node(node)

    def purify_graph(self):
        first_edge = None
        for edge in self.graph.edges:
            if self.graph.out_degree(edge[0]) == 1 and self.graph.in_degree(edge[1]) == 1:
                first_edge = edge
                break
        assert first_edge is not None
        first_edge_properties = self.graph.edges[first_edge]
        self.graph.remove_edge(*first_edge)

        while True:
            coverages = get_coverage(self.graph)
            min_coverage = np.inf
            edges = list(self.graph.edges)
            edge2remove = None
            for edge in edges:
                gr = self.graph.copy()
                gr.remove_edge(*edge)
                if nx.is_weakly_connected(gr):
                    if coverages[edge] < min_coverage:
                        edge2remove = edge
                        min_coverage = coverages[edge]
            if edge2remove is not None:
                print(edge2remove)
                self.graph.remove_edge(*edge2remove)
                assert nx.is_weakly_connected(self.graph)
                self.graph.remove_nodes_from(list(nx.isolates(self.graph)))
                #self.collapse_nonbranching_paths()
            else:
                break

        self.graph.add_edge(*first_edge, **first_edge_properties)
        self.remove_tips()
        self.collapse_nonbranching_paths(respect_color=False)


def get_coverage(graph):
    coverages = {}
    for edge in graph.edges:
        coverages[edge] = np.min(graph.edges[edge]['coverages'])
    return coverages


def get_kmer_counts_reads(ncrf_report, k=19):
    kmer_counts_reads = defaultdict(int)
    for i, (r_id, record) in enumerate(ncrf_report.records.items()):
        if (i % 100 == 0):
            print(i + 1, len(ncrf_report.records))
        r_al = record.r_al
        r_al = r_al.replace('-', '')
        read_kmers = defaultdict(int)
        for i in range(len(r_al)-k+1):
            kmer = r_al[i:i+k]
            kmer_counts_reads[kmer] += 1
    return kmer_counts_reads

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads-ncrf", help="NCRF report on centromeric reads", required=True)
    parser.add_argument("--unit", help="Initial unit sequence of centromeric read", required=True)
    parser.add_argument("-k", type=int, default=30)
    parser.add_argument("--output", help="Output file for polished unit", required=True)
    params = parser.parse_args()

    reads_ncrf_report = NCRF_Report(params.reads_ncrf)
    unit_seq = read_bio_seq(params.unit)

    kmer_counts_reads = get_kmer_counts_reads(reads_ncrf_report, k=params.k)
    unit_double_seq = unit_seq + unit_seq
    unit_kmers = set(unit_double_seq[i:i+params.k] for i in range(len(unit_seq)))
    len(unit_kmers)

    most_frequent_kmers = heapq.nlargest(n=int(len(unit_kmers)*3),
                                         iterable=kmer_counts_reads,
                                         key=kmer_counts_reads.get)
    most_frequent_kmers_cnts = sorted([kmer_counts_reads[kmer] for kmer in most_frequent_kmers], reverse=True)

    debr = DeBruijnGraph(k=params.k)
    most_frequent_kmers = set(most_frequent_kmers)
    unit_kmers = set(unit_kmers)
    common_kmers = most_frequent_kmers & unit_kmers

    #debr.add_kmers(common_kmers, 'green', kmer_counts_reads)
    #debr.add_kmers(unit_kmers - common_kmers, 'blue', kmer_counts_reads)
    #debr.add_kmers(most_frequent_kmers - common_kmers, 'red', kmer_counts_reads)
    debr.add_kmers(most_frequent_kmers, 'red', kmer_counts_reads)

    debr.collapse_nonbranching_paths()
    debr.remove_tips()
    debr.collapse_nonbranching_paths()
    debr.purify_graph()

    edge = list(debr.graph.edges)[0]
    new_unit = debr.graph.edges[edge]['edge_kmer']

    trim = 1
    while new_unit[:trim] != new_unit[-trim:]:
        trim += 1
    new_unit = new_unit[:-trim]


    doubled_unit = new_unit + new_unit
    alignment_locations = edlib.align(unit_seq, doubled_unit, mode='HW', task='locations')['locations'][0]
    new_unit = doubled_unit[alignment_locations[0]:alignment_locations[0]+len(new_unit)]

    write_bio_seqs(params.output, {'updated_unit': new_unit})


if __name__ == "__main__":
    main()
