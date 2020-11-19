# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
from collections import defaultdict
import logging
import math
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

import ahocorasick
from ahocorapy.keywordtree import KeywordTree
import networkx as nx
import edlib

from sequence_graph.path_graph import IDBMappings
from standard_logger import get_logger
from subprocess import call
from utils.bio import read_bio_seqs, read_bio_seq, compress_homopolymer, RC
from utils.git import get_git_revision_short_hash
from utils.os_utils import smart_makedirs, expandpath
from utils.various import fst_iterable


logger = logging.getLogger("centroFlye.sequence_graph.path_graph_multik")


class PathMultiKGraph:
    SATURATING_K = 40002

    def __init__(self, nx_graph,
                 edge2seq, edge2index, index2edge, node2len,
                 max_edge_index, max_node_index,
                 idb_mappings,
                 init_k,
                 unique_edges=None):
        self.nx_graph = nx_graph
        self.edge2seq = edge2seq
        self.edge2index = edge2index
        self.index2edge = index2edge
        self.node2len = node2len
        self.max_edge_index = max_edge_index
        self.max_node_index = max_node_index
        self.idb_mappings = idb_mappings
        self.init_k = init_k
        if unique_edges is None:
            unique_edges = set()
        self.unique_edges = unique_edges

        self.niter = 0

        self.unresolved = set()
        self._update_unresolved_vertices()

        self.assert_validity()

    def assert_validity(self):
        self.max_edge_index == 1 + max(self.index2edge)
        self.max_node_index == 1 + max(self.nx_graph.nodes)
        edges = set(self.nx_graph.edges(keys=True))
        assert edges == set(self.edge2index.keys())
        assert edges == set(self.index2edge.values())
        assert all(edge == self.index2edge[self.edge2index[edge]]
                   for edge in self.edge2index)

        ac = self.idb_mappings.get_active_connections()
        for i, j in ac:
            _, u, _ = self.index2edge[i]
            v, _, _ = self.index2edge[j]
            assert u == v

        for node in self.nx_graph.nodes:
            assert node in self.node2len
            nlen = self.node2len[node]
            assert self.nx_graph.in_degree(node) != 1 or \
                self.nx_graph.out_degree(node) != 1
            for in_edge in self.nx_graph.in_edges(node, keys=True):
                e_index = self.edge2index[in_edge]
                in_seq = self.edge2seq[e_index]
                insuf = in_seq[-nlen:]
                for out_edge in self.nx_graph.out_edges(node, keys=True):
                    e_outdex = self.edge2index[out_edge]
                    out_seq = self.edge2seq[e_outdex]
                    outpref = out_seq[:nlen]
                    assert insuf == outpref

    @classmethod
    def fromDB(cls, db, string_set, neutral_symbs=None, raw_mappings=None):
        if raw_mappings is None:
            raw_mappings = db.map_strings(string_set,
                                          only_unique_paths=True,
                                          neutral_symbs=neutral_symbs)

        nx_graph = nx.MultiDiGraph()
        edge2seq = {}
        edge_index = 0
        edge2index = {}
        index2edge = {}
        node2len = {}

        for u in db.nx_graph.nodes():
            node2len[u] = db.k-1

        for u, v, key, data in db.nx_graph.edges(data=True, keys=True):
            nx_graph.add_edge(u, v, key)
            edge2index[(u, v, key)] = edge_index
            index2edge[edge_index] = (u, v, key)
            edge2seq[edge_index] = list(data['string'])
            edge_index += 1

        mappings = {}
        for r_id, (raw_mapping, _, _) in raw_mappings.items():
            mappings[r_id] = [edge2index[edge] for edge in raw_mapping]

        idb_mappings = IDBMappings(mappings)
        max_node_index = 1 + max(nx_graph.nodes)

        return cls(nx_graph=nx_graph,
                   edge2seq=edge2seq,
                   edge2index=edge2index,
                   index2edge=index2edge,
                   node2len=node2len,
                   max_edge_index=edge_index,
                   max_node_index=max_node_index,
                   idb_mappings=idb_mappings,
                   init_k=db.k)

    @classmethod
    def from_mono_db(cls, db, monostring_set,
                     mappings=None):
        monostring = fst_iterable(monostring_set.values())
        neutral_symbs = set([monostring.gap_symb])
        return cls.fromDB(db=db,
                          string_set=monostring_set,
                          neutral_symbs=neutral_symbs,
                          raw_mappings=mappings)

    @classmethod
    def fromDR(cls, db_fn, align_fn, k):
        class PerfectHash:
            hash2index = {}
            next_key = 0

            def __getitem__(self, key):
                if key in self.hash2index:
                    return self.hash2index[key]
                self.hash2index[key] = self.next_key
                self.next_key += 1
                return self.hash2index[key]

        db = read_bio_seqs(db_fn)
        nx_graph = nx.MultiDiGraph()
        edge2seq = {}
        edge2index = {}
        index2edge = {}
        node2len = {}
        ph = PerfectHash()
        max_edge_index = 0
        vertex_nucl2edgeindex = {}
        edge2cov = {}
        for e_id, seq in db.items():
            split_id = e_id.split('_')
            index, s, e, _, cov = split_id
            index, s, e = int(index), int(s), int(e)
            cov = float(cov)
            s = ph[s]
            e = ph[e]
            key = nx_graph.add_edge(s, e)
            edge = (s, e, key)
            edge2seq[index] = list(seq)
            edge2index[edge] = index
            index2edge[index] = edge
            edge2cov[index] = cov
            max_edge_index = max(index, max_edge_index)
            vertex_nucl2edgeindex[(s, seq[k-1])] = index
            node2len[s] = k - 1
            node2len[e] = k - 1

        max_edge_index += 1

        average_cov = sum(edge2cov[index] * len(edge2seq[index])
                          for index in edge2cov) / \
            sum(len(edge2seq[index]) for index in edge2cov)
        logger.info(f'Average coverage {average_cov}')
        unique_edges = set([index for index, cov in edge2cov.items()
                            if cov <= average_cov * 1.2])

        mappings = {}
        with open(align_fn) as f:
            for i, line in enumerate(f):
                s_id, u, nucl = line.strip().split()
                u = ph[int(u)]
                mapping = []
                for c in nucl:
                    e = vertex_nucl2edgeindex[(u, c)]
                    mapping.append(e)
                    u = index2edge[e][1]
                mappings[s_id] = mapping

        idb_mappings = IDBMappings(mappings)

        return cls(nx_graph=nx_graph,
                   edge2seq=edge2seq,
                   edge2index=edge2index,
                   index2edge=index2edge,
                   node2len=node2len,
                   max_edge_index=max_edge_index,
                   max_node_index=ph.next_key,
                   init_k=k,
                   idb_mappings=idb_mappings,
                   unique_edges=unique_edges)

    def move_edge(self, e1_st, e1_en, e1_key,
                  e2_st, e2_en, e2_key=None):
        old_edge = (e1_st, e1_en, e1_key)
        i = self.edge2index[old_edge]
        self.nx_graph.remove_edge(*old_edge)
        e2_key = self.nx_graph.add_edge(e2_st, e2_en, key=e2_key)
        new_edge = (e2_st, e2_en, e2_key)
        self.edge2index[new_edge] = i
        del self.edge2index[old_edge]
        self.index2edge[i] = new_edge

    def remove_edge(self, edge=None, index=None, moving=True):
        assert (edge is None) != (index is None)
        if edge is None:
            edge = self.index2edge[index]
        else:
            index = self.edge2index[edge]
        self.idb_mappings.remove(index)
        self.nx_graph.remove_edge(*edge)
        del self.edge2index[edge]
        del self.index2edge[index]
        del self.edge2seq[index]

        self.unique_edges.discard(index)

        if moving:
            for e in list(self.nx_graph.in_edges(edge[1], keys=True)):
                self.move_edge(*e, e[0], edge[0])

            for e in list(self.nx_graph.out_edges(edge[1], keys=True)):
                self.move_edge(*e, edge[0], e[1])

    def get_new_vertex_index(self):
        index = self.max_node_index
        self.node2len[index] = self.init_k - 1 + self.niter
        self.max_node_index += 1
        return index

    def merge_edges(self, e1, e2):
        assert self.nx_graph.degree(e1[1]) == 1
        assert self.nx_graph.degree(e2[0]) == 1
        i = self.edge2index[e1]
        j = self.edge2index[e2]
        self.idb_mappings.merge(i, j)
        if j in self.unique_edges:
            self.unique_edges.add(i)
        in_seq = self.edge2seq[i]
        out_seq = self.edge2seq[j]
        nlen = self.node2len[e2[0]]-1  # need -1 since new nodes have ++
        assert in_seq[-nlen:] == out_seq[:nlen]
        seq = in_seq + out_seq[nlen:]
        self.edge2seq[i] = seq
        self.move_edge(*e1, e1[0], e2[1])
        self.remove_edge(edge=e2, moving=False)

    def add_edge(self, i, j, seq):
        in_edge = self.index2edge[i]
        out_edge = self.index2edge[j]
        new_edge = (in_edge[1], out_edge[0])
        key = self.nx_graph.add_edge(*new_edge)
        new_edge = (*new_edge, key)
        self.edge2index[new_edge] = self.max_edge_index
        self.index2edge[self.max_edge_index] = new_edge
        self.edge2seq[self.max_edge_index] = seq
        self.idb_mappings.add(i, j, self.max_edge_index)
        self.max_edge_index += 1

    def _update_unresolved_vertices(self):
        for u in self.nx_graph.nodes:
            if u in self.unresolved:
                continue
            in_indexes = set(
                [self.edge2index[e_in]
                 for e_in in self.nx_graph.in_edges(u, keys=True)])
            out_indexes = set(
                [self.edge2index[e_out]
                 for e_out in self.nx_graph.out_edges(u, keys=True)])
            indegree = self.nx_graph.in_degree(u)
            outdegree = self.nx_graph.out_degree(u)
            if indegree >= 2 and outdegree >= 2:
                # do not process anything at all
                # self.unresolved.add(u)

                # process only fully resolved vertices
                # all_ac = self.idb_mappings.get_active_connections()
                # pairs = set()
                # for e_in in in_indexes:
                #     for e_out in out_indexes:
                #         if (e_in, e_out) in all_ac:
                #             pairs.add((e_in, e_out))
                # all_pairs = set((e_in, e_out)
                #                 for e_in in in_indexes
                #                 for e_out in out_indexes)
                # if len(all_pairs - pairs):
                #     self.unresolved.add(u)

                # initial heuristic
                paired_in = set()
                paired_out = set()
                loops = in_indexes & out_indexes
                if len(loops) == 1:
                    loop = fst_iterable(loops)
                    if loop in self.unique_edges:
                        rest_in = in_indexes - loops
                        rest_out = out_indexes - loops
                        if len(rest_in) == 1:
                            in_index = fst_iterable(rest_in)
                            paired_in.add(in_index)
                            paired_out.add(loop)
                        if len(rest_out) == 1:
                            out_index = fst_iterable(rest_out)
                            paired_in.add(loop)
                            paired_out.add(out_index)
                for e_in in in_indexes:
                    for e_out in out_indexes:
                        if (e_in, e_out) in self.idb_mappings.pairindex2pos:
                            paired_in.add(e_in)
                            paired_out.add(e_out)
                tips = (in_indexes - paired_in) | (out_indexes - paired_out)

                if len(tips):
                    self.unresolved.add(u)

        prev, new = self.unresolved, set()
        while len(prev):
            for u in prev:
                for edge in self.nx_graph.in_edges(u, keys=True):
                    index = self.edge2index[edge]
                    seq = self.edge2seq[index]
                    v = edge[0]
                    if v in self.unresolved:
                        continue
                    if self.node2len[v] + 1 == len(seq):
                        new.add(v)
                for edge in self.nx_graph.out_edges(u, keys=True):
                    index = self.edge2index[edge]
                    seq = self.edge2seq[index]
                    v = edge[1]
                    if v in self.unresolved:
                        continue
                    if self.node2len[v] + 1 == len(seq):
                        new.add(v)
            self.unresolved |= new
            prev, new = new, set()

    def __process_vertex(self, u):
        def process_simple():
            if indegree == 1 and outdegree == 1:
                # node on nonbranching path - should not be happening
                assert False

            if indegree == 0 and outdegree == 0:
                # isolate - should be removed
                self.nx_graph.remove_node(u)
                del self.node2len[u]
                return

            elif indegree == 0 and outdegree > 0:
                # starting vertex
                for j in out_indexes[1:]:
                    old_edge = self.index2edge[j]
                    new_edge = (self.get_new_vertex_index(), old_edge[1], 0)
                    self.move_edge(*old_edge, *new_edge)

            elif indegree > 0 and outdegree == 0:
                # ending vertex
                for i in in_indexes[1:]:
                    old_edge = self.index2edge[i]
                    new_edge = (old_edge[0], self.get_new_vertex_index(), 0)
                    self.move_edge(*old_edge, *new_edge)

            elif indegree == 1 and outdegree > 1:
                # simple 1-in vertex
                assert len(in_indexes) == 1
                in_index = in_indexes[0]
                in_seq = self.edge2seq[in_index]
                c = in_seq[-nlen-1]
                for j in out_indexes:
                    assert self.edge2seq[j][:nlen] == in_seq[-nlen:]
                    self.edge2seq[j].insert(0, c)

            elif indegree > 1 and outdegree == 1:
                # simple 1-out vertex
                assert len(out_indexes) == 1
                out_index = out_indexes[0]
                out_seq = self.edge2seq[out_index]
                c = out_seq[nlen]
                for i in in_indexes:
                    assert self.edge2seq[i][-nlen:] == out_seq[:nlen]
                    self.edge2seq[i].append(c)
            self.node2len[u] += 1

        def process_complex():
            # complex vertex
            for i in in_indexes:
                old_edge = self.index2edge[i]
                new_edge = (old_edge[0], self.get_new_vertex_index(), 0)
                self.move_edge(*old_edge, *new_edge)

            for j in out_indexes:
                old_edge = self.index2edge[j]
                new_edge = (self.get_new_vertex_index(), old_edge[1], 0)
                self.move_edge(*old_edge, *new_edge)

            ac_s2e = defaultdict(set)
            ac_e2s = defaultdict(set)
            for e_in in in_indexes:
                for e_out in out_indexes:
                    if (e_in, e_out) in self.idb_mappings.pairindex2pos:
                        ac_s2e[e_in].add(e_out)
                        ac_e2s[e_out].add(e_in)

            loops = set(in_indexes) & set(out_indexes)
            if len(loops) == 1:
                loop = fst_iterable(loops)
                if loop in self.unique_edges:
                    rest_in = set(in_indexes) - loops
                    rest_out = set(out_indexes) - loops
                    if len(rest_in) == 1:
                        in_index = fst_iterable(rest_in)
                        ac_s2e[in_index].add(loop)
                        ac_e2s[loop].add(in_index)
                    if len(rest_out) == 1:
                        out_index = fst_iterable(rest_out)
                        ac_s2e[loop].add(out_index)
                        ac_e2s[out_index].add(loop)

            merged = {}
            for i in ac_s2e:
                for j in ac_s2e[i]:
                    if i in merged:
                        i = merged[i]
                    if j in merged:
                        j = merged[j]
                    e_i = self.index2edge[i]
                    e_j = self.index2edge[j]
                    in_seq = self.edge2seq[i]
                    out_seq = self.edge2seq[j]
                    assert in_seq[-nlen:] == out_seq[:nlen]
                    if len(ac_s2e[i]) == len(ac_e2s[j]) == 1:
                        self.merge_edges(e_i, e_j)
                        merged[j] = i
                    elif len(ac_s2e[i]) >= 2 and len(ac_e2s[j]) >= 2:
                        seq = in_seq[-nlen-1:] + [out_seq[nlen]]
                        assert len(seq) == nlen + 2
                        self.add_edge(i, j, seq)
                    elif len(ac_s2e[i]) == 1 and len(ac_e2s[j]) >= 2:
                        # extend left edge to the right
                        self.move_edge(*e_i, e_i[0], e_j[0])
                        seq = in_seq + [out_seq[nlen]]
                        self.edge2seq[i] = seq
                    elif len(ac_e2s[j]) == 1 and len(ac_s2e[i]) >= 2:
                        # extend right edge to the left
                        self.move_edge(*e_j, e_i[1], e_j[1])
                        seq = [in_seq[-nlen-1]] + out_seq
                        self.edge2seq[j] = seq
                    else:
                        assert False

            assert self.nx_graph.in_degree(u) == 0
            assert self.nx_graph.out_degree(u) == 0
            self.nx_graph.remove_node(u)
            del self.node2len[u]

        in_indexes = [self.edge2index[e_in]
                      for e_in in self.nx_graph.in_edges(u, keys=True)]
        out_indexes = [self.edge2index[e_out]
                       for e_out in self.nx_graph.out_edges(u, keys=True)]

        indegree = self.nx_graph.in_degree(u)
        outdegree = self.nx_graph.out_degree(u)

        nlen = self.node2len[u]

        if indegree >= 2 and outdegree >= 2:
            process_complex()
        else:
            process_simple()

    def transform_single(self):
        if self.unresolved == set(self.nx_graph.nodes):
            return True

        self.niter += 1
        for u in list(self.nx_graph.nodes):
            if u not in self.unresolved:
                self.__process_vertex(u)
        self.finalize_transformation()
        return False

    def transform(self, N):
        for _ in range(N):
            self.transform_single()

    def transform_until_saturated(self):
        while not self.transform_single():
            pass

    def get_niter_wo_complex(self):
        n_iter_wo_complex = math.inf
        for u in list(self.nx_graph.nodes):
            if u in self.unresolved:
                continue
            indegree = self.nx_graph.in_degree(u)
            outdegree = self.nx_graph.out_degree(u)
            if (indegree >= 2 and outdegree >= 2) or \
                    (indegree == 0 and outdegree >= 2) or \
                    (indegree >= 2 and outdegree == 0):
                n_iter_wo_complex = 0
                break
            nlen = self.node2len[u]
            in_indexes = [self.edge2index[e_in]
                          for e_in in self.nx_graph.in_edges(u, keys=True)]
            out_indexes = [self.edge2index[e_out]
                           for e_out in self.nx_graph.out_edges(u, keys=True)]
            if indegree == 1:
                index = in_indexes[0]
                fin_node = self.index2edge[index][0]
            elif outdegree == 1:
                index = out_indexes[0]
                fin_node = self.index2edge[index][1]
            else:
                assert False
            seq = self.edge2seq[index]
            n_iter_node = len(seq) - nlen
            if fin_node in self.unresolved:
                n_iter_node -= 1
            n_iter_wo_complex = min(n_iter_wo_complex, n_iter_node)
        n_iter_wo_complex = min(n_iter_wo_complex,
                                self.SATURATING_K - self.init_k - self.niter)
        # k + n + N == K => N = K - k - n
        return n_iter_wo_complex

    def _transform_simple_N(self, N):
        if N == 0:
            return
        for u in list(self.nx_graph.nodes):
            if u in self.unresolved:
                continue
            in_indexes = [self.edge2index[e_in]
                          for e_in in self.nx_graph.in_edges(u, keys=True)]
            out_indexes = [self.edge2index[e_out]
                           for e_out in self.nx_graph.out_edges(u, keys=True)]

            indegree = self.nx_graph.in_degree(u)
            outdegree = self.nx_graph.out_degree(u)
            nlen = self.node2len[u]
            if indegree == 0 and outdegree == 1:
                pass
            elif indegree == 1 and outdegree == 0:
                pass
            elif indegree == 1 and outdegree > 1:
                # simple 1-in vertex
                assert len(in_indexes) == 1
                in_index = in_indexes[0]
                in_seq = self.edge2seq[in_index]
                prefix = in_seq[-nlen-N:-nlen]
                for j in out_indexes:
                    assert self.edge2seq[j][:nlen] == in_seq[-nlen:]
                    self.edge2seq[j] = prefix + self.edge2seq[j]

            elif indegree > 1 and outdegree == 1:
                # simple 1-out vertex
                assert len(out_indexes) == 1
                out_index = out_indexes[0]
                out_seq = self.edge2seq[out_index]
                suffix = out_seq[nlen:nlen+N]
                for i in in_indexes:
                    assert self.edge2seq[i][-nlen:] == out_seq[:nlen]
                    self.edge2seq[i] += suffix
            else:
                assert False
            self.node2len[u] += N
        self.niter += N
        self.finalize_transformation()

    def finalize_transformation(self):
        collapsed_edges = []
        for edge in self.nx_graph.edges:
            index = self.edge2index[edge]
            seq = self.edge2seq[index]
            u, v, _ = edge
            if len(seq) == self.node2len[u] or len(seq) == self.node2len[v]:
                assert self.node2len[u] == self.node2len[v]
                assert u not in self.unresolved and v not in self.unresolved
                collapsed_edges.append(index)
        # remove collapsed edges
        [self.remove_edge(index=index) for index in collapsed_edges]

        self.nx_graph.remove_nodes_from(list(nx.isolates(self.nx_graph)))
        self._update_unresolved_vertices()
        self.assert_validity()

    def transform_single_fast(self):
        if self.unresolved == set(self.nx_graph.nodes):
            return True

        n_iter_wo_complex = self.get_niter_wo_complex()
        logger.info(f'iter={self.niter}, simple_iter={n_iter_wo_complex}')

        if n_iter_wo_complex > 0:
            self._transform_simple_N(N=n_iter_wo_complex)
        else:
            self.transform_single()
        return False

    def transform_fast_until_saturated(self):
        while self.init_k + self.niter < self.SATURATING_K and \
                not self.transform_single_fast():
            pass
        self.assert_validity()
        K = self.init_k+self.niter
        logger.info(f'Graph saturated, niter={self.niter}, K={K}')

    def estimate_lower_mult(self):
        mult = {edge: 1 for edge in self.index2edge}
        changed = True
        while changed:
            changed = False
            for u in self.nx_graph.nodes():
                in_indexes = \
                    [self.edge2index[e_in]
                     for e_in in self.nx_graph.in_edges(u, keys=True)]
                out_indexes = \
                    [self.edge2index[e_out]
                     for e_out in self.nx_graph.out_edges(u, keys=True)]
                indegree = self.nx_graph.in_degree(u)
                outdegree = self.nx_graph.out_degree(u)

                in_mult = sum(mult[edge] for edge in in_indexes)
                out_mult = sum(mult[edge] for edge in out_indexes)
                if indegree == 1 and in_mult < out_mult:
                    mult[in_indexes[0]] = out_mult
                    changed = True
                elif outdegree == 1 and out_mult < in_mult:
                    mult[out_indexes[0]] = in_mult
                    changed = True
        return mult

    def write_dot(self, outfile, reffn=None, compact=False, export_pdf=True,):
        if reffn is not None:
            ref = read_bio_seq(reffn)
            ref = compress_homopolymer(ref)
            A = ahocorasick.Automaton()
            # A = KeywordTree()
            # seq2index = {}
            for index, seq in self.edge2seq.items():
                seq = ''.join(seq)
                A.add_word(''.join(seq), index)
                # A.add(seq)
                # seq2index[seq] = index
            A.make_automaton()
            # A.finalize()
            print('made automaton')

            mult = defaultdict(lambda: [0, 0])
            for end_pos, index in A.iter(ref):
            # for seq, start_pos in A.search_all(ref):
                # index = seq2index[seq]
                mult[index][0] += 1
            for end_pos, index in A.iter(RC(ref)):
            # for seq, start_pos in A.search_all(RC(ref)):
                # index = seq2index[seq]
                mult[index][1] += 1

        if outfile[-3:] == 'dot':
            outfile = outfile[:-4]
        graph = nx.MultiDiGraph()
        for node in self.nx_graph.nodes():
            graph.add_node(node, label=f'{node} len={self.node2len[node]}')
        for edge in self.nx_graph.edges(keys=True):
            index = self.edge2index[edge]
            seq = self.edge2seq[index] if not compact else None
            seqlen = len(self.edge2seq[index])
            label = f'index={index}\nlen={seqlen}'
            if reffn is not None:
                # print(mult[index], mult_est[index])
                # assert mult[index] == 0 or mult[index] >= mult_est[index]
                label += f'\nmult_real={mult[index]}'
            graph.add_edge(*edge,
                           label=label,
                           seq=seq)
        dotfile = f'{outfile}.dot'
        nx.drawing.nx_pydot.write_dot(graph, dotfile)
        if export_pdf:
            pdffile = f'{outfile}.pdf'
            # https://stackoverflow.com/a/3516106
            cmd = ['dot', '-Tpdf', dotfile, '-o', pdffile]
            call(cmd)

    def align_ref2(self, ref):
        A = ahocorasick.Automaton()
        for index, seq in self.edge2seq.items():
            A.add_word(''.join(seq), index)
        A.make_automaton()
        print('made automaton')

        mult = defaultdict(lambda: [0, 0])
        for end_pos, index in A.iter(ref):
            start_pos = end_pos - len(self.edge2seq[index]) + 1
            mult[index][0] += 1
        for end_pos, index in A.iter(RC(ref)):
            start_pos = end_pos - len(self.edge2seq[index]) + 1
            mult[index][1] += 1
        return mult

    def align_ref(self, ref, MAX_TIP_LEN=200000):
        ref = list(ref)
        path = []
        coord = 0
        coords = defaultdict(list)
        for node in self.nx_graph.nodes():
            in_degree = self.nx_graph.in_degree(node)
            out_degree = self.nx_graph.out_degree(node)
            if in_degree == 0 and out_degree == 1:
                e_out = self.nx_graph.out_edges(node, keys=True)
                e_out = list(e_out)[0]
                end_nodelen = self.node2len[e_out[1]]
                out_index = self.edge2index[e_out]
                seq = self.edge2seq[out_index]
                alignment = edlib.align(seq, ref, k=0,
                                        task='locations',
                                        mode='HW')
                if len(alignment['locations']) == 1:
                    path.append(out_index)
                    st, en = alignment['locations'][0]
                    en += 1
                    print(st, en, en - st, len(seq))
                    coord = en - end_nodelen
                    coords[out_index].append((st, en))
                    break
        if len(path) == 0:
            return
        while True:
            last_index = path[-1]
            _, u, _ = self.index2edge[last_index]
            nodelen = self.node2len[u]
            out_edges = self.nx_graph.out_edges(u, keys=True)
            if len(out_edges) == 0:
                break
            extended = False
            for edge in out_edges:
                index = self.edge2index[edge]
                seq = self.edge2seq[index]
                c = seq[nodelen]
                if c == ref[coord+nodelen]:
                    # print(index)
                    refseq = ref[coord:coord+len(seq)]
                    seq = seq[:len(refseq)]
                    if refseq != seq:
                        print(path)
                        print(last_index, index)
                        # print(''.join(refseq))
                        # print(''.join(seq))
                        print(len(refseq))
                        print(len(seq))
                        print(coord)
                        for i, (a, b) in enumerate(zip(refseq, seq)):
                            if a != b:
                                print(i, a, b)
                    assert refseq == seq
                    path.append(index)
                    end_nodelen = self.node2len[edge[1]]
                    old_coord = coord
                    coord += len(seq) - end_nodelen
                    coords[index].append((old_coord, coord))
                    extended = True
                    break
            if not extended:
                break

        print(len(ref) - end_nodelen - coord)
        iscomplete = (len(ref) - end_nodelen - coord) < MAX_TIP_LEN
        return path, iscomplete, coords

    def _get_path(self, list_indexes):
        if len(list_indexes) == 0:
            return tuple()

        path = []
        fst_index = list_indexes[0]
        path += self.edge2seq[fst_index]
        for index in list_indexes[1:]:
            string = self.edge2seq[index]
            in_node = self.index2edge[index][0]
            len_node = self.node2len[in_node]
            assert path[-len_node:] == string[:len_node]
            path += string[len_node:]
        return tuple(path)

    def get_paths_for_mappings(self):
        paths = {}
        for r_id, mapping in self.idb_mappings.mappings.items():
            paths[r_id] = self._get_path(mapping)
        return paths

    def assert_validity_of_paths(self, reads):
        paths = self.get_paths_for_mappings()
        # print(set(paths.keys()) - set(reads.keys()))
        print(set(reads.keys()) - set(paths.keys()))
        # assert set(paths.keys()) == set(reads.keys())
        for r_id in paths:
            if r_id not in reads:
                continue
            path = ''.join(paths[r_id])
            read = reads[r_id]
            alignment = edlib.align(read, path, k=0, mode='HW')
            if alignment['editDistance'] == -1:
                print(r_id)
                print(len(path), len(read), self.idb_mappings.mappings[r_id])
                print(path)
                print(read)
                print(alignment)
                break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--dbg", required=True,
                        help="Directory with DBG output")
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--ref")
    params = parser.parse_args()

    params.dbg = expandpath(params.dbg)
    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    logfn = os.path.join(params.outdir, 'inc_k.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='centroFlye: inc_k')
    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    db_fn = os.path.join(params.dbg, 'graph.fasta')
    align_fn = os.path.join(params.dbg, 'alignments.txt')
    dbg_log_fn = os.path.join(params.dbg, 'dbg.log')
    with open(dbg_log_fn) as f:
        cmd = f.readline().strip().split(' ')
        i = 0
        while cmd[i] != '-k':
            i += 1
        k = int(cmd[i+1]) + 1
    logger.info(f'init k = {k}')
    logger.info(f'Reading DBG output from {params.dbg}')
    lpdb = PathMultiKGraph.fromDR(db_fn=db_fn, align_fn=align_fn, k=k)
    logger.info(f'# vertices = {nx.number_of_nodes(lpdb.nx_graph)}')
    logger.info(f'# edges = {nx.number_of_edges(lpdb.nx_graph)}')
    logger.info(f'Finished reading DBG output')
    logger.info(f'Starting increasing k')
    lpdb.transform_fast_until_saturated()
    logger.info(f'Finished increasing k')
    logger.info(f'# vertices = {nx.number_of_nodes(lpdb.nx_graph)}')
    logger.info(f'# edges = {nx.number_of_edges(lpdb.nx_graph)}')


    outac = os.path.join(params.outdir, f'active_connections.txt')
    logger.info(f'Active connections output to {outac}')
    with open(outac, 'w') as f:
        ac = lpdb.idb_mappings.get_active_connections()
        ac = sorted(list(ac))
        for i, j in ac:
            print(f'{i} {j}', file=f)

    outuniquedges = os.path.join(params.outdir, f'unique_edges.txt')
    logger.info(f'Unique edges output to {outuniquedges}')
    with open(outuniquedges, 'w') as f:
        for index in sorted(list(lpdb.unique_edges)):
            print(index, file=f)

    outdot = os.path.join(params.outdir, f'dbg_{k}-{lpdb.init_k+lpdb.niter}')
    logger.info(f'Writing final graph to {outdot}')
    lpdb.write_dot(outdot, compact=True, reffn=params.ref)
    logger.info(f'Finished writing final graph')


if __name__ == "__main__":
    main()


# class DB1:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 2, string='CCT')
#     nx_graph.add_edge(1, 2, string='GACT')
#     nx_graph.add_edge(2, 3, string='CTA')
#     nx_graph.add_edge(3, 4, string='TAG')
#     nx_graph.add_edge(3, 5, string='TAA')
#     nx_graph.add_edge(2, 4, string='CTGT')
#
#
# pg = PathMultiKGraph.fromDB(DB1(),
#                             string_set={},
#                             raw_mappings={0: ([(0, 2, 0), (2, 3, 0)], 0, 0),
#                                           1: ([(0, 2, 0), (2, 4, 0)], 0, 0)
#                                           })
#
# # graph starting vertex
#
#
# class DBStVertex:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 1, string='AAAAA')
#     nx_graph.add_edge(0, 2, string='AAACA')
#     nx_graph.add_edge(0, 3, string='AAA')
#
#
# pg = PathMultiKGraph.fromDB(DBStVertex(),
#                             string_set={},
#                             raw_mappings={0: ([(0, 1, 0)], 0, 0),
#                                           1: ([(0, 2, 0)], 0, 0),
#                                           2: ([(0, 3, 0)], 0, 0)})
# assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: [2]}
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == [((0, 1, 0), 0), ((4, 2, 0), 1)]
# assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: []}
#
#
# # graph ending vertex
# class DBEnVertex:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 3, string='AAAA')
#     nx_graph.add_edge(1, 3, string='AACA')
#     nx_graph.add_edge(2, 3, string='AAA')
#
#
# pg = PathMultiKGraph.fromDB(DBEnVertex(),
#                             string_set={},
#                             raw_mappings={0: ([(0, 3, 0)], 0, 0),
#                                           1: ([(1, 3, 0)], 0, 0),
#                                           2: ([(2, 3, 0)], 0, 0)})
# assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: [2]}
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == [((0, 3, 0), 0), ((1, 4, 0), 1)]
# assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: []}
#
#
# # graph 1-in >1-out
# class DB1inVertex:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 1, string='AACAG')
#     nx_graph.add_edge(1, 2, string='AGACC')
#     nx_graph.add_edge(1, 3, string='AGATT')
#     nx_graph.add_edge(1, 4, string='AGAGG')
#
#
# pg = PathMultiKGraph.fromDB(
#     DB1inVertex(), string_set={}, raw_mappings={})
# assert[pg.edge2seq[pg.edge2index[edge]]
#        for edge in pg.nx_graph.edges] == [list('AACAG'),
#                                           ['A', 'G', 'A', 'C', 'C'],
#                                           ['A', 'G', 'A', 'T', 'T'],
#                                           ['A', 'G', 'A', 'G', 'G']]
# pg.transform_single()
# assert[pg.edge2seq[pg.edge2index[edge]]
#        for edge in pg.nx_graph.edges] == [list('AACAG'),
#                                           ['C', 'A', 'G', 'A', 'C', 'C'],
#                                           ['C', 'A', 'G', 'A', 'T', 'T'],
#                                           ['C', 'A', 'G', 'A', 'G', 'G']]
#
#
# # graph >1-in 1-out
# class DB1outVertex:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 3, string='CCAGA')
#     nx_graph.add_edge(1, 3, string='TTAGA')
#     nx_graph.add_edge(2, 3, string='GGAGA')
#     nx_graph.add_edge(3, 4, string='GAAAA')
#
#
# pg = PathMultiKGraph.fromDB(
#     DB1outVertex(), string_set={}, raw_mappings={})
# assert sorted(
#     [pg.edge2seq[pg.edge2index[edge]]
#      for edge in pg.nx_graph.edges]) == sorted(
#     [list('CCAGA'),
#      list('TTAGA'),
#      list('GGAGA'),
#      list('GAAAA')])
# pg.transform_single()
# assert sorted(
#     [pg.edge2seq[pg.edge2index[edge]]
#      for edge in pg.nx_graph.edges]) == sorted(
#     [list('CCAGAA'),
#      list('TTAGAA'),
#      list('GGAGAA'),
#      list('GAAAA')])
#
#
# # graph with a complex vertex
#
# class DBComplexVertex1:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 2, string='ACAAA')
#     nx_graph.add_edge(1, 2, string='GGAAA')
#     nx_graph.add_edge(2, 3, string='AATGC')
#     nx_graph.add_edge(2, 4, string='AATT')
#
#
# pg = PathMultiKGraph.fromDB(
#     DBComplexVertex1(),
#     string_set={},
#     raw_mappings={0: ([(0, 2, 0), (2, 3, 0)], 0, 0),
#                   1: ([(0, 2, 0), (2, 4, 0)], 0, 0),
#                   2: ([(1, 2, 0), (2, 4, 0)], 0, 0)})
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#          [((0, 2, 0), 0), ((2, 3, 0), 1), ((2, 4, 0), 2), ((1, 2, 0), 3)]
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#     [((0, 5, 0), 0),
#      ((1, 8, 0), 3),
#      ((5, 3, 0), 1),
#      ((5, 8, 0), 4),
#      ((8, 4, 0), 2)]
# assert pg.edge2seq == \
#     {0: ['A', 'C', 'A', 'A', 'A'],
#      1: ['A', 'A', 'A', 'T', 'G', 'C'],
#      2: ['A', 'A', 'T', 'T'],
#      3: ['G', 'G', 'A', 'A', 'A', 'T'],
#      4: ['A', 'A', 'A', 'T']}
#
#
# # graph with a complex vertex
#
# class DBComplexVertex2:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 2, string='ACAAA')
#     nx_graph.add_edge(1, 2, string='GGAAA')
#     nx_graph.add_edge(2, 3, string='AATGC')
#     nx_graph.add_edge(2, 4, string='AATT')
#
#
# pg = PathMultiKGraph.fromDB(
#     DBComplexVertex2(),
#     string_set={},
#     raw_mappings={0: ([(0, 2, 0), (2, 3, 0)], 0, 0),
#                   1: ([(0, 2, 0), (2, 4, 0)], 0, 0),
#                   2: ([(1, 2, 0), (2, 4, 0)], 0, 0),
#                   3: ([(1, 2, 0), (2, 3, 0)], 0, 0)})
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#          [((0, 2, 0), 0), ((2, 3, 0), 1), ((2, 4, 0), 2), ((1, 2, 0), 3)]
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#     [((0, 5, 0), 0), ((1, 6, 0), 3), ((5, 7, 0), 4),
#      ((5, 8, 0), 5), ((6, 7, 0), 6),
#      ((6, 8, 0), 7), ((7, 3, 0), 1), ((8, 4, 0), 2)]
# assert pg.edge2seq == \
#     {0: ['A', 'C', 'A', 'A', 'A'],
#      1: list('AATGC'),
#      2: list('AATT'),
#      3: ['G', 'G', 'A', 'A', 'A'],
#      4: list('AAAT'),
#      5: list('AAAT'),
#      6: list('AAAT'),
#      7: list('AAAT')}
#
#
# # graph with a complex vertex
#
# class DBComplexVertex3:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 2, string='CCAC')
#     nx_graph.add_edge(1, 2, string='TTAC')
#     nx_graph.add_edge(2, 5, string='ACG')
#     nx_graph.add_edge(3, 5, string='AACG')
#     nx_graph.add_edge(5, 4, string='CGTA')
#     nx_graph.add_edge(5, 6, string='CGA')
#     nx_graph.add_edge(6, 7, string='GACC')
#     nx_graph.add_edge(6, 8, string='GATT')
#
#
# pg = PathMultiKGraph.fromDB(
#     DBComplexVertex3(),
#     string_set={},
#     raw_mappings={0: ([(0, 2, 0),
#                        (2, 5, 0),
#                        (5, 6, 0),
#                        (6, 7, 0)],
#                       0, 0),
#                   1: ([(1, 2, 0),
#                        (2, 5, 0),
#                        (5, 4, 0)],
#                       0, 0),
#                   2: ([(3, 5, 0),
#                        (5, 6, 0),
#                        (6, 8, 0)],
#                       0, 0)})
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#     [((0, 2, 0), 0), ((2, 5, 0), 1),
#      ((5, 4, 0), 3), ((5, 6, 0), 4), ((1, 2, 0), 2),
#      ((6, 7, 0), 6), ((6, 8, 0), 7), ((3, 5, 0), 5)]
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#     [((0, 2, 0), 0), ((2, 4, 0), 3),
#      ((2, 12, 0), 8), ((1, 2, 0), 2), ((3, 12, 0), 5),
#      ((12, 7, 0), 6), ((12, 8, 0), 7)]
# assert pg.edge2seq == \
#     {0: ['C', 'C', 'A', 'C', 'G'],
#      2: ['T', 'T', 'A', 'C', 'G'],
#      3: ['A', 'C', 'G', 'T', 'A'],
#      5: ['A', 'A', 'C', 'G', 'A'],
#      6: ['C', 'G', 'A', 'C', 'C'],
#      7: ['C', 'G', 'A', 'T', 'T'],
#      8: ['A', 'C', 'G', 'A']}
#
#
# # graph with a loop
#
# class DBLoop:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 1, string='CCAA')
#     nx_graph.add_edge(1, 1, string='AACAA')
#     nx_graph.add_edge(1, 2, string='AATT')
#
#
# pg = PathMultiKGraph.fromDB(
#     DBLoop(), string_set={},
#     raw_mappings={0: ([(0, 1, 0), (1, 1, 0)], 0, 0),
#                   1: ([(1, 1, 0), (1, 2, 0)], 0, 0)})
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == [((0, 2, 0), 0)]
# pg.edge2seq
# assert pg.edge2seq == \
#     {0: ['C', 'C', 'A', 'A', 'C', 'A', 'A', 'T', 'T']}
#
#
# # graph with a loop
#
# class DBLoop2:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 1, string='CCAA')
#     nx_graph.add_edge(1, 1, string='AACAA')
#     nx_graph.add_edge(1, 1, string='AAGAA')
#     nx_graph.add_edge(1, 2, string='AATT')
#
#
# pg = PathMultiKGraph.fromDB(
#     DBLoop2(), string_set={},
#     raw_mappings={0: ([(0, 1, 0), (1, 1, 0), (1, 2, 0)], 0, 0),
#                   1: ([(0, 1, 0), (1, 1, 1), (1, 2, 0)], 0, 0)})
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#          [((0, 3, 0), 0), ((3, 8, 0), 1), ((3, 8, 1), 2), ((8, 2, 0), 3)]
# pg.edge2seq
# assert pg.edge2seq == \
#     {0: ['C', 'C', 'A', 'A'],
#      1: ['C', 'A', 'A', 'C', 'A', 'A', 'T'],
#      2: ['C', 'A', 'A', 'G', 'A', 'A', 'T'],
#      3: ['A', 'A', 'T', 'T']}
#
#
# # graph with a loop
#
# class DBLoop3:
#     k = 3
#     nx_graph = nx.MultiDiGraph()
#     nx_graph.add_edge(0, 1, string='CCAA')
#     nx_graph.add_edge(1, 1, string='AACAA')
#     nx_graph.add_edge(1, 1, string='AAGAA')
#     nx_graph.add_edge(1, 2, string='AATT')
#     nx_graph.add_edge(3, 1, string='GGAA')
#     nx_graph.add_edge(1, 4, string='AARR')
#
#
# pg = PathMultiKGraph.fromDB(
#     DBLoop3(),
#     string_set={},
#     raw_mappings={0: ([(0, 1, 0),
#                        (1, 1, 0),
#                        (1, 2, 0)],
#                       0, 0),
#                   1: ([(3, 1, 0),
#                        (1, 1, 1),
#                        (1, 1, 1),
#                        (1, 4, 0)],
#                       0, 0),
#                   2: ([(3, 1, 0),
#                        (1, 2, 0)],
#                       0, 0)})
# pg.transform_single()
# assert [(edge, pg.edge2index[edge])
#         for edge in pg.nx_graph.edges] == \
#     [((0, 11, 0), 0), ((3, 8, 0), 5), ((7, 10, 0), 6),
#      ((7, 4, 0), 4), ((8, 10, 0), 7), ((8, 11, 0), 8),
#      ((10, 7, 0), 2), ((11, 2, 0), 3)]
# pg.edge2seq
# assert pg.edge2seq == \
#     {0: ['C', 'C', 'A', 'A', 'C', 'A', 'A', 'T'],
#      2: ['A', 'A', 'G', 'A', 'A'],
#      3: ['A', 'A', 'T', 'T'],
#      4: ['G', 'A', 'A', 'R', 'R'],
#      5: ['G', 'G', 'A', 'A'],
#      6: ['G', 'A', 'A', 'G'],
#      7: ['G', 'A', 'A', 'G'],
#      8: ['G', 'A', 'A', 'T']}
