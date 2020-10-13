# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging
from collections import defaultdict

import networkx as nx

from sequence_graph.path_graph import IDBMappings
from subprocess import call
from utils.various import fst_iterable


logger = logging.getLogger("centroFlye.sequence_graph.path_graph_multik")


class PathMultiKGraph:
    def __init__(self, nx_graph,
                 edge2seq, edge2index, index2edge, node2len,
                 max_edge_index, max_node_index,
                 idb_mappings,
                 init_k):
        self.nx_graph = nx_graph
        self.edge2seq = edge2seq
        self.edge2index = edge2index
        self.index2edge = index2edge
        self.node2len = node2len
        self.max_edge_index = max_edge_index
        self.max_node_index = max_node_index
        self.idb_mappings = idb_mappings
        self.init_k = init_k

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
                all_ac = self.idb_mappings.get_active_connections()
                paired_indexes = set()
                for e_in in in_indexes:
                    for e_out in out_indexes:
                        if (e_in, e_out) in all_ac:
                            paired_indexes.add(e_in)
                            paired_indexes.add(e_out)
                indexes = (in_indexes | out_indexes) - paired_indexes

                if len(indexes):
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
                        break
                for edge in self.nx_graph.out_edges(u, keys=True):
                    index = self.edge2index[edge]
                    seq = self.edge2seq[index]
                    v = edge[1]
                    if v in self.unresolved:
                        continue
                    if self.node2len[v] + 1 == len(seq):
                        new.add(v)
                        break
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

            all_ac = self.idb_mappings.get_active_connections()
            ac_s2e = defaultdict(set)
            ac_e2s = defaultdict(set)
            for e_in in in_indexes:
                for e_out in out_indexes:
                    if (e_in, e_out) in all_ac:
                        ac_s2e[e_in].add(e_out)
                        ac_e2s[e_out].add(e_in)

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
        self.niter += 1
        for u in list(self.nx_graph.nodes):
            if u not in self.unresolved:
                self.__process_vertex(u)

        collapsed_edges = []
        for edge in self.nx_graph.edges:
            index = self.edge2index[edge]
            seq = self.edge2seq[index]
            u, v, _ = edge
            if u not in self.unresolved and \
                    v not in self.unresolved and \
                    len(seq) == self.node2len[u] == self.node2len[v]:
                collapsed_edges.append(index)
        # remove collapsed edges
        [self.remove_edge(index=index) for index in collapsed_edges]

        self._update_unresolved_vertices()
        self.assert_validity()

    def transform(self, N):
        for _ in range(N):
            self.transform_single()

    def write_dot(self, outfile, compact=False, export_pdf=True):
        if outfile[-3:] == 'dot':
            outfile = outfile[:-4]
        graph = nx.MultiDiGraph()
        for node in self.nx_graph.nodes():
            graph.add_node(node, label=f'{node} len={self.node2len[node]}')
        for edge in self.nx_graph.edges(keys=True):
            index = self.edge2index[edge]
            seq = self.edge2seq[index] if not compact else None
            seqlen = len(self.edge2seq[index])
            graph.add_edge(*edge,
                           label=f'index={index}\nlen={seqlen}',
                           seq=seq)
        dotfile = f'{outfile}.dot'
        nx.drawing.nx_pydot.write_dot(graph, dotfile)
        if export_pdf:
            pdffile = f'{outfile}.pdf'
            print(pdffile)
            # https://stackoverflow.com/a/3516106
            cmd = ['dot', '-Tpdf', dotfile, '-o', pdffile]
            call(cmd)


class DB1:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 2, string='CCT')
    nx_graph.add_edge(1, 2, string='GACT')
    nx_graph.add_edge(2, 3, string='CTA')
    nx_graph.add_edge(3, 4, string='TAG')
    nx_graph.add_edge(3, 5, string='TAA')
    nx_graph.add_edge(2, 4, string='CTGT')


pg = PathMultiKGraph.fromDB(DB1(),
                            string_set={},
                            raw_mappings={0: ([(0, 2, 0), (2, 3, 0)], 0, 0),
                                          1: ([(0, 2, 0), (2, 4, 0)], 0, 0)
                                          })

# graph starting vertex


class DBStVertex:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 1, string='AAAAA')
    nx_graph.add_edge(0, 2, string='AAACA')
    nx_graph.add_edge(0, 3, string='AAA')


pg = PathMultiKGraph.fromDB(DBStVertex(),
                            string_set={},
                            raw_mappings={0: ([(0, 1, 0)], 0, 0),
                                          1: ([(0, 2, 0)], 0, 0),
                                          2: ([(0, 3, 0)], 0, 0)})
assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: [2]}
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == [((0, 1, 0), 0), ((4, 2, 0), 1)]
assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: []}


# graph ending vertex
class DBEnVertex:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 3, string='AAAA')
    nx_graph.add_edge(1, 3, string='AACA')
    nx_graph.add_edge(2, 3, string='AAA')


pg = PathMultiKGraph.fromDB(DBEnVertex(),
                            string_set={},
                            raw_mappings={0: ([(0, 3, 0)], 0, 0),
                                          1: ([(1, 3, 0)], 0, 0),
                                          2: ([(2, 3, 0)], 0, 0)})
assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: [2]}
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == [((0, 3, 0), 0), ((1, 4, 0), 1)]
assert pg.idb_mappings.mappings == {0: [0], 1: [1], 2: []}


# graph 1-in >1-out
class DB1inVertex:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 1, string='AACAG')
    nx_graph.add_edge(1, 2, string='AGACC')
    nx_graph.add_edge(1, 3, string='AGATT')
    nx_graph.add_edge(1, 4, string='AGAGG')


pg = PathMultiKGraph.fromDB(
    DB1inVertex(), string_set={}, raw_mappings={})
assert[pg.edge2seq[pg.edge2index[edge]]
       for edge in pg.nx_graph.edges] == [list('AACAG'),
                                          ['A', 'G', 'A', 'C', 'C'],
                                          ['A', 'G', 'A', 'T', 'T'],
                                          ['A', 'G', 'A', 'G', 'G']]
pg.transform_single()
assert[pg.edge2seq[pg.edge2index[edge]]
       for edge in pg.nx_graph.edges] == [list('AACAG'),
                                          ['C', 'A', 'G', 'A', 'C', 'C'],
                                          ['C', 'A', 'G', 'A', 'T', 'T'],
                                          ['C', 'A', 'G', 'A', 'G', 'G']]


# graph >1-in 1-out
class DB1outVertex:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 3, string='CCAGA')
    nx_graph.add_edge(1, 3, string='TTAGA')
    nx_graph.add_edge(2, 3, string='GGAGA')
    nx_graph.add_edge(3, 4, string='GAAAA')


pg = PathMultiKGraph.fromDB(
    DB1outVertex(), string_set={}, raw_mappings={})
assert sorted(
    [pg.edge2seq[pg.edge2index[edge]]
     for edge in pg.nx_graph.edges]) == sorted(
    [list('CCAGA'),
     list('TTAGA'),
     list('GGAGA'),
     list('GAAAA')])
pg.transform_single()
assert sorted(
    [pg.edge2seq[pg.edge2index[edge]]
     for edge in pg.nx_graph.edges]) == sorted(
    [list('CCAGAA'),
     list('TTAGAA'),
     list('GGAGAA'),
     list('GAAAA')])


# graph with a complex vertex

class DBComplexVertex1:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 2, string='ACAAA')
    nx_graph.add_edge(1, 2, string='GGAAA')
    nx_graph.add_edge(2, 3, string='AATGC')
    nx_graph.add_edge(2, 4, string='AATT')


pg = PathMultiKGraph.fromDB(
    DBComplexVertex1(),
    string_set={},
    raw_mappings={0: ([(0, 2, 0), (2, 3, 0)], 0, 0),
                  1: ([(0, 2, 0), (2, 4, 0)], 0, 0),
                  2: ([(1, 2, 0), (2, 4, 0)], 0, 0)})
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
         [((0, 2, 0), 0), ((2, 3, 0), 1), ((2, 4, 0), 2), ((1, 2, 0), 3)]
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
    [((0, 5, 0), 0),
     ((1, 8, 0), 3),
     ((5, 3, 0), 1),
     ((5, 8, 0), 4),
     ((8, 4, 0), 2)]
assert pg.edge2seq == \
    {0: ['A', 'C', 'A', 'A', 'A'],
     1: ['A', 'A', 'A', 'T', 'G', 'C'],
     2: ['A', 'A', 'T', 'T'],
     3: ['G', 'G', 'A', 'A', 'A', 'T'],
     4: ['A', 'A', 'A', 'T']}


# graph with a complex vertex

class DBComplexVertex2:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 2, string='ACAAA')
    nx_graph.add_edge(1, 2, string='GGAAA')
    nx_graph.add_edge(2, 3, string='AATGC')
    nx_graph.add_edge(2, 4, string='AATT')


pg = PathMultiKGraph.fromDB(
    DBComplexVertex2(),
    string_set={},
    raw_mappings={0: ([(0, 2, 0), (2, 3, 0)], 0, 0),
                  1: ([(0, 2, 0), (2, 4, 0)], 0, 0),
                  2: ([(1, 2, 0), (2, 4, 0)], 0, 0),
                  3: ([(1, 2, 0), (2, 3, 0)], 0, 0)})
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
         [((0, 2, 0), 0), ((2, 3, 0), 1), ((2, 4, 0), 2), ((1, 2, 0), 3)]
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
    [((0, 5, 0), 0), ((1, 6, 0), 3), ((5, 7, 0), 4),
     ((5, 8, 0), 5), ((6, 7, 0), 6),
     ((6, 8, 0), 7), ((7, 3, 0), 1), ((8, 4, 0), 2)]
assert pg.edge2seq == \
    {0: ['A', 'C', 'A', 'A', 'A'],
     1: list('AATGC'),
     2: list('AATT'),
     3: ['G', 'G', 'A', 'A', 'A'],
     4: list('AAAT'),
     5: list('AAAT'),
     6: list('AAAT'),
     7: list('AAAT')}


# graph with a complex vertex

class DBComplexVertex3:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 2, string='CCAC')
    nx_graph.add_edge(1, 2, string='TTAC')
    nx_graph.add_edge(2, 5, string='ACG')
    nx_graph.add_edge(3, 5, string='AACG')
    nx_graph.add_edge(5, 4, string='CGTA')
    nx_graph.add_edge(5, 6, string='CGA')
    nx_graph.add_edge(6, 7, string='GACC')
    nx_graph.add_edge(6, 8, string='GATT')


pg = PathMultiKGraph.fromDB(
    DBComplexVertex3(),
    string_set={},
    raw_mappings={0: ([(0, 2, 0),
                       (2, 5, 0),
                       (5, 6, 0),
                       (6, 7, 0)],
                      0, 0),
                  1: ([(1, 2, 0),
                       (2, 5, 0),
                       (5, 4, 0)],
                      0, 0),
                  2: ([(3, 5, 0),
                       (5, 6, 0),
                       (6, 8, 0)],
                      0, 0)})
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
    [((0, 2, 0), 0), ((2, 5, 0), 1),
     ((5, 4, 0), 3), ((5, 6, 0), 4), ((1, 2, 0), 2),
     ((6, 7, 0), 6), ((6, 8, 0), 7), ((3, 5, 0), 5)]
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
    [((0, 2, 0), 0), ((2, 4, 0), 3),
     ((2, 12, 0), 8), ((1, 2, 0), 2), ((3, 12, 0), 5),
     ((12, 7, 0), 6), ((12, 8, 0), 7)]
assert pg.edge2seq == \
    {0: ['C', 'C', 'A', 'C', 'G'],
     2: ['T', 'T', 'A', 'C', 'G'],
     3: ['A', 'C', 'G', 'T', 'A'],
     5: ['A', 'A', 'C', 'G', 'A'],
     6: ['C', 'G', 'A', 'C', 'C'],
     7: ['C', 'G', 'A', 'T', 'T'],
     8: ['A', 'C', 'G', 'A']}


# graph with a loop

class DBLoop:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 1, string='CCAA')
    nx_graph.add_edge(1, 1, string='AACAA')
    nx_graph.add_edge(1, 2, string='AATT')


pg = PathMultiKGraph.fromDB(
    DBLoop(), string_set={},
    raw_mappings={0: ([(0, 1, 0), (1, 1, 0)], 0, 0),
                  1: ([(1, 1, 0), (1, 2, 0)], 0, 0)})
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == [((0, 2, 0), 0)]
pg.edge2seq
assert pg.edge2seq == \
    {0: ['C', 'C', 'A', 'A', 'C', 'A', 'A', 'T', 'T']}


# graph with a loop

class DBLoop2:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 1, string='CCAA')
    nx_graph.add_edge(1, 1, string='AACAA')
    nx_graph.add_edge(1, 1, string='AAGAA')
    nx_graph.add_edge(1, 2, string='AATT')


pg = PathMultiKGraph.fromDB(
    DBLoop2(), string_set={},
    raw_mappings={0: ([(0, 1, 0), (1, 1, 0), (1, 2, 0)], 0, 0),
                  1: ([(0, 1, 0), (1, 1, 1), (1, 2, 0)], 0, 0)})
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
         [((0, 3, 0), 0), ((3, 8, 0), 1), ((3, 8, 1), 2), ((8, 2, 0), 3)]
pg.edge2seq
assert pg.edge2seq == \
    {0: ['C', 'C', 'A', 'A'],
     1: ['C', 'A', 'A', 'C', 'A', 'A', 'T'],
     2: ['C', 'A', 'A', 'G', 'A', 'A', 'T'],
     3: ['A', 'A', 'T', 'T']}


# graph with a loop

class DBLoop3:
    k = 3
    nx_graph = nx.MultiDiGraph()
    nx_graph.add_edge(0, 1, string='CCAA')
    nx_graph.add_edge(1, 1, string='AACAA')
    nx_graph.add_edge(1, 1, string='AAGAA')
    nx_graph.add_edge(1, 2, string='AATT')
    nx_graph.add_edge(3, 1, string='GGAA')
    nx_graph.add_edge(1, 4, string='AARR')


pg = PathMultiKGraph.fromDB(
    DBLoop3(),
    string_set={},
    raw_mappings={0: ([(0, 1, 0),
                       (1, 1, 0),
                       (1, 2, 0)],
                      0, 0),
                  1: ([(3, 1, 0),
                       (1, 1, 1),
                       (1, 1, 1),
                       (1, 4, 0)],
                      0, 0),
                  2: ([(3, 1, 0),
                       (1, 2, 0)],
                      0, 0)})
pg.transform_single()
assert [(edge, pg.edge2index[edge])
        for edge in pg.nx_graph.edges] == \
    [((0, 11, 0), 0), ((3, 8, 0), 5), ((7, 10, 0), 6),
     ((7, 4, 0), 4), ((8, 10, 0), 7), ((8, 11, 0), 8),
     ((10, 7, 0), 2), ((11, 2, 0), 3)]
pg.edge2seq
assert pg.edge2seq == \
    {0: ['C', 'C', 'A', 'A', 'C', 'A', 'A', 'T'],
     2: ['A', 'A', 'G', 'A', 'A'],
     3: ['A', 'A', 'T', 'T'],
     4: ['G', 'A', 'A', 'R', 'R'],
     5: ['G', 'G', 'A', 'A'],
     6: ['G', 'A', 'A', 'G'],
     7: ['G', 'A', 'A', 'G'],
     8: ['G', 'A', 'A', 'T']}
