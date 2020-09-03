# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging
from collections import Counter, defaultdict, deque
import itertools
import os

import networkx as nx
import numpy as np

from config.config import config
from sequence_graph.db_graph import DeBruijnGraph
from sequence_graph.db_graph_3col import DeBruijnGraph3Color
from utils.various import fst_iterable, filter_sublsts_n2_dict
from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.sequence_graph.path_graph")


class IDBMappings:
    def __init__(self, mappings, paths):
        self.mappings = mappings
        self.paths = paths
        self._update_active_connections()

    def _update_active_connections(self, min_con=1):
        self.mappings = {s_id: mapping
                         for s_id, mapping in self.mappings.items()
                         if len(mapping) > 1}
        connections = Counter()
        self.pairs_coords = defaultdict(list)
        for s_id, mapping in self.mappings.items():
            it1, it2 = itertools.tee(mapping)
            next(it2, None)
            for i, p in enumerate(zip(it1, it2)):
                connections[p] += 1
                self.pairs_coords[p].append((s_id, i))
        for p in self.pairs_coords:
            self.pairs_coords[p].sort(reverse=True)
        self.active_connections = set(k for k, v in connections.items()
                                      if v >= min_con)

    @classmethod
    def from_graph(cls, graph, string_set,
                   raw_mappings=None,
                   neutral_symbs=None):
        if neutral_symbs is None:
            neutral_symbs = set()
        if raw_mappings is None:
            raw_mappings = graph.map_strings(string_set,
                                             only_unique_paths=True,
                                             neutral_symbs=neutral_symbs)
        paths = {r_id: graph.get_path(edge_path, st, en)
                 for r_id, (edge_path, st, en) in raw_mappings.items()}

        unfilt_mappings = {}
        for s_id, (edge_list, _, _) in raw_mappings.items():
            if len(edge_list) == 0:
                continue
            mapping = [graph.nx_graph.get_edge_data(*edge)['edge_index']
                       for edge in edge_list]
            unfilt_mappings[s_id] = mapping

        mappings = filter_sublsts_n2_dict(unfilt_mappings)
        mappings = {k: deque(v) for k, v in mappings.items()}

        return cls(mappings=mappings, paths=paths)

    def remove_edge(self, edge_index):
        for s_id, mapping in self.mappings.items():
            mapping = deque(filter(lambda a: a != edge_index, mapping))
            self.mappings[s_id] = mapping
        self._update_active_connections()

    def remove_edges(self, edge_indexes):
        edge_indexes = set(edge_indexes)
        for s_id, mapping in self.mappings.items():
            mapping = deque(filter(lambda a: a not in edge_indexes, mapping))
            self.mappings[s_id] = mapping
        self._update_active_connections()

    def insert_edge(self, edge_to_insert, edge_pair):
        for (s_id, i) in self.pairs_coords[edge_pair]:
            self.mappings[s_id].insert(i+1, edge_to_insert)
        self._update_active_connections()

    def insert_edges(self, edge_pairs_edges_to_insert):
        for ep, edge in edge_pairs_edges_to_insert:
            self.insert_edge(edge, ep)

    def extend_nonbranching_mappings(self, db):
        def get_inedges(edge_index):
            s, _, _ = db.edge_index2edge[edge_index]
            out_degree = db.nx_graph.out_degree(s)
            inedges = list(db.nx_graph.in_edges(s, keys=True))
            inedges = [db.nx_graph.get_edge_data(*edge)['edge_index']
                       for edge in inedges]
            return out_degree, inedges

        for s_id, mapping in self.mappings.items():
            if len(mapping) == 0:
                continue
            out_degree, inedges = get_inedges(mapping[0])
            while out_degree == 1 and len(inedges) == 1:
                mapping.insert(0, inedges[0])
                out_degree, inedges = get_inedges(mapping[0])
        self._update_active_connections()

    def get_mapping_lens(self, db):
        lens = {}
        for s_id, mapping in self.mappings.items():
            edge_path = [db.edge_index2edge[edge_index]
                         for edge_index in mapping]
            path = db.get_path(edge_path, respect_cycled=False)
            lens[s_id] = len(path)
        return lens

    def remove_connections(self, removed_connections):
        coords = defaultdict(lambda: [0])
        for edge_pair in removed_connections:
            for (s_id, i) in self.pairs_coords[edge_pair]:
                coords[s_id].append(i)
        for s_id, coords_s_id in coords.items():
            coords_s_id.append(len(self.mappings[s_id]))
            for i, (s, e) in enumerate(zip(coords_s_id[::2],
                                           coords_s_id[1::2])):
                self.mappings[f'{s_id}_{i}'] = \
                    deque(list(self.mappings[s_id])[s:e])
            del self.mappings[s_id]
        self._update_active_connections()


class PathDeBruijnGraph(DeBruijnGraph):
    def __init__(self, db, idbmappings):
        self.__dict__.update(db.__dict__)
        self.idbmappings = idbmappings

    @classmethod
    def __increase_k_by_one(cls, db, idbmappings):
        tr_nx_db = nx.MultiDiGraph()
        nodeindex2label = {}
        nodelabel2index = {}

        k = db.k
        node_cnt = 0

        eindex2st = {}
        eindex2en = {}
        max_edge_index = 0
        edge2node = []

        for st, en, key, edge_data in db.nx_graph.edges(keys=True, data=True):
            max_edge_index = max(max_edge_index, edge_data[cls.edge_index])
            if k < len(edge_data[cls.string]):
                pref = edge_data[cls.string][:k]
                suff = edge_data[cls.string][-k:]
                nodeindex2label[node_cnt] = pref
                nodeindex2label[node_cnt+1] = suff
                nodelabel2index[pref] = node_cnt
                nodelabel2index[suff] = node_cnt+1

                length = edge_data[cls.length]-1
                label = DeBruijnGraph._generate_label(
                    {DeBruijnGraph.length: length,
                     DeBruijnGraph.coverage: np.mean(edge_data[cls.coverage]),
                     DeBruijnGraph.edge_index: edge_data[cls.edge_index]}
                )
                edge_index = edge_data[cls.edge_index]
                cov = edge_data[cls.coverage][:-1]
                tr_nx_db.add_edge(node_cnt, node_cnt+1,
                                  string=edge_data[cls.string],
                                  coverage=cov,
                                  label=label,
                                  length=length,
                                  color=edge_data[cls.color],
                                  edge_index=edge_index)
                eindex2st[edge_index] = node_cnt
                eindex2en[edge_index] = node_cnt + 1
                node_cnt += 2
            else:
                label = edge_data[cls.string]
                nodeindex2label[node_cnt] = label
                nodelabel2index[label] = node_cnt
                tr_nx_db.add_node(node_cnt)
                edge_index = edge_data[cls.edge_index]
                edge2node.append(edge_index)
                eindex2st[edge_index] = node_cnt
                eindex2en[edge_index] = node_cnt
                node_cnt += 1

        insertions = []
        removed_connections = []
        for node in db.nodeindex2label:
            in_degree = db.nx_graph.in_degree(node)
            out_degree = db.nx_graph.out_degree(node)
            simple_vertex = in_degree <= 1 or out_degree <= 1
            for _, _, _, st_data in db.nx_graph.in_edges(node,
                                                         data=True,
                                                         keys=True):
                st_edge_index = st_data[cls.edge_index]
                st = eindex2en[st_edge_index]
                for _, _, _, en_data in db.nx_graph.out_edges(node,
                                                              data=True,
                                                              keys=True):
                    en_edge_index = en_data[cls.edge_index]
                    en = eindex2st[en_edge_index]
                    connection = (st_edge_index, en_edge_index)
                    is_active = connection in idbmappings.active_connections
                    if simple_vertex or is_active:
                        max_edge_index += 1
                        last_chr = nodeindex2label[en][-1]
                        string = nodeindex2label[st] + (last_chr,)
                        label_params = {'length': 1,
                                        'coverage': 1.00,
                                        'edge_index': max_edge_index}
                        label = DeBruijnGraph._generate_label(label_params)
                        insertions.append(((st_edge_index, en_edge_index),
                                           max_edge_index))
                        tr_nx_db.add_edge(st, en,
                                          string=string,
                                          coverage=[1],
                                          label=label,
                                          length=1,
                                          color='black',
                                          edge_index=max_edge_index)
                    else:
                        removed_connections.append(connection)

        mapping_lens = idbmappings.get_mapping_lens(db)
        db = DeBruijnGraph(nx_graph=tr_nx_db,
                           nodeindex2label=nodeindex2label,
                           nodelabel2index=nodelabel2index,
                           k=k+1, collapse=False)

        idbmappings.remove_connections(removed_connections)
        idbmappings.insert_edges(insertions)
        idbmappings.remove_edges(edge2node)
        idbmappings.extend_nonbranching_mappings(db)

        edges_to_remove = db.collapse_nonbranching_paths()
        idbmappings.remove_edges(edges_to_remove)

        new_mapping_lens = idbmappings.get_mapping_lens(db)
        for s_id, new_len in new_mapping_lens.items():
            if s_id in mapping_lens:
                old_len = mapping_lens[s_id]
                if old_len > new_len:
                    logger.error(s_id, old_len, new_len)
                assert old_len <= new_len

        return db, idbmappings

    @classmethod
    def from_db(cls, db, string_set, k,
                remap_iter=config['path_db']['remap_iter'],
                neutral_symbs=None,
                mappings=None,
                outdir=None,
                assembly=None):
        if neutral_symbs is None:
            neutral_symbs = set()
        logger.info('Constructing a path graph')
        niter = k - db.k
        logger.info(f'Starting from DB(k) with k = {db.k}')
        logger.info(f'Goal: k = {k}. # iterations = {niter}')
        for i in range(niter):
            if i % remap_iter == 0:
                idbmappings = \
                    IDBMappings.from_graph(graph=db,
                                           string_set=string_set,
                                           raw_mappings=mappings,
                                           neutral_symbs=neutral_symbs)
            logger.info(f'Iteration #{i+1}, k = {db.k+1}')
            db, idbmappings = cls.__increase_k_by_one(db, idbmappings)

        obj = cls(db=db, idbmappings=idbmappings)

        if outdir is not None:
            smart_makedirs(outdir)
            dot_file = os.path.join(outdir, f'db_k{k}.dot')
            db.write_dot(outfile=dot_file, export_pdf=False)

            dot_compact_file = os.path.join(outdir, f'db_k{k}_compact.dot')
            db.write_dot(outfile=dot_compact_file,
                         export_pdf=True,
                         compact=True)
            pickle_file = os.path.join(outdir, f'db_k{k}.pickle')
            db.pickle_dump(pickle_file)
            if assembly is not None:
                DeBruijnGraph3Color.from_read_db_and_assembly(gr_reads=db,
                                                              assembly=assembly,
                                                              outdir=outdir)
        return obj

    @classmethod
    def from_mono_db(cls, db, monostring_set, k,
                     mappings=None,
                     outdir=None,
                     assembly=None):
        monostring = fst_iterable(monostring_set.values())
        neutral_symbs = set([monostring.gap_symb])
        return cls.from_db(db=db,
                           string_set=monostring_set,
                           k=k,
                           neutral_symbs=neutral_symbs,
                           mappings=mappings,
                           outdir=outdir,
                           assembly=assembly)

    def get_paths(self):
        return self.idbmappings.paths
