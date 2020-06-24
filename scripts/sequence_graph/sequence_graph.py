# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from abc import ABC, abstractmethod
from collections import defaultdict, namedtuple
import logging
import pickle
from subprocess import Popen

import networkx as nx

from utils.various import fst_iterable

logger = logging.getLogger("centroFlye.sequence_graph.sequence_graph")


class SequenceGraph(ABC):
    label = 'label'
    color = 'color'
    string = 'string'
    length = 'length'

    def __init__(self,
                 nx_graph=None,
                 nodeindex2label=None,
                 nodelabel2index=None,
                 collapse=True):
        all_None = all((nx_graph is None,
                        nodeindex2label is None,
                        nodelabel2index is None))
        all_notNone = all((nx_graph is not None,
                           nodeindex2label is not None,
                           nodelabel2index is not None))
        assert all_None or all_notNone

        if nx_graph is None:  # and thus all parameters are None
            self.nx_graph = nx.MultiDiGraph()
            self.nodeindex2label = {}
            self.nodelabel2index = {}
        else:
            self.nx_graph = nx_graph
            self.nodeindex2label = nodeindex2label
            self.nodelabel2index = nodelabel2index
            self._assert_nx_graph_validity()

        self.db_index = None
        if collapse:
            self.collapse_nonbranching_paths()

    @classmethod
    def from_pickle(cls, fn):
        with open(fn, 'rb') as f:
            db_graph = pickle.load(f)
        return db_graph

    @classmethod
    @abstractmethod
    def _generate_label(par_dict):
        pass

    @abstractmethod
    def _add_edge(self, node, color, string,
                  in_node, out_node,
                  in_data, out_data,
                  edge_len):
        pass

    def _assert_nx_graph_validity(self):
        for n_ind1, n_ind2, string in self.nx_graph.edges.data(self.string):
            n_label1 = self.nodeindex2label[n_ind1]
            n_label2 = self.nodeindex2label[n_ind2]
            assert string[:len(n_label1)] == n_label1
            assert string[-len(n_label2):] == n_label2
        assert len(self.nx_graph.nodes) == len(self.nodeindex2label)
        assert len(self.nx_graph.nodes) == len(self.nodelabel2index)

    def collapse_nonbranching_paths(self):
        def node_on_nonbranching_path(node):
            if self.nx_graph.in_degree(node) != 1 or \
               self.nx_graph.out_degree(node) != 1:
                return False

            in_edge = fst_iterable(self.nx_graph.in_edges(node, keys=True))
            out_edge = fst_iterable(self.nx_graph.out_edges(node, keys=True))

            in_edge_color = self.nx_graph.get_edge_data(*in_edge)[self.color]
            out_edge_color = self.nx_graph.get_edge_data(*out_edge)[self.color]

            # if in_edge == out_edge this is a loop and should not be removed
            return in_edge != out_edge and in_edge_color == out_edge_color

        nodes = list(self.nx_graph)
        for node in nodes:
            if node_on_nonbranching_path(node):
                in_edge = fst_iterable(self.nx_graph.in_edges(node, keys=True))
                out_edge = fst_iterable(self.nx_graph.out_edges(node,
                                                                keys=True))
                in_node, out_node = in_edge[0], out_edge[1]

                in_data = self.nx_graph.get_edge_data(*in_edge)
                out_data = self.nx_graph.get_edge_data(*out_edge)

                in_color = in_data[self.color]
                out_color = out_data[self.color]
                assert in_color == out_color
                color = in_color

                in_string = in_data[self.string]
                out_string = out_data[self.string]
                len_node = len(self.nodeindex2label[node])
                string = in_string + out_string[len_node:]
                string = tuple(string)
                edge_len = len(string) - len_node

                self._add_edge(node=node, color=color, string=string,
                               in_node=in_node, out_node=out_node,
                               in_data=in_data, out_data=out_data,
                               edge_len=edge_len)

                self.nx_graph.remove_node(node)
                label = self.nodeindex2label[node]
                del self.nodeindex2label[node]
                del self.nodelabel2index[label]

        self._assert_nx_graph_validity()
        self.index_edges()

    def get_path(self, list_edges):
        for e1, e2 in zip(list_edges[:-1], list_edges[1:]):
            assert e1[1] == e2[0]

        path = []
        fst_edge = list_edges[0]
        fst_edge_data = self.nx_graph.get_edge_data(*fst_edge)
        path += fst_edge_data[self.string]
        for edge in list_edges[1:]:
            edge_data = self.nx_graph.get_edge_data(*edge)
            edge_string = edge_data[self.string]
            in_node = edge[0]
            len_node = len(self.nodeindex2label[in_node])
            assert tuple(path[-len_node:]) == edge_string[:len_node]
            path += edge_string[len_node:]

        # process cycled path
        fst_node = fst_edge[0]
        lst_edge = list_edges[-1]
        lst_node = lst_edge[1]
        if fst_node == lst_node:
            len_node = len(self.nodeindex2label[fst_node])
            path = path[:-len_node]
        return tuple(path)

    def index_edges(self):
        if self.db_index is not None:
            return self.db_index
        all_index = {}

        # in case of db this is just min_k == max_k == k
        label_lens = [len(label) for label in self.nodelabel2index.keys()]
        min_k = 1 + min(label_lens)
        max_k = 1 + max(label_lens)
        for k in range(min_k, max_k+1):
            index = defaultdict(list)
            for e_ind, edge in enumerate(self.nx_graph.edges(keys=True)):
                e_data = self.nx_graph.get_edge_data(*edge)
                e_str = e_data[self.string]
                for i in range(len(e_str)-k+1):
                    kmer = e_str[i:i+k]
                    index[kmer].append((e_ind, i))
            index_unique = {kmer: pos[0]
                            for kmer, pos in index.items()
                            if len(pos) == 1}
            all_index[k] = index_unique
        self.db_index = all_index
        return all_index

    def get_contigs(self):
        # this function does not respect color of edges
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
                    # all outpaths for edges AFTER edge are computed
                    for i, e in enumerate(outpaths[edge][1:]):
                        outpaths[e] = outpaths[edge][i+1:]
            return outpaths

        outpaths = get_longest_valid_outpaths(self.nx_graph)

        rev_graph = self.nx_graph.reverse()
        rev_inpaths = get_longest_valid_outpaths(rev_graph)
        # need to reverse rev_inpaths
        inpaths = {}
        for rev_edge in rev_inpaths:
            edge = (rev_edge[1], rev_edge[0], rev_edge[2])
            inpaths[edge] = rev_inpaths[rev_edge][::-1]
            inpaths[edge] = [(e[1], e[0], e[2]) for e in inpaths[edge]]

        # unite inpaths and outpaths
        valid_paths = []
        for edge in outpaths:
            valid_path = inpaths[edge]
            seen_edges = set(inpaths[edge])
            for e in outpaths[edge][1:]:
                if e not in seen_edges:
                    valid_path.append(e)
                    seen_edges.add(e)
                else:
                    # cycle
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

        # convert edge paths into strings
        contigs = []
        for path in valid_paths:
            contigs.append(self.get_path(path))
        contigs = list(set(contigs))
        return contigs, valid_paths

    def map_strings(self, strings):
        Mapping = namedtuple('Mapping', ['s_st', 's_en', 'valid', 'epath'])
        index = self.index_edges()
        mappings = {}
        db_edges = list(self.nx_graph.edges(keys=True))

        # for db k == self.k
        k = max(index.keys())
        index = index[k]
        for s_id, string in strings.items():
            str_coords = []
            for i in range(len(string)-k+1):
                kmer = tuple(string[i:i+k])
                if kmer in index:
                    str_coords.append((index[kmer], i))

            if len(str_coords) == 0:
                mappings[s_id] = None
                continue

            (edge_ind, edge_coord), _ = str_coords[0]
            edge_path = [edge_ind]
            for (cur_edge_ind, cur_edge_coord), _ in str_coords[1:]:
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
            mappings[s_id] = Mapping(str_coords[0], str_coords[-1],
                                     valid_path, path)
        return mappings

    def write_dot(self, outfile, compact=False, export_pdf=True):
        if outfile[-3:] == 'dot':
            outfile = outfile[:-4]
        if not compact:
            graph = self.nx_graph
        else:
            graph = nx.MultiDiGraph()
            for s, e, key, data in self.nx_graph.edges(keys=True, data=True):
                graph.add_edge(s, e, key,
                               color=data[self.color],
                               label=data[self.label])
        dotfile = f'{outfile}.dot'
        nx.drawing.nx_pydot.write_dot(graph, dotfile)
        if export_pdf:
            pdffile = f'{outfile}.pdf'
            # https://stackoverflow.com/a/3516106
            cmd = ['dot','-Tpdf', dotfile,'-o', pdffile]
            proc = Popen(cmd,
                         shell=False,
                         stdin=None, stdout=None, stderr=None,
                         close_fds=True)

    def pickle_dump(self, fn):
        with open(fn, 'wb') as f:
            pickle.dump(self, f)
