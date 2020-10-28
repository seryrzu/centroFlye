# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging
import os

import networkx as nx
import numpy as np

from sequence_graph.seq_graph import SequenceGraph
from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.sequence_graph.db_graph_comparison")


class DeBruijnGraphComparison(SequenceGraph):
    cov_asm1 = 'cov_asm1'
    cov_asm2 = 'cov_asm2'
    col_asm1 = 'blue'
    col_asm2 = 'red'
    col_both_samecov = 'green'
    col_both_diffcov = 'orange'

    def __init__(self, nx_graph, nodeindex2label, nodelabel2index, k,
                 name_asm1, name_asm2,
                 collapse=True):
        self.k = k  # length of an edge in the uncompressed graph
        self.name_asm1 = name_asm1
        self.name_asm2 = name_asm2
        super().__init__(nx_graph=nx_graph,
                         nodeindex2label=nodeindex2label,
                         nodelabel2index=nodelabel2index,
                         collapse=collapse)

    @classmethod
    def _generate_label(cls, par_dict):
        length = par_dict[cls.length]
        cov_asm1 = par_dict[cls.cov_asm1]
        cov_asm2 = par_dict[cls.cov_asm2]
        name_asm1 = par_dict['name_asm1']
        name_asm2 = par_dict['name_asm2']
        if cov_asm1 is None:
            mean_cov_asm2 = int(np.mean(cov_asm2))
            label = f'len={length}\n{name_asm2}_Cov={mean_cov_asm2}'
        elif cov_asm2 is None:
            mean_cov_asm1 = int(np.mean(cov_asm1))
            label = f'len={length}\n{name_asm1}_Cov={mean_cov_asm1}'
        else:
            assert cov_asm1 is not None and cov_asm2 is not None
            mean_cov_asm1 = int(np.mean(cov_asm1))
            mean_cov_asm2 = int(np.mean(cov_asm2))
            label = f'len={length}\n{name_asm1}_Cov={mean_cov_asm1}\n' + \
                    f'{name_asm2}_Cov={mean_cov_asm2}'
        return label

    @classmethod
    def from_monoassemblies(cls, monoasm1, monoasm2, k, collapse=True,
                            outdir=None):
        def add_kmer(kmer, cov_asm1, cov_asm2, color):
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
            if cov_asm1 is not None:
                cov_asm1 = [cov_asm1]
            if cov_asm2 is not None:
                cov_asm2 = [cov_asm2]
            label = \
                cls._generate_label({cls.length: length,
                                     cls.cov_asm1: cov_asm1,
                                     cls.cov_asm2: cov_asm2,
                                     'name_asm1': monoasm1.seq_id,
                                     'name_asm2': monoasm2.seq_id})
            nx_graph.add_edge(prefix_node_ind, suffix_node_ind,
                              string=kmer,
                              length=length,
                              cov_asm1=cov_asm1,
                              cov_asm2=cov_asm2,
                              label=label,
                              color=color)

        kmers_cnt1 = monoasm1.get_kmer_index(maxk=k, mink=k, positions=False)[k]
        kmers_cnt2 = monoasm2.get_kmer_index(maxk=k, mink=k, positions=False)[k]

        kmers1 = set(kmers_cnt1.keys())
        kmers2 = set(kmers_cnt2.keys())
        kmers_both = kmers1 & kmers2
        kmers1 = kmers1 - kmers_both
        kmers2 = kmers2 - kmers_both

        nx_graph = nx.MultiDiGraph()
        nodeindex2label = {}
        nodelabel2index = {}

        for kmer in kmers_both:
            cov_asm1 = kmers_cnt1[kmer]
            cov_asm2 = kmers_cnt2[kmer]
            color = cls.col_both_samecov if cov_asm1 == cov_asm2 else cls.col_both_diffcov
            add_kmer(kmer=kmer,
                     cov_asm1=cov_asm1,
                     cov_asm2=cov_asm2,
                     color=color)
        for kmer in kmers1:
            add_kmer(kmer=kmer,
                     cov_asm1=kmers_cnt1[kmer],
                     cov_asm2=None,
                     color=cls.col_asm1)
        for kmer in kmers2:
            add_kmer(kmer=kmer,
                     cov_asm1=None,
                     cov_asm2=kmers_cnt2[kmer],
                     color=cls.col_asm2)

        color3graph = cls(nx_graph=nx_graph,
                          nodeindex2label=nodeindex2label,
                          nodelabel2index=nodelabel2index,
                          k=k,
                          collapse=collapse,
                          name_asm1=monoasm1.seq_id,
                          name_asm2=monoasm2.seq_id)
        if outdir is not None:
            smart_makedirs(outdir)
            dot_file = os.path.join(outdir, f'c3g_k{k}.dot')
            color3graph.write_dot(outfile=dot_file,
                                  export_pdf=True,
                                  compact=True)
            c3g_pickle_file = os.path.join(outdir, f'c3g_k{k}.pickle')
            color3graph.pickle_dump(c3g_pickle_file)

        return color3graph

    def _add_edge(self, node, color, string,
                  in_node, out_node,
                  in_data, out_data,
                  edge_len):
        in_cov_asm1 = in_data[self.cov_asm1]
        out_cov_asm1 = out_data[self.cov_asm1]

        in_cov_asm2 = in_data[self.cov_asm2]
        out_cov_asm2 = out_data[self.cov_asm2]

        assert (in_cov_asm1 is None) == (out_cov_asm1 is None)
        assert (in_cov_asm2 is None) == (out_cov_asm2 is None)
        assert not (in_cov_asm1 is None and out_cov_asm1 is None and
                    in_cov_asm2 is None and out_cov_asm2 is None)
        if in_cov_asm1 is not None:
            cov_asm1 = sorted(in_cov_asm1 + out_cov_asm1)
            assert len(cov_asm1) == edge_len
        else:
            cov_asm1 = None

        if in_cov_asm2 is not None:
            cov_asm2 = sorted(in_cov_asm2 + out_cov_asm2)
            assert len(cov_asm2) == edge_len
        else:
            cov_asm2 = None

        label = \
            self._generate_label({self.length: edge_len,
                                  self.cov_asm1: cov_asm1,
                                  self.cov_asm2: cov_asm2,
                                  'name_asm1': self.name_asm1,
                                  'name_asm2': self.name_asm2})
        self.nx_graph.add_edge(in_node, out_node,
                               string=string,
                               length=edge_len,
                               cov_asm1=cov_asm1,
                               cov_asm2=cov_asm2,
                               label=label,
                               color=color)
