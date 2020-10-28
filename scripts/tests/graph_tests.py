import os
import logging
import unittest

import networkx as nx
from sequence_graph.seq_graph import SequenceGraph
from sequence_graph.db_graph import DeBruijnGraph

logger = logging.getLogger("centroFlye.monomers_tests")

this_dirname = os.path.dirname(os.path.realpath(__file__))


class DBGraphTests(unittest.TestCase):
    def test_collapsing(self):
        kmers = ['ACGT','ACGA', 'ACGC', 'CGCC', 'GCCC', 'CCCC']
        db_graph = DeBruijnGraph.from_kmers(kmers, collapse=False)
        self.assertDictEqual(db_graph.nodeindex2label,
                             {0: tuple('ACG'), 1: tuple('CGT'),
                              2: tuple('CGA'), 3: tuple('CGC'),
                              4: tuple('GCC'), 5: tuple('CCC')})

        self.assertDictEqual(db_graph.nodelabel2index,
                             {tuple('ACG'): 0, tuple('CGT'): 1,
                              tuple('CGA'): 2, tuple('CGC'): 3,
                              tuple('GCC'): 4, tuple('CCC'): 5})


        self.assertListEqual(list(db_graph.nx_graph.edges),
                             [(0, 1, 0), (0, 2, 0), (0, 3, 0), (3, 4, 0),
                              (4, 5, 0), (5, 5, 0)])

        self.assertTupleEqual(tuple('ACGCCC'),
                              db_graph.get_path([(0, 3, 0),
                                                 (3, 4, 0),
                                                 (4, 5, 0)]))
        db_graph.collapse_nonbranching_paths()

        self.assertListEqual(list(db_graph.nx_graph.edges),
                             [(0, 1, 0), (0, 2, 0), (0, 5, 0), (5, 5, 0)])

        self.assertTupleEqual(tuple('ACGCCC'),
                              db_graph.get_path([(0, 5, 0)]))

        self.assertTupleEqual(tuple('C'),
                              db_graph.get_path([(5, 5, 0)]))
        edges = list(db_graph.nx_graph.edges(keys=True))
        self.assertListEqual(edges,
                             [(0, 1, 0), (0, 2, 0), (0, 5, 0), (5, 5, 0)])
        k = len(kmers[0])
        index = db_graph.index_edges()
        self.assertDictEqual(index[k],
                             {tuple('ACGT'): (0, 0), tuple('ACGA'): (1, 0),
                              tuple('ACGC'): (2, 0), tuple('CGCC'): (2, 1),
                              tuple('GCCC'): (2, 2), tuple('CCCC'): (3, 0)})
        mappings = db_graph.map_strings({0: 'ACGCCCCCC'}, neutral_symbs=set())

    def test_contigs(self):
        kmers = ['TAA', 'GAC', 'TCA', 'AAC', 'ACA', 'CAA']
        db_graph = DeBruijnGraph.from_kmers(kmers)
        contigs, valid_paths = db_graph.get_contigs()
        self.assertListEqual(sorted(contigs),
                             sorted([tuple('TAACAA'),
                                     tuple('GACAAC'),
                                     tuple('TCAACA')]))
