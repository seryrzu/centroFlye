import os
import logging
import unittest

from monomers.monomer_db import MonomerDB

logger = logging.getLogger("centroFlye.monomers_tests")

this_dirname = os.path.dirname(os.path.realpath(__file__))
monomers_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                           'test_dataset', 'ont', 'monomers.fasta')


class MonomersTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.monomer_db = MonomerDB.from_fasta_file(fn=monomers_fn)
        super(MonomersTests, self).__init__(*args, **kwargs)

    def test_monomer_db(self):
        self.assertEqual(self.monomer_db.get_size(), 12)
        mono_id = 'A_0_DXZ1*_doubled/1978_2147/R'
        self.assertEqual(
            self.monomer_db.id2index[mono_id], 0)
        self.assertEqual(mono_id, self.monomer_db.index2id[0])
        self.assertEqual(len(self.monomer_db.get_ids()),
                         self.monomer_db.get_size())
        self.monomer_db.get_seq_by_index(0)
        self.monomer_db.get_seq_by_id(mono_id)
        self.monomer_db.get_ids()

        for monomer_id in self.monomer_db.get_ids():
            index = self.monomer_db.id2index[monomer_id]
            monomer_id2 = self.monomer_db.index2id[index]
            self.assertEqual(monomer_id, monomer_id2)
