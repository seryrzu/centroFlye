import logging
import unittest

from monomers.monomers import MonomerDB

logger = logging.getLogger("centroFlye.monomers_tests")

MONOMERS_FN = '../../test_dataset/ont/monomers.fasta'


class MonomersTests(unittest.TestCase):
    def test_monomer_db(self):
        monomer_db = MonomerDB.from_fasta_file(fn=MONOMERS_FN)
        self.assertEqual(monomer_db.get_total_monomers(), 12)
        self.assertEqual(
            monomer_db.names2ident['A_0_DXZ1*_doubled/1978_2147/R'], 0)
        self.assertEqual(
            'A_0_DXZ1*_doubled/1978_2147/R', monomer_db.ident2names[0])
        self.assertEqual(
            len(monomer_db.get_names()), monomer_db.get_total_monomers())
