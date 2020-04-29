import os
import logging
import unittest

from monomers.monomer_db import MonomerDB
from sd_parser.sd_parser import SD_Report
from utils.bio import RC


this_dirname = os.path.dirname(os.path.realpath(__file__))
monomers_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                           'test_dataset', 'assembly', 'monomers.fasta')
sd_report_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                           'test_dataset', 'assembly',
                            'assembly_final_decomposition.tsv')
sequences_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                            'test_dataset', 'assembly',
                            'T2TX7_hpc.fasta')

sd_report_wo_hpc_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                                   'test_dataset', 'hifi_hpc',
                                   'final_decomposition_without_hpc.tsv')
monomers_for_report_wo_hpc_fn = \
    os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                 'test_dataset', 'hifi_hpc',
                 'monomers_for_report_without_hpc.fasta')

sequences_for_report_wo_hpc_fn = \
    os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                 'test_dataset', 'hifi_hpc',
                 'centromeric_reads_for_report_wo_hpc.fasta')

class SDParserTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_fn,
                                   monomers_fn=monomers_fn,
                                   sequences_fn=sequences_fn)
        super(SDParserTests, self).__init__(*args, **kwargs)

    def test_sd_parser(self):
        monostring = list(self.sd_report.monostrings.values())[0]
        self.assertEqual(monostring.index2monoinstance[0].index, 9)
        max_index = max(monostring.index2monoinstance)
        self.assertEqual(monostring.index2monoinstance[max_index].index, 1)
        self.assertEqual(monostring.index2monoinstance[0].st, 125)
        self.assertEqual(monostring.index2monoinstance[0].en, 239)
        self.assertEqual(len(monostring.index2monoinstance[0].segment), 239-125)


class SDParserWOHPCTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_wo_hpc_fn,
                                   monomers_fn=monomers_for_report_wo_hpc_fn,
                                   sequences_fn=sequences_for_report_wo_hpc_fn,
                                   hpc=False)
        super(SDParserWOHPCTests, self).__init__(*args, **kwargs)

    def test_sd_parser(self):
        monostring = list(self.sd_report.monostrings.values())[0]
        nucl_sequence = monostring.nucl_sequence
        self.assertEqual(monostring.index2monoinstance[0].index, 541)
        self.assertEqual(monostring.string[0], 1506)
        st = monostring.index2monoinstance[0].st
        en = monostring.index2monoinstance[0].en
        segment = monostring.index2monoinstance[0].segment
        self.assertEqual(nucl_sequence[st:en], RC(segment), msg=None)
