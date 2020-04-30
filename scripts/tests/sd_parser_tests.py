import os
import unittest

from sd_parser.sd_parser import SD_Report


this_dirname = os.path.dirname(os.path.realpath(__file__))
monomers_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                           'test_dataset', 'assembly', 'monomers.fasta')
sd_report_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                            'test_dataset', 'assembly',
                            'assembly_final_decomposition.tsv')
sequences_fn = os.path.join(this_dirname, os.path.pardir, os.path.pardir,
                            'test_dataset', 'assembly',
                            'T2TX7_hpc.fasta')

sd_report_wo_hpc_fn = \
    os.path.join(this_dirname, os.path.pardir, os.path.pardir,
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
        monostring = list(self.sd_report.monostring_set.monostrings.values())[0]
        monoinstances = monostring.monoinstances
        self.assertEqual(monoinstances[0].index, 9)
        self.assertEqual(monoinstances[-1].index, 10)
        self.assertEqual(monoinstances[0].st, 125)
        self.assertEqual(monoinstances[0].en, 239)
        self.assertEqual(len(monoinstances[0].nucl_segment), 239 - 125)
        self.assertFalse(monostring.is_reversed)
        self.assertEqual(monostring.pref_cut, 1)
        self.assertEqual(monostring.suf_cut, 4)


class SDParserWOHPCTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_wo_hpc_fn,
                                   monomers_fn=monomers_for_report_wo_hpc_fn,
                                   sequences_fn=sequences_for_report_wo_hpc_fn,
                                   hpc=False)
        super(SDParserWOHPCTests, self).__init__(*args, **kwargs)

    def test_sd_parser(self):
        monostring = list(self.sd_report.monostring_set.monostrings.values())[0]
        nucl_sequence = monostring.nucl_sequence
        monoinstances = monostring.monoinstances
        self.assertEqual(monostring.pref_cut, 1)
        self.assertEqual(monostring.suf_cut, 3)
        self.assertEqual(monoinstances[0].st, 464, msg=None)
        self.assertEqual(monoinstances[0].en, 635, msg=None)
        self.assertEqual(monoinstances[0].index, 884)
        self.assertEqual(monostring.string[0], 884)
        st = monostring.monoinstances[0].st
        en = monostring.monoinstances[0].en
        segment = monoinstances[0].nucl_segment
        self.assertEqual(nucl_sequence[st:en], segment, msg=None)
        self.assertTrue(monostring.is_reversed)
