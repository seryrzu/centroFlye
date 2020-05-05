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

    def test_sd_parser_general(self):
        monostring = list(
            self.sd_report.monostring_set.monostrings.values())[0]
        monoinstances = monostring.monoinstances
        self.assertEqual(monoinstances[0].get_monoindex(), 9)
        self.assertEqual(monoinstances[-3].get_monoindex(), 10)
        self.assertEqual(monoinstances[0].st, 125)
        self.assertEqual(monoinstances[0].en, 239)
        self.assertEqual(len(monoinstances[0].nucl_segment), 239 - 125)
        self.assertFalse(monostring.is_reversed)

    def test_monomerinstances_classification(self):
        monostring = list(
            self.sd_report.monostring_set.monostrings.values())[0]
        monomerinstances_dict = monostring.classify_monomerinstances()
        self.assertEqual(len(monomerinstances_dict[0]), 1505)
        self.assertEqual(len(monomerinstances_dict[9]), 1517)
        monomerinstances_0 = monostring.get_monomerinstances_by_monoindex(0)
        self.assertEqual(len(monomerinstances_0), 1505)
        monomerinstances_0_withUnrel = \
            monostring.get_monomerinstances_by_monoindex(mono_index=0,
                                                         only_reliable=False)
        self.assertEqual(len(monomerinstances_0_withUnrel), 1508)

        monostring_set = self.sd_report.monostring_set
        monomerinstances_dict = monostring_set.classify_monomerinstances()
        self.assertEqual(len(monomerinstances_dict[0]), 1505)
        self.assertEqual(len(monomerinstances_dict[9]), 1517)
        monomerinstances_0_withUnrel = \
            monostring_set.get_monomerinstances_by_monoindex(mono_index=0,
                                                             only_reliable=False)
        self.assertEqual(len(monomerinstances_0_withUnrel), 1508)


class SDParserWOHPCTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_wo_hpc_fn,
                                   monomers_fn=monomers_for_report_wo_hpc_fn,
                                   sequences_fn=sequences_for_report_wo_hpc_fn,
                                   hpc=False)
        super(SDParserWOHPCTests, self).__init__(*args, **kwargs)

    def test_sd_parser(self):
        monostring = list(
            self.sd_report.monostring_set.monostrings.values())[0]
        nucl_sequence = monostring.nucl_sequence
        monoinstances = monostring.monoinstances
        mi = monoinstances[2]
        self.assertEqual(mi.st, 464, msg=None)
        self.assertEqual(mi.en, 635, msg=None)
        self.assertEqual(mi.get_monoindex(), 884)
        self.assertEqual(monostring.raw_monostring[2], 884)
        self.assertEqual(nucl_sequence[mi.st:mi.en], mi.nucl_segment, msg=None)
        self.assertTrue(monostring.is_reversed)
