# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
from collections import Counter
from string import ascii_uppercase, ascii_lowercase

import numpy as np
import pandas as pd

from utils.bio import read_bio_seqs, compress_homopolymer


class SD_Report:
    class SD_Record:
        def __init__(self, r_id, pandas_df,
                     max_gap, mean_monomer_len, gap_symb):
            self.r_id = r_id
            self.monomers = pandas_df.monomer.to_list()
            self.r_st = pandas_df.r_st.to_list()
            self.r_en = pandas_df.r_en.to_list()
            self.triples = list(zip(self.monomers, self.r_st, self.r_en))
            self.score = pandas_df.score.to_list()
            self.alt_call = pandas_df.alt_call.to_list()

            self.string = \
                [self.monomers[0][0] if self.alt_call[0] == 'None'
                 else gap_symb]
            self.gaps = []
            for tr1, tr2, rel2 in \
                    zip(self.triples[:-1],
                        self.triples[1:],
                        self.alt_call[1:]):
                gap_len = tr2[1] - tr1[2]
                if gap_len > max_gap:
                    self.gaps.append((tr1[2], tr2[1]))
                    self.string.append(gap_symb *
                                       int(round(gap_len / mean_monomer_len)))
                if rel2 == 'None':
                    self.string.append(tr2[0])
                else:
                    self.string.append(gap_symb)

            self.string = ''.join(self.string)
            perc_lower_case = np.mean([c.islower() for c in self.string
                                       if c.lower() != c.upper()])
            if perc_lower_case > 0.5:
                self.string = self.string[::-1].swapcase()
                self.strand = '-'
            else:
                self.strand = '+'

            self.string = self.string.strip(gap_symb)
            self.split_strings = self.string.split(gap_symb)
            self.split_strings = [s for s in self.split_strings if len(s)]

    def __init__(self, SD_report_fn, monomers_fn, max_gap=100, gap_symb='?'):
        self.gap_symb = gap_symb
        monomers = read_bio_seqs(monomers_fn)
        mean_monomer_len = \
            np.mean([len(monomer) for monomer in monomers.values()])

        self.monomer_names_map = {}
        for monomer_name, ucode, lcode in \
                zip(monomers.keys(), ascii_uppercase, ascii_lowercase):
            self.monomer_names_map[monomer_name] = ucode
            self.monomer_names_map[monomer_name + "'"] = lcode

        self.records = {}
        df = pd.read_csv(SD_report_fn, sep='\t',
                         header=None,
                         names=['r_id', 'monomer',
                                'r_st', 'r_en',
                                'score',
                                'alt_call',
                                'alt_score'
                                ])
        df.monomer = df.monomer.apply(lambda x: self.monomer_names_map[x])
        df = df.groupby('r_id')
        for r_id, group in df:
            self.records[r_id] = \
                self.SD_Record(r_id, group,
                               max_gap=max_gap,
                               mean_monomer_len=mean_monomer_len,
                               gap_symb=self.gap_symb)

    def get_monomer_strings(self):
        return {r_id: record.string for r_id, record in self.records.items()
                if len(record.string)}


def get_ngap_symbols(monostrings, compr_hmp=False, gap_symb='?'):
    cnt = 0
    for monostring in monostrings.values():
        monostring = monostring.tostring()
        if compr_hmp:
            monostring = compress_homopolymer(monostring)
        cnt += Counter(monostring)[gap_symb]
    return cnt


def get_stats(monostrings, verbose=False, return_stats=True):
    stats = {}
    monostrings_lens = [len(monostr) for r_id, monostr in monostrings.items()]
    stats['ntranslations'] = len(monostrings_lens)
    stats['min_len'] = np.min(monostrings_lens)
    stats['max_len'] = np.max(monostrings_lens)
    stats['mean_len'] = np.mean(monostrings_lens)
    stats['tot_len'] = np.sum(monostrings_lens)
    stats['ngaps'] = get_ngap_symbols(monostrings)
    stats['pgaps'] = stats['ngaps'] / stats['tot_len']
    stats['ngap_runs'] = get_ngap_symbols(monostrings, compr_hmp=True)

    if verbose:
        print(f'Number of translations: {stats["ntranslations"]}')
        print(f'Min length = {stats["min_len"]}')
        print(f'Mean length = {stats["mean_len"]}')
        print(f'Max length = {stats["max_len"]}')
        print(f'Total length = {stats["tot_len"]}')

        print(f'#(%) Gap symbols = {stats["ngaps"]} ({stats["pgaps"]})')
        print(f'#Gap runs = {stats["ngap_runs"]}')
    if return_stats:
        return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Report of SD", required=True)
    parser.add_argument("-m",
                        "--monomers",
                        help="Fasta with monomers",
                        required=True)
    params = parser.parse_args()
    sd_report = SD_Report(params.input, params.monomers)
    get_stats(monostrings=sd_report.get_monomer_strings(),
              verbose=True, return_stats=False)


if __name__ == "__main__":
    main()
