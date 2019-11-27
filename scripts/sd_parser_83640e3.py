import argparse
from collections import Counter
from string import ascii_uppercase, ascii_lowercase

import numpy as np
import pandas as pd

from utils.bio import read_bio_seqs, compress_homopolymer


class MonoString:
    def __init__(self, gap_symb):
        self.string = []
        self.mono2nucl = {}
        self.gap_symb = gap_symb
        self.trimmed_len = 0
        self.strand = '+'

    def add_gap(self, length):
        self.string += [self.gap_symb] * length

    def add_monomer(self, monomer, st, en):
        self.mono2nucl[len(self.string)] = (monomer, st, en)
        self.string.append(monomer)

    def assert_validity(self):
        for coord, (c, _, _) in self.mono2nucl.items():
            # print(coord, c, self.string[coord])
            assert c == self.string[coord]

    def check_reverse(self):
        perc_lower_case = np.mean([c.islower() for c in self.string
                                   if c.lower() != c.upper()])
        if perc_lower_case > 0.5:
            self.string = self.string[::-1]
            self.string = [m.swapcase() for m in self.string]
            self.strand = '-'
            rev_mono2nucl = {}
            for coord, (monomer, st, en) in self.mono2nucl.items():
                rev_coord = len(self.string) - coord - 1
                rev_mono2nucl[rev_coord] = (monomer.swapcase(), en, st)
            self.mono2nucl = rev_mono2nucl
        self.assert_validity()

    def trim_read(self, left, right):
        self.string = self.string[left:right]
        self.mono2nucl = \
            {k - left: v for k, v in self.mono2nucl.items()
             if left <= k < right}
        self.assert_validity()

    def strip(self):
        i, j = 0, len(self.string) - 1
        while i < len(self.string) and self.string[i] == self.gap_symb:
            i += 1
        while j >= 0 and self.string[j] == self.gap_symb:
            j -= 1
        self.trim_read(i, j+1)

    def tostring(self):
        return ''.join(self.string)

    def __len__(self):
        return self.string.__len__()

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.string[sub.start:sub.stop:sub.step]
            sublist = ''.join(sublist)
            return sublist
        return self.string[sub]


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

            self.string = MonoString(gap_symb=gap_symb)

            if self.alt_call[0] == '+':
                monomer, st, en = self.triples[0]
                self.string.add_monomer(monomer=monomer, st=st, en=en)
            else:
                self.string.add_gap(1)

            self.gaps = []
            for tr1, tr2, rel2 in zip(self.triples[:-1],
                                      self.triples[1:],
                                      self.alt_call[1:]):
                gap_len = tr2[1] - tr1[2]
                if gap_len > max_gap:
                    self.gaps.append((tr1[2], tr2[1]))
                    int_gap_len = int(round(gap_len / mean_monomer_len))
                    self.string.add_gap(int_gap_len)
                if rel2 == '+':
                    self.string.add_monomer(*tr2)
                else:
                    self.string.add_gap(1)

            self.string.strip()
            self.string.check_reverse()

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
                                'alt_call'
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
        return {r_id: record.string
                for r_id, record in self.records.items() if len(record.string)}


def get_ngap_symbols(monostrings, compr_hmp=False, gap_symb='?'):
    cnt = 0
    for monostring in monostrings.values():
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
