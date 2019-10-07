import argparse
from string import ascii_uppercase, ascii_lowercase
import numpy as np
import pandas as pd

from utils.bio import read_bio_seqs


class LRD_Report:
    class LRD_Record:
        def __init__(self, r_id, pandas_df, max_gap, mean_monomer_len):
            self.r_id = r_id
            self.monomers = pandas_df.monomer.to_list()
            self.r_st = pandas_df.r_st.to_list()
            self.r_en = pandas_df.r_en.to_list()
            self.triples = list(zip(self.monomers, self.r_st, self.r_en))
            self.score = pandas_df.score.to_list()
            self.reliability = pandas_df.reliability.to_list()

            self.string = [self.monomers[0][0] if self.reliability[0] == '+' else '=']
            self.gaps = []
            for tr1, tr2, rel2 in \
                    zip(self.triples[:-1],
                        self.triples[1:],
                        self.reliability[1:]):
                gap_len = tr2[1] - tr1[2]
                if gap_len > max_gap:
                    self.gaps.append((tr1[2], tr2[1]))
                    self.string.append('=' * int(round(gap_len / mean_monomer_len)))
                if rel2 == '+':
                    self.string.append(tr2[0])
                else:
                    self.string.append('=')

            self.string = ''.join(self.string)
            perc_lower_case = np.mean([c.islower() for c in self.string \
                                       if c.lower() != c.upper()])
            if perc_lower_case > 0.5:
                self.string = self.string[::-1].swapcase()
                self.strand = '-'
            else:
                self.strand = '+'

            self.string = self.string.strip('=')
            self.split_strings = self.string.strip('=').split('=')
            self.split_strings = [s for s in self.split_strings if len(s)]


    def __init__(self, lrd_report_fn, monomers_fn, max_gap=100):
        monomers = read_bio_seqs(monomers_fn)
        mean_monomer_len = \
            np.mean([len(monomer) for monomer in monomers.values()])

        self.monomer_names_map = {}
        for monomer_name, ucode, lcode in \
                zip(monomers.keys(), ascii_uppercase, ascii_lowercase):
            self.monomer_names_map[monomer_name] = ucode
            self.monomer_names_map[monomer_name + "'"] = lcode

        self.records = {}
        df = pd.read_csv(lrd_report_fn, sep='\t',
                         header=None,
                         names=['r_id', 'monomer',
                                'r_st', 'r_en',
                                'score',
                                'reliability'])
        df.monomer = df.monomer.apply(lambda x: self.monomer_names_map[x])
        df = df.groupby('r_id')
        for r_id, group in df:
            self.records[r_id] = \
                self.LRD_Record(r_id, group,
                                max_gap=max_gap,
                                mean_monomer_len=mean_monomer_len)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Report of LRD", required=True)
    parser.add_argument("-m",
                        "--monomers",
                        help="Fasta with monomers",
                        required=True)
    params = parser.parse_args()
    lrd_report = LRD_Report(params.input, params.monomers)


if __name__ == "__main__":
    main()
