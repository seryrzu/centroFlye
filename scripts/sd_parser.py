# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
from collections import Counter
# from string import ascii_uppercase, ascii_lowercase
import unicodedata as ud

from scipy.stats import binom_test
import numpy as np
import pandas as pd

from utils.bio import read_bio_seqs, compress_homopolymer


class MonoString:
    # string is stored as a list
    def __init__(self, name, rev_monomer_index,
                 string=list(), mono2nucl=dict(),
                 gap_symb='?', strand='+'):
        self.name = name
        self.rev_monomer_index = rev_monomer_index
        self.string = string.copy()
        self.mono2nucl = mono2nucl.copy()
        self.gap_symb = gap_symb
        self.strand = strand

    @classmethod
    def FromSDRecord(cls, name, monomers, starts, ends, reliability,
                     max_gap, mean_monomer_len, gap_symb, rev_monomer_index):

        monostring = cls(name=name, rev_monomer_index=rev_monomer_index)

        # Ignore the first and last monomer, as they may be incomplete
        monomers, starts, ends = monomers[1:-1], starts[1:-1], ends[1:-1]
        
        if reliability[0] == '+':
            monostring.add_monomer(monomer=monomers[0],
                                   st=starts[0],
                                   en=ends[0])
        else:
            monostring.add_gap(1)

        triple1 = zip(monomers[:-1], starts[:-1], ends[:-1])
        triple2 = zip(monomers[1:], starts[1:], ends[1:])
        for (m1, s1, e1), (m2, s2, e2), rel2 in \
                zip(triple1, triple2, reliability[1:]):
            gap_len = s2 - e1
            if gap_len > max_gap:
                int_gap_len = int(round(gap_len / mean_monomer_len))
                monostring.add_gap(int_gap_len)
            if rel2 == '+':
                monostring.add_monomer(m2, s2, e2)
            else:
                monostring.add_gap(1)

        monostring.assert_validity()
        monostring.strip()
        monostring.check_reverse()
        return monostring

    # def tostring(self):
    #     return ''.join(self.string)

    def __len__(self):
        return self.string.__len__()

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.string[sub.start:sub.stop:sub.step]
            # sublist = ''.join(sublist)
            return sublist
        return self.string[sub]

    def __setitem__(self, sub, value):
        if isinstance(sub, slice):
            if isinstance(value, str):
                value = list(value)
            self.string[sub.start:sub.stop:sub.step] = value
        else:
            self.string[sub] = value


    def assert_validity(self):
        for coord, (c, _, _) in self.mono2nucl.items():
            if coord < 0:
                print(coord, c)
            if coord >= len(self.string):
                print(coord, c, len(self.stirng))
            assert c == self.string[coord]

    def add_monomer(self, monomer, st, en):
        self.mono2nucl[len(self.string)] = \
            (monomer, st, en)
        self.string.append(monomer)

    def add_gap(self, length):
        self.string += [self.gap_symb] * length

    def check_reverse(self):
        def swapcase_monomer(m):
            if not isinstance(m, int):
                return m
            if m >= self.rev_monomer_index:
                return m - self.rev_monomer_index
            return m + self.rev_monomer_index

        perc_lower_case = np.mean([c > self.rev_monomer_index for c in self.string
                                   if isinstance(c, int)])
        if perc_lower_case > 0.5:
            self.string = [swapcase_monomer(m) for m in self.string[::-1]]
            self.strand = '-'
            rev_mono2nucl = {}
            for coord, (monomer, st, en) in self.mono2nucl.items():
                rev_coord = len(self.string) - coord - 1
                rev_mono2nucl[rev_coord] = (swapcase_monomer(monomer), en, st)
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

    def split(self, c, min_length):
        resulting_monostrings = {}

        split_strings = []
        i = 0
        while i < len(self.string):
            j = i
            while j < len(self.string) and self.string[j] != c:
                j += 1
            if j > i:
                split_strings.append(self.string[i:j])
            i = j + 1

        pos = list(self.mono2nucl)
        coord_pos = 0
        cumm_len = 0

        for i, split_string in enumerate(split_strings):
            if len(split_string) < min_length:
                cumm_len += len(split_string) + 1
                continue

            split_mono2nucl = {}
            while coord_pos < len(pos) and pos[coord_pos] < cumm_len:
                coord_pos += 1

            while pos[coord_pos] < + len(split_string):
                p = pos[coord_pos]
                assert split_string[p - cumm_len] == self.mono2nucl[p][0]
                split_mono2nucl[p - cumm_len] = self.mono2nucl[p]
                coord_pos += 1

            new_monostring = MonoString(name=self.name,
                                        rev_monomer_index=self.rev_monomer_index,
                                        string=list(split_string),
                                        mono2nucl=split_mono2nucl,
                                        gap_symb=self.gap_symb,
                                        strand=self.strand)
            new_monostring.assert_validity()

            resulting_monostrings[(self.name, i)] = new_monostring
            cumm_len += len(split_string) + 1
        return resulting_monostrings


def ident_hybrids(df, monomer_names_map,
                  max_ident_diff=5, min_occ=500):
    def count_pairs(df):
        pairs = []
        for m1, m2, rel in df[['monomer',
                               'sec_monomer',
                               'reliability']].to_numpy():
            if rel == '+':
                m1, m2 = m1.upper(), m2.upper()
                m1, m2 = min(m1, m2), max(m1, m2)
                pairs.append((m1, m2))
        cnt = Counter(pairs)
        return cnt

    def get_med_pair_p(all_cnt, pairing_cnt, min_occ=min_occ):
        ps = [pairing_cnt[pair] / all_cnt[pair]
              for pair in all_cnt.keys()
              if all_cnt[pair] >= min_occ]
        return np.median(ps)

    def recr_sign_pairs(all_cnt, pairing_cnt, med_pair_p,
                        alpha=0.05, min_occ=min_occ, min_pc=5):
        pairs = []
        for pair, ac in all_cnt.items():
            pc = pairing_cnt[pair]
            if ac < min_occ or pc < min_pc:
                continue
            pv = binom_test(pc, ac, p=med_pair_p, alternative='greater')

            if pv * len(all_cnt) < alpha:
                print(pv, ac, pc, pair)
                pairs.append((pv, pair))
        pairs.sort()

        max_char = max(monomer_names_map.values())
        max_index = ascii_lowercase.index(max_char)
        left_letters = len(ascii_lowercase) - max_index

        print(pairs)
        pairs = [pair for (pv, pair) in pairs[:left_letters]]
        return pairs

    def get_pair2char(pairs):
        pair2char = {}
        lower_alphabet = get_unicode_lowercase()
        upper_alphabet = get_unicode_uppercase()
        max_char = max(monomer_names_map.values())
        index = lower_alphabet.index(max_char) + 1
        for pair in pairs:
            pair2char[pair] = upper_alphabet[index]
            pair2char[pair[::-1]] = upper_alphabet[index]

            lower_pair = (pair[0].lower(), pair[1].lower())
            pair2char[lower_pair] = lower_alphabet[index]
            pair2char[lower_pair[::-1]] = lower_alphabet[index]
            index += 1
        return pair2char

    def transform_df(df_pairing, pair2char):
        for i, row in df_pairing.iterrows():
            pair = (row.monomer, row.sec_monomer)
            if pair in pair2char:
                df.at[i, 'monomer'] = pair2char[pair]
        return df

    all_cnt = count_pairs(df)
    close_identity = (abs(df.identity - df.sec_identity) < max_ident_diff) & (df.reliability == '+')
    pairing = close_identity
    df_pairing = df[pairing]
    pairing_cnt = count_pairs(df_pairing)
    med_pair_p = get_med_pair_p(all_cnt, pairing_cnt)
    pairs = recr_sign_pairs(all_cnt, pairing_cnt, med_pair_p)
    pair2char = get_pair2char(pairs)

    transformed_df = transform_df(df[pairing], pair2char)
    return transformed_df, pair2char


def get_unicode():
    return ''.join(chr(i) for i in range(65536))


def get_unicode_uppercase():
    all_unicode = get_unicode()
    uppercase = ''.join(c for c in all_unicode
                        if ud.category(c)=='Lu')
    uppercase = uppercase.lower().upper()
    uppercase = Counter(uppercase)
    uppercase = ''.join(c for c in uppercase if uppercase[c] == 1)
    return uppercase


def get_unicode_lowercase():
    unicode_uppercase = get_unicode_uppercase()
    return unicode_uppercase.lower()


class SD_Report:
    def __init__(self, SD_report_fn, monomers_fn,
                 max_gap=100,
                 gap_symb='?',
                 hpc=True):
        self.gap_symb = gap_symb
        self.monomers = read_bio_seqs(monomers_fn)
        mean_monomer_len = \
            np.mean([len(monomer) for monomer in self.monomers.values()])

        self.monomer_names_map = {}
        self.rev_monomer_names_map = {}
        for i, monomer_name in enumerate(self.monomers.keys()):
            self.monomer_names_map[monomer_name] = i 
            self.monomer_names_map[monomer_name + "'"] = i + len(self.monomers)
            self.rev_monomer_names_map[i] = monomer_name
        self.monomer_names_map['None'] = '*' # '?' is reserved for '*'

        self.monostrings = {}
        if hpc:
            names = ['r_id', 'monomer',
                     'r_st', 'r_en',
                     'identity',
                     'sec_monomer', 'sec_identity',
                     'homo_monomer', 'homo_identity',
                     'sec_homo_monomer', 'sec_homo_identity',
                     'reliability'
                    ]
        else:
            names = ['r_id', 'monomer',
                     'r_st', 'r_en',
                     'identity',
                     'sec_monomer', 'sec_identity',
                     'reliability'
                    ]
            
        df = pd.read_csv(SD_report_fn, sep='\t',
                         header=None,
                         names=names)
        df.monomer = df.monomer.apply(lambda x: self.monomer_names_map[x])
        df.sec_monomer = df.sec_monomer.apply(lambda x: self.monomer_names_map[x])

        for r_id, group in df.groupby('r_id'):
            group = group.sort_values(by=['r_st'])
            self.monostrings[r_id] = \
                MonoString.FromSDRecord(name=r_id,
                                        monomers=group.monomer.to_list(),
                                        starts=group.r_st.to_list(),
                                        ends=group.r_en.to_list(),
                                        reliability=group.reliability.to_list(),
                                        max_gap=max_gap,
                                        mean_monomer_len=mean_monomer_len,
                                        gap_symb=self.gap_symb,
                                        rev_monomer_index=len(self.monomers))
        self.df = df


def get_ngap_symbols(monostrings, compr_hmp=False, gap_symb='?'):
    cnt = 0
    for monostring in monostrings.values():
        if compr_hmp:
            monostring = compress_homopolymer(monostring, return_list=True)
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
    get_stats(monostrings=sd_report.monostrings,
              verbose=True, return_stats=False)


if __name__ == "__main__":
    main()
