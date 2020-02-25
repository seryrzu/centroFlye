import argparse
from collections import Counter
from string import ascii_lowercase

from joblib import Parallel, delayed
import numpy as np
from scipy.stats import binom_test
from statsmodels.stats.proportion import proportions_ztest

from sd_parser import SD_Report
from utils.bio import read_bio_seqs, hybrid_alignment, RC, calc_identity, write_bio_seqs


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd", required=True, help="SD report")
    parser.add_argument("--monomers", required=True, help="Monomers")
    parser.add_argument("--reads", required=True, help="Cen reads")
    parser.add_argument("--max-ident-diff", type=int, default=5)
    parser.add_argument("--outfile", required=True)
    params = parser.parse_args()
    return params


def read_input(params):
    sd_report = SD_Report(SD_report_fn=params.sd,
                          monomers_fn=params.monomers,
                          ident_hybrid=False)
    reads = read_bio_seqs(params.reads)
    return sd_report, reads


def ident_putative_pairs(df, monomer_names_map,
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
                pairs.append((pv, pair))
        pairs.sort()

        max_char = max(monomer_names_map.values())
        max_index = ascii_lowercase.index(max_char)
        left_letters = len(ascii_lowercase) - max_index

        pairs = [pair for (pv, pair) in pairs[:left_letters]]
        return pairs

    all_cnt = count_pairs(df)
    pairing = \
        (abs(df.identity - df.sec_identity) < max_ident_diff) & \
        (df.reliability == '+')
    df_pairing = df[pairing]
    pairing_cnt = count_pairs(df_pairing)
    med_pair_p = get_med_pair_p(all_cnt, pairing_cnt)
    pairs = recr_sign_pairs(all_cnt, pairing_cnt, med_pair_p)
    return pairs


def extract_hybrid(pair, df, reads, max_ident_diff, A, B, sign_lev=0.05):
    def extract_read_segms(df, reads):
        segms = []
        for i, row in df.iterrows():
            r_id, st, en = row.r_id, row.r_st, row.r_en + 1
            monomer = row.monomer
            segm = reads[r_id][st:en]
            if monomer.islower():
                segm = RC(segm)
            segms.append(segm)
        return segms

    def get_cuts(segms, A, B, t=50):
        results = \
            Parallel(n_jobs=t)(delayed(hybrid_alignment)(segm, A, B)
                               for segm in segms)
        jumps = []
        for _, _, max_sc, sec_max_sc, jump, orient in results:
            if max_sc * 0.9 > sec_max_sc:
                jumps.append((jump[1:], orient))
        return jumps

    def n_improved_pairs(segms, A, B, hybrid):
        n_improved = 0
        for segm in segms:
            identA = calc_identity(segm, A)
            identB = calc_identity(segm, B)
            ident_hybrid = calc_identity(segm, hybrid)
            if ident_hybrid > max(identA, identB):
                n_improved += 1
        return n_improved

    bases = list(pair)
    bases += [b.lower() for b in bases]
    AB = np.where(df.monomer.isin(bases) &
                  df.sec_monomer.isin(bases) &
                  (abs(df.identity - df.sec_identity) < max_ident_diff) &
                  (df.sec_identity > 85) &
                  df.reliability.isin(['+']))[0]
    dfAB = df.iloc[AB]

    read_segms = extract_read_segms(dfAB, reads)
    np.random.shuffle(read_segms)
    if len(read_segms) < 100:
        return None, None

    nreads = min(1000, len(read_segms)) // 2

    train, test = read_segms[:nreads], read_segms[nreads:nreads*2]
    jumps = Counter(get_cuts(train, A, B))
    mc_jump = jumps.most_common(1)[0]
    print(jumps)
    if mc_jump[1] < 5:
        return None, None
    mc_jump, orient = mc_jump[0]

    if orient == '>':
        hybrid = A[:mc_jump[0]] + B[mc_jump[1]:]
    else:
        hybrid = B[:mc_jump[0]] + A[mc_jump[1]:]

    train_n_improved = n_improved_pairs(train, A, B, hybrid)
    test_n_improved = n_improved_pairs(test, A, B, hybrid)

    stat, pval = proportions_ztest([train_n_improved, test_n_improved],
                                   [len(train), len(test)])
    print(mc_jump, orient)
    print(train_n_improved, test_n_improved, len(train), len(test))
    print(stat, pval)
    if np.isnan(pval) or pval < sign_lev:
        return None, None

    return hybrid, (mc_jump, orient)


def main():
    params = parse_args()
    sd_report, reads = read_input(params)
    putative_pairs = ident_putative_pairs(sd_report.df,
                                          sd_report.monomer_names_map,
                                          max_ident_diff=params.max_ident_diff,
                                          min_occ=500)
    print(putative_pairs)
    hybrids = {}
    for pair in putative_pairs:
        A = sd_report.monomers[sd_report.rev_monomer_names_map[pair[0]]]
        B = sd_report.monomers[sd_report.rev_monomer_names_map[pair[1]]]
        hybrid, jump = extract_hybrid(pair, sd_report.df, reads,
                                      params.max_ident_diff,
                                      A, B)
        if hybrid is not None:
            print(pair, hybrid, jump)
            hybrids[f'{pair[0]}_{pair[1]}|{jump[0][0]}_{jump[0][1]}|{jump[1]}'] = hybrid
    write_bio_seqs(params.outfile, hybrids)


if __name__ == "__main__":
    main()
